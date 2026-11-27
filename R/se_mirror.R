# =============================================================================
# SummarizedExperiment mirror: uns["se"]
#
# build_explore_h5ad(mirror = se) stores a whole SummarizedExperiment in the
# .h5ad, so that se_from_h5ad() rebuilds it: Every value at its path, with its
# class and attributes. Spec: inst/spec/exploreDE-h5ad-1.0.md, section 11.
#
# The value tree is walked depth-first. Every node (list, S4 object, atomic
# vector, environment, function, ...) is one row of uns$se$manifest, in the
# order of the walk, with its parent's id. Nodes are taken apart generically:
# Lists into their elements, S4 objects into their slots, and every object's
# attributes (names, dim, class, row.names, slots) into child rows. Only bare
# vectors (no attributes) hold data; each is placed by its shape:
#
#   genes x samples matrix          -> layers[<label>]
#   one value per gene / genes x k  -> varm[<label>] (data frame / matrix)
#   one value per sample / samples x k -> obsm[<label>]
#   anything else                   -> uns$se$values, nested like the R paths
#
# A bare vector identical() to one already stored (or to an exploreDE/1.0
# slot: counts, obs/var columns, DE columns, enrichment rows) is not stored
# again: Its manifest row points at that location. A node identical() to an
# earlier node is a "ref" to it. Lists of many bare vectors of one type (gene
# sets) become one long vector plus lengths.
# =============================================================================

.MIRROR_LONG_MIN <- 20L # a list of at least this many bare vectors is stored long

# An empty symbol (a missing argument: formals without a default, unforced
# frames) cannot be passed as an argument, so it is tested where it sits and
# passed around as this marker.
.MIRROR_MISSING <- structure(list(), class = "se_mirror_missing")
.mirror_elem <- function(lst, i) {
  if (identical(.subset2(lst, i), quote(expr = ))) .MIRROR_MISSING else .subset2(lst, i)
}
.mirror_binding <- function(e, nm) {
  tryCatch(
    {
      v <- get(nm, envir = e, inherits = FALSE)
      if (identical(v, quote(expr = ))) .MIRROR_MISSING else v
    },
    error = function(err) .MIRROR_MISSING
  )
}

# ---- small helpers ------------------------------------------------------------

.mirror_enc <- function(x) {
  x <- gsub("%", "%25", x, fixed = TRUE)
  x <- gsub("/", "%2F", x, fixed = TRUE)
  x <- gsub("|", "%7C", x, fixed = TRUE)
  sub("^([@$^\\[])", "%\\1", x) # element names never look like markers
}
.mirror_dec <- function(x) utils::URLdecode(sub("^%([@$^\\[])", "\\1", x))

# Path segment of a child (for uns nesting and labels).
.mirror_seg <- function(kind, key) {
  switch(
    kind,
    elem = if (grepl("^\\[[0-9]+\\]$", key)) key else .mirror_enc(key),
    attr = paste0("@", .mirror_enc(key)),
    binding = paste0("$", .mirror_enc(key)),
    paste0("^", kind)
  )
}

# Readable label for a layers/varm/obsm key: the named path, without the
# structural parts of Bioconductor containers.
.mirror_label <- function(segs) {
  s <- vapply(segs, .mirror_dec, "")
  s <- sub("^[@$^]", "", s)
  s <- s[!s %in% c("data", "listData", "metadata", "values")]
  if (!length(s)) "root" else paste(s, collapse = ".")
}

# R's human-readable path, for the manifest's path column.
.mirror_path <- function(parent_path, kind, key) {
  switch(
    kind,
    root = "se",
    elem = if (grepl("^\\[[0-9]+\\]$", key)) paste0(parent_path, key) else paste0(parent_path, "$", key),
    attr = paste0(parent_path, "@", key),
    binding = paste0(parent_path, "$", key),
    paste0(parent_path, " ", kind)
  )
}

.mirror_special_env <- function(e) {
  if (identical(e, globalenv())) {
    return("R_GlobalEnv")
  }
  if (identical(e, baseenv())) {
    return("base")
  }
  if (identical(e, emptyenv())) {
    return("R_EmptyEnv")
  }
  if (isNamespace(e)) {
    return(paste0("namespace:", getNamespaceName(e)))
  }
  nm <- environmentName(e)
  if (startsWith(nm, "package:")) {
    return(nm)
  }
  NULL
}

.mirror_env_from_name <- function(nm) {
  switch(
    nm,
    R_GlobalEnv = globalenv(),
    base = baseenv(),
    R_EmptyEnv = emptyenv(),
    if (startsWith(nm, "namespace:")) {
      asNamespace(sub("^namespace:", "", nm))
    } else if (startsWith(nm, "package:")) {
      as.environment(nm)
    } else {
      stop("se_from_h5ad: unknown environment ", nm)
    }
  )
}

# Exact text of a call: hexadecimal numbers, kept integers and NAs, but no
# quote() wrapper (which deparse(control = "exact") adds).
.mirror_deparse <- function(x) {
  paste(
    deparse(x, control = c("keepNA", "keepInteger", "niceNames", "hexNumeric"), width.cutoff = 500L),
    collapse = "\n"
  )
}
# Parsed, never evaluated: A file must not be able to run code when read. An
# expression() is rebuilt from its parsed arguments.
.mirror_parse <- function(txt, type) {
  tryCatch(
    {
      x <- str2lang(txt)
      if (type == "expression") {
        if (!is.call(x) || !identical(x[[1]], as.name("expression"))) {
          stop("not an expression() call")
        }
        as.expression(as.list(x)[-1])
      } else {
        x
      }
    },
    error = function(e) NULL
  )
}

# A builtin is stored as .Primitive("<name>"); only that form is accepted.
.mirror_builtin <- function(txt) {
  name <- sub('^\\.Primitive\\("(.*)"\\)$', "\\1", txt)
  if (identical(name, txt)) {
    stop("se_from_h5ad: unexpected builtin ", txt)
  }
  .Primitive(name)
}

# ---- encoder ----------------------------------------------------------------

# Encoder state: manifest rows, the slot contents, and the lookup tables for
# values already stored (pool), nodes already walked (nodes) and environments.
.mirror_state <- function(G, S, var_names, obs_names) {
  st <- new.env(parent = emptyenv())
  st$G <- G
  st$S <- S
  st$var_names <- var_names
  st$obs_names <- obs_names
  st$rows <- vector("list", 4096L)
  st$n <- 0L
  st$uns <- list() # list of list(segs, value)
  st$layers <- list()
  st$varm <- list()
  st$obsm <- list()
  st$varm_df <- list()
  st$obsm_df <- list()
  st$taken <- list(layers = character(), varm = character(), obsm = character())
  st$groups <- list(varm = list(), obsm = list()) # one key per data frame label
  st$pool <- list() # "typeof:length" -> list of list(loc, detail, value)
  st$nodes <- list() # "typeof:class:length" -> list of list(id, value)
  st$envs <- list() # list of list(env, id)
  st$env_queue <- list()
  st
}

.mirror_row <- function(st, parent, kind, key, path, x, how = "", location = "", detail = "") {
  st$n <- st$n + 1L
  if (st$n > length(st$rows)) {
    length(st$rows) <- 2L * length(st$rows)
  }
  st$rows[[st$n]] <- list(
    id = st$n,
    parent = parent,
    kind = kind,
    key = key,
    path = path,
    type = typeof(x),
    class = paste(class(x), collapse = "/"),
    s4 = isS4(x),
    how = how,
    location = location,
    detail = detail
  )
  st$n
}

.mirror_set <- function(st, id, ...) {
  v <- list(...)
  for (k in names(v)) {
    st$rows[[id]][[k]] <- v[[k]]
  }
}

.mirror_pool_key <- function(v) paste(typeof(v), length(v), sep = ":")

.mirror_pool_find <- function(st, v) {
  for (cand in st$pool[[.mirror_pool_key(v)]]) {
    if (identical(cand$value, v)) {
      return(cand)
    }
  }
  NULL
}

.mirror_pool_add <- function(st, v, loc, detail = "") {
  k <- .mirror_pool_key(v)
  st$pool[[k]] <- c(st$pool[[k]], list(list(loc = loc, detail = detail, value = v)))
}

.mirror_take <- function(st, slot, label) {
  key <- .mirror_unique_key(st$taken[[slot]], label)
  st$taken[[slot]] <- c(st$taken[[slot]], key)
  key
}
.mirror_group <- function(st, slot, label) {
  key <- st$groups[[slot]][[label]]
  if (is.null(key)) {
    key <- .mirror_take(st, slot, label)
    st$groups[[slot]][[label]] <- key
  }
  key
}

.mirror_unique_key <- function(existing, key) {
  if (!key %in% existing) {
    return(key)
  }
  i <- 2L
  while (paste0(key, "_", i) %in% existing) {
    i <- i + 1L
  }
  paste0(key, "_", i)
}

# HDF5 layers keep NaN but not NA_real_: Record where a double vector headed
# for a layer has NA_real_, so the reader can put them back.
.mirror_na_detail <- function(st, v, layer_key) {
  if (!is.double(v)) {
    return(character())
  }
  na_pos <- which(is.na(v) & !is.nan(v))
  if (!length(na_pos)) {
    return(character())
  }
  st$uns[[length(st$uns) + 1L]] <- list(segs = c("^na", .mirror_enc(layer_key)), value = na_pos)
  paste0("na=", .mirror_enc(layer_key))
}

# Where a bare vector goes, by its shape (spec section 11).
.mirror_slot_for <- function(v, dim, G, S) {
  is_mat <- length(dim) == 2L
  d1 <- if (is_mat) dim[1] else length(v)
  d2 <- if (is_mat) dim[2] else NA_integer_
  # in order of precedence; anything else goes to uns
  rules <- c(
    layer = is_mat & !is.character(v) & d1 == G & d2 == S,
    varm_matrix = is_mat & d1 == G & d2 != S,
    obsm_matrix = is_mat & d1 == S & d2 != G,
    varm_column = !is_mat & G != S & d1 == G,
    obsm_column = !is_mat & G != S & S > 1L & d1 == S
  )
  c(names(rules)[rules %in% TRUE], "uns")[1]
}

.mirror_as_matrix <- function(v, nr) {
  m <- matrix(v, nrow = nr)
  if (is.logical(m)) {
    storage.mode(m) <- "integer" # HDF5 matrices cannot hold logicals; "lgl" restores them
  }
  m
}

.mirror_store_layer <- function(st, v, segs) {
  key <- .mirror_take(st, "layers", .mirror_label(segs))
  st$layers[[key]] <- t(.mirror_as_matrix(v, st$G))
  list(loc = paste0("layers|", key), detail = .mirror_na_detail(st, v, key))
}

.mirror_store_matrix <- function(st, v, segs, slot, nr) {
  key <- .mirror_take(st, slot, .mirror_label(segs))
  st[[slot]][[key]] <- .mirror_as_matrix(v, nr)
  list(loc = paste0(slot, "|", key), detail = character())
}

.mirror_store_column <- function(st, v, segs, slot) {
  df_slot <- paste0(slot, "_df")
  grp <- .mirror_group(st, slot, .mirror_label(utils::head(segs, -1L)))
  col <- .mirror_unique_key(names(st[[df_slot]][[grp]]), .mirror_label(utils::tail(segs, 1L)))
  st[[df_slot]][[grp]][[col]] <- v
  list(loc = paste(slot, grp, col, sep = "|"), detail = character())
}

.mirror_store_uns <- function(st, v, segs) {
  val <- v
  detail <- character()
  # Long text with repeats (gene sets, go2gene) is stored as categorical.
  if (is.character(v) && length(v) >= 1000L && length(unique(v)) < 0.5 * length(v)) {
    val <- factor(v)
    detail <- "fct"
  }
  st$uns[[length(st$uns) + 1L]] <- list(segs = segs, value = val)
  list(loc = paste0("uns|", paste(segs, collapse = "/")), detail = detail)
}

.mirror_store <- function(st, v, slot, segs) {
  switch(
    slot,
    layer = .mirror_store_layer(st, v, segs),
    varm_matrix = .mirror_store_matrix(st, v, segs, "varm", st$G),
    obsm_matrix = .mirror_store_matrix(st, v, segs, "obsm", st$S),
    varm_column = .mirror_store_column(st, v, segs, "varm"),
    obsm_column = .mirror_store_column(st, v, segs, "obsm"),
    .mirror_store_uns(st, v, segs)
  )
}

# A value already stored: Point at it (a layer view still records its NAs).
.mirror_put_hit <- function(st, v, hit) {
  detail <- hit$detail
  if (startsWith(hit$loc, "layers|") && !grepl("na=", detail, fixed = TRUE)) {
    extra <- .mirror_na_detail(st, v, sub("^layers\\|", "", hit$loc))
    detail <- paste(c(if (nzchar(detail)) detail, extra), collapse = ";")
  }
  list(location = hit$loc, detail = detail)
}

# Store one bare vector (no attributes), placed by its shape; returns
# list(location, detail). `dim` is the dim attribute of the node it came from.
.mirror_put <- function(st, v, dim, segs, force_uns = FALSE) {
  hit <- .mirror_pool_find(st, v)
  if (!is.null(hit)) {
    return(.mirror_put_hit(st, v, hit))
  }
  slot <- if (force_uns) "uns" else .mirror_slot_for(v, dim, st$G, st$S)
  stored <- .mirror_store(st, v, slot, segs)
  detail <- paste(c(if (is.logical(v)) "lgl", stored$detail), collapse = ";")
  .mirror_pool_add(st, v, stored$loc, detail)
  list(location = stored$loc, detail = detail)
}

# Segments of a child: Attributes of an atomic node sit beside it (an HDF5
# dataset has no children), everything else inside its group.
.mirror_child_segs <- function(segs, x, kind, key) {
  seg <- .mirror_seg(kind, key)
  if (is.atomic(x) && kind == "attr" && length(segs)) {
    c(utils::head(segs, -1L), paste0(utils::tail(segs, 1L), seg))
  } else {
    c(segs, seg)
  }
}

.mirror_attrs <- function(st, x, id, path, segs) {
  at <- attributes(x)
  for (a in names(at)) {
    .mirror_node(st, at[[a]], id, "attr", a, path, .mirror_child_segs(segs, x, "attr", a))
  }
}

.mirror_is_long <- function(x) {
  if (length(x) < .MIRROR_LONG_MIN) {
    return(FALSE)
  }
  at <- attributes(x)
  if (length(setdiff(names(at), "names"))) {
    return(FALSE)
  }
  ty <- typeof(.subset2(x, 1L))
  if (!ty %in% c("logical", "integer", "double", "character")) {
    return(FALSE)
  }
  all(vapply(
    seq_along(x),
    function(i) {
      e <- .subset2(x, i)
      typeof(e) == ty && is.null(attributes(e))
    },
    logical(1)
  ))
}

.mirror_node <- function(st, x, parent, kind, key, parent_path, segs) {
  path <- .mirror_path(parent_path, kind, key)
  if (inherits(x, "se_mirror_missing")) {
    id <- .mirror_row(st, parent, kind, key, path, NULL, how = "missing")
    st$rows[[id]]$type <- "symbol"
    st$rows[[id]]$class <- "name"
    return(invisible(id))
  }
  id <- .mirror_row(st, parent, kind, key, path, x)
  # A node identical() to one already walked is a reference to it.
  bucket <- .mirror_bucket(x)
  if (!is.null(bucket)) {
    hit <- .mirror_find_node(st, bucket, x)
    if (!is.null(hit)) {
      .mirror_set(st, id, how = "ref", location = as.character(hit))
      return(invisible(id))
    }
  }
  .mirror_encoder(x)(list(st = st, id = id, path = path, segs = segs, kind = kind), x)
  if (!is.null(bucket)) {
    st$nodes[[bucket]] <- c(st$nodes[[bucket]], list(list(id = id, value = x)))
  }
  invisible(id)
}

.mirror_bucket <- function(x) {
  ty <- typeof(x)
  if (ty %in% c("list", "S4", "closure", "language") || (is.atomic(x) && !is.null(attributes(x)))) {
    paste(ty, paste(class(x), collapse = "/"), length(x), sep = ":")
  }
}

.mirror_find_node <- function(st, bucket, x) {
  for (cand in st$nodes[[bucket]]) {
    if (identical(cand$value, x)) {
      return(cand$id)
    }
  }
  NULL
}

.mirror_child <- function(ctx, value, kind, key, segs) {
  .mirror_node(ctx$st, value, ctx$id, kind, key, ctx$path, segs)
}

.mirror_encoder <- function(x) {
  switch(
    typeof(x),
    "NULL" = function(ctx, x) .mirror_set(ctx$st, ctx$id, how = "null"),
    environment = function(ctx, x) .mirror_env(ctx$st, x, ctx$id, ctx$path, ctx$segs),
    symbol = function(ctx, x) .mirror_set(ctx$st, ctx$id, how = "symbol", detail = as.character(x)),
    language = ,
    expression = .mirror_enc_language,
    closure = .mirror_enc_closure,
    builtin = ,
    special = function(ctx, x) .mirror_set(ctx$st, ctx$id, how = "builtin", detail = paste(deparse(x), collapse = "")),
    list = ,
    pairlist = .mirror_enc_list,
    S4 = function(ctx, x) {
      .mirror_set(ctx$st, ctx$id, how = "s4")
      .mirror_attrs(ctx$st, x, ctx$id, ctx$path, ctx$segs)
    },
    if (is.atomic(x)) .mirror_enc_atomic else .mirror_enc_unsupported
  )
}

.mirror_enc_language <- function(ctx, x) {
  bare <- x
  attributes(bare) <- NULL
  txt <- .mirror_deparse(bare)
  if (identical(.mirror_parse(txt, typeof(x)), bare)) {
    put <- .mirror_put(ctx$st, txt, NULL, c(ctx$segs, "^text"), force_uns = TRUE)
    .mirror_set(ctx$st, ctx$id, how = "language", location = put$location)
  } else {
    # Text does not give it back exactly (an embedded value, say): Store the
    # call's parts instead.
    .mirror_set(ctx$st, ctx$id, how = "call")
    parts <- as.list(bare)
    for (i in seq_along(parts)) {
      k <- paste0("[", i, "]")
      .mirror_child(ctx, .mirror_elem(parts, i), "elem", k, c(ctx$segs, k))
    }
  }
  .mirror_attrs(ctx$st, x, ctx$id, ctx$path, ctx$segs)
}

.mirror_enc_closure <- function(ctx, x) {
  .mirror_set(ctx$st, ctx$id, how = "closure")
  .mirror_child(ctx, as.list(formals(x)), "formals", "", c(ctx$segs, "^formals"))
  .mirror_child(ctx, body(x), "body", "", c(ctx$segs, "^body"))
  .mirror_child(ctx, environment(x), "env", "", c(ctx$segs, "^env"))
  .mirror_attrs(ctx$st, x, ctx$id, ctx$path, ctx$segs)
}

# An element's key: Its name, or [i] when it has none or shares it.
.mirror_elem_key <- function(nms, i) {
  nm <- nms[i]
  if (is.null(nm) || is.na(nm) || !nzchar(nm) || sum(nms == nm, na.rm = TRUE) > 1L) paste0("[", i, "]") else nm
}

.mirror_enc_list <- function(ctx, x) {
  if (typeof(x) == "pairlist") {
    x <- as.list(x)
    .mirror_set(ctx$st, ctx$id, detail = "pairlist")
  }
  xx <- unclass(x) # walk the list itself, not an S3 `[[` method (package_version)
  if (.mirror_is_long(xx)) {
    .mirror_set(ctx$st, ctx$id, how = "long")
    .mirror_child(ctx, unlist(xx, use.names = FALSE), "values", "", c(ctx$segs, "^values"))
    .mirror_child(ctx, lengths(xx), "lengths", "", c(ctx$segs, "^lengths"))
  } else {
    .mirror_set(ctx$st, ctx$id, how = "list")
    nms <- names(xx)
    for (i in seq_along(xx)) {
      k <- .mirror_elem_key(nms, i)
      .mirror_child(ctx, .mirror_elem(xx, i), "elem", k, c(ctx$segs, .mirror_seg("elem", k)))
    }
  }
  .mirror_attrs(ctx$st, x, ctx$id, ctx$path, ctx$segs)
}

.mirror_enc_atomic <- function(ctx, x) {
  bare <- x
  attributes(bare) <- NULL
  put <- .mirror_put(ctx$st, bare, attr(x, "dim"), ctx$segs, force_uns = ctx$kind %in% c("values", "lengths"))
  .mirror_set(ctx$st, ctx$id, how = "atomic", location = put$location, detail = put$detail)
  .mirror_attrs(ctx$st, x, ctx$id, ctx$path, ctx$segs)
}

.mirror_enc_unsupported <- function(ctx, x) {
  warning(sprintf(
    "build_explore_h5ad: %s (%s) cannot be stored; se_from_h5ad() gives NULL there",
    ctx$path,
    typeof(x)
  ))
  .mirror_set(ctx$st, ctx$id, how = "unsupported")
}

.mirror_env <- function(st, e, id, path, segs) {
  special <- .mirror_special_env(e)
  if (!is.null(special)) {
    .mirror_set(st, id, how = "envname", detail = special)
    return(invisible())
  }
  for (seen in st$envs) {
    if (identical(seen$env, e)) {
      .mirror_set(st, id, how = "envref", location = as.character(seen$id))
      return(invisible())
    }
  }
  st$envs[[length(st$envs) + 1L]] <- list(env = e, id = id)
  .mirror_set(st, id, how = "env")
  # Walked after the rest of the tree: An environment a function or formula
  # captured often holds copies of values that have a proper path elsewhere
  # (DESeq2's design formula holds the dds itself). Walking it last makes those
  # copies refs to the proper path.
  st$env_queue[[length(st$env_queue) + 1L]] <- list(env = e, id = id, path = path, segs = segs)
}

.mirror_env_contents <- function(st, item) {
  e <- item$env
  for (nm in ls(e, all.names = TRUE, sorted = FALSE)) {
    .mirror_node(
      st,
      .mirror_binding(e, nm),
      item$id,
      "binding",
      nm,
      item$path,
      c(item$segs, .mirror_seg("binding", nm))
    )
  }
  .mirror_node(st, parent.env(e), item$id, "parent", "", item$path, c(item$segs, "^parent"))
  .mirror_attrs(st, e, item$id, item$path, item$segs)
}

# Values the exploreDE/1.0 parts already hold, as the reader decodes them, so
# that a bare vector of the SE equal to one is stored once.
.mirror_seed_pool <- function(st, parts) {
  for (k in names(parts$layers)) {
    .mirror_pool_add(st, as.vector(t(parts$layers[[k]])), paste0("layers|", k))
  }
  .mirror_pool_add(st, rownames(parts$obs), "obs_names")
  .mirror_pool_add(st, rownames(parts$var), "var_names")
  .mirror_seed_frame(st, parts$obs, "obs")
  .mirror_seed_frame(st, parts$var, "var")
  for (k in names(parts$varm)) {
    .mirror_seed_frame(st, parts$varm[[k]], paste0("varm|", k), factors = FALSE)
  }
  enr <- parts$uns$enrichment %||% list()
  for (slot in names(enr)) {
    for (m in names(enr[[slot]])) {
      .mirror_seed_table(st, enr[[slot]][[m]], paste("enrichment", .mirror_enc(slot), .mirror_enc(m), sep = "/"))
    }
  }
}

.mirror_seed_frame <- function(st, df, where, factors = TRUE) {
  for (cl in colnames(df)) {
    v <- df[[cl]]
    if (factors && is.factor(v)) {
      .mirror_pool_add(st, as.character(v), paste0(where, "|", cl), "chr")
    } else {
      .mirror_pool_add(st, as.vector(v), paste0(where, "|", cl))
    }
  }
}

# An enrichment table: Each column of each ontology's rows.
.mirror_seed_table <- function(st, tab, where) {
  if (!is.data.frame(tab) || !"ontology" %in% colnames(tab)) {
    return(invisible())
  }
  for (o in unique(tab$ontology)) {
    rows <- tab$ontology == o
    for (cl in colnames(tab)) {
      loc <- paste("unsrows", where, .mirror_enc(cl), .mirror_enc(o), sep = "|")
      .mirror_pool_add(st, as.vector(tab[[cl]][rows]), loc)
    }
  }
}

# Nest list(segs, value) items into one list tree.
.mirror_nest <- function(items) {
  tree <- list()
  put <- function(tree, segs, value) {
    if (length(segs) == 1L) {
      tree[[segs]] <- value
      return(tree)
    }
    sub <- tree[[segs[1]]] %||% list()
    tree[[segs[1]]] <- put(sub, segs[-1], value)
    tree
  }
  for (it in items) {
    tree <- put(tree, it$segs, it$value)
  }
  tree
}

# Add the mirror of `se` to the exploreDE/1.0 parts.
.mirror_add <- function(parts, se) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("build_explore_h5ad: `mirror` must be a SummarizedExperiment")
  }
  G <- nrow(parts$var)
  S <- nrow(parts$obs)
  if (nrow(se) != G || ncol(se) != S) {
    stop(sprintf("build_explore_h5ad: `mirror` is %d x %d, the file %d x %d", nrow(se), ncol(se), G, S))
  }
  st <- .mirror_state(G, S, rownames(parts$var), rownames(parts$obs))
  st$taken <- list(layers = names(parts$layers), varm = names(parts$varm), obsm = names(parts$obsm))
  .mirror_seed_pool(st, parts)
  .mirror_node(st, se, 0L, "root", "", "", character())
  while (length(st$env_queue)) {
    item <- st$env_queue[[1]]
    st$env_queue[[1]] <- NULL
    .mirror_env_contents(st, item)
  }

  rows <- st$rows[seq_len(st$n)]
  col <- function(k, fn) fn(vapply(rows, function(r) r[[k]], fn(NA)))
  manifest <- data.frame(
    id = col("id", as.integer),
    parent = col("parent", as.integer),
    kind = col("kind", as.character),
    key = col("key", as.character),
    how = col("how", as.character),
    type = col("type", as.character),
    class = col("class", as.character),
    s4 = col("s4", as.logical),
    location = col("location", as.character),
    detail = col("detail", as.character),
    path = col("path", as.character),
    stringsAsFactors = FALSE
  )
  parts$layers <- c(parts$layers, st$layers)
  varm_df <- lapply(st$varm_df, function(cols) {
    d <- as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE, optional = TRUE)
    rownames(d) <- st$var_names
    d
  })
  obsm_df <- lapply(st$obsm_df, function(cols) {
    d <- as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE, optional = TRUE)
    rownames(d) <- st$obs_names
    d
  })
  parts$varm <- c(parts$varm, st$varm, varm_df)
  parts$obsm <- c(parts$obsm %||% list(), st$obsm, obsm_df)
  parts$uns$se <- list(
    schema = "exploreDE-se/1",
    manifest = manifest,
    values = .mirror_nest(st$uns)
  )
  parts
}

# ---- decoder ----------------------------------------------------------------

.mirror_drop_1d <- function(v) {
  if (length(dim(v)) == 1L) as.vector(v) else v
}

.mirror_uns_get <- function(uns, segs) {
  v <- uns
  for (s in segs) {
    v <- v[[s]]
    if (is.null(v)) stop("se_from_h5ad: missing uns entry ", paste(segs, collapse = "/"))
  }
  v
}

# The bare vector at a manifest location, as stored.
.mirror_fetch <- function(adata, loc, detail) {
  p <- strsplit(loc, "|", fixed = TRUE)[[1]]
  det <- if (nzchar(detail)) strsplit(detail, ";", fixed = TRUE)[[1]] else character()
  v <- switch(
    p[1],
    layers = as.vector(t(as.matrix(adata$layers[[p[2]]]))),
    varm = if (length(p) == 2L) as.vector(as.matrix(adata$varm[[p[2]]])) else adata$varm[[p[2]]][[p[3]]],
    obsm = if (length(p) == 2L) as.vector(as.matrix(adata$obsm[[p[2]]])) else adata$obsm[[p[2]]][[p[3]]],
    obs = adata$obs[[p[2]]],
    var = adata$var[[p[2]]],
    obs_names = adata$obs_names,
    var_names = adata$var_names,
    uns = .mirror_uns_get(adata$uns$se$values, strsplit(p[2], "/", fixed = TRUE)[[1]]),
    unsrows = {
      tab <- .mirror_uns_get(adata$uns, vapply(strsplit(p[2], "/", fixed = TRUE)[[1]], .mirror_dec, ""))
      cl <- .mirror_dec(p[3])
      tab[[cl]][tab$ontology == .mirror_dec(p[4])]
    },
    stop("se_from_h5ad: unknown location ", loc)
  )
  v <- .mirror_drop_1d(v)
  if (is.factor(v)) {
    v <- as.character(v)
  }
  attributes(v) <- NULL
  if ("lgl" %in% det) {
    v <- as.logical(v)
  }
  na <- grep("^na=", det, value = TRUE)
  if (length(na)) {
    pos <- .mirror_drop_1d(adata$uns$se$values[["^na"]][[sub("^na=", "", na[1])]])
    v[pos] <- NA_real_
  }
  v
}

#' Rebuild the SummarizedExperiment stored in an exploreDE `.h5ad`
#'
#' Files written by [build_explore_h5ad()] with `mirror = se` hold the whole
#' SummarizedExperiment in `uns["se"]` (spec section 11). This rebuilds it:
#' Every value at its original path, with its class and attributes, including
#' S4 objects such as `DESeqDataSet`, `DGEGLM`, `enrichResult`, and
#' `gseaResult` (their packages must be installed), functions, and the
#' environments they refer to.
#'
#' @param x Path to an `.h5ad` file, or an anndataR `AnnData` read from one.
#' @return The SummarizedExperiment.
#' @export
se_from_h5ad <- function(x) {
  adata <- if (is.character(x)) anndataR::read_h5ad(x, as = "InMemoryAnnData") else x
  m <- adata$uns$se$manifest
  if (is.null(m)) {
    stop(
      "se_from_h5ad: this .h5ad holds no SummarizedExperiment (no uns$se); ",
      "write it with build_explore_h5ad(mirror = se)"
    )
  }
  m <- as.data.frame(lapply(m, .mirror_drop_1d), stringsAsFactors = FALSE)
  m <- m[order(m$id), , drop = FALSE]
  # Decoder state: The manifest, each node's children, decoded values, and the
  # environments built so far (filled after the rest of the tree, as written).
  dc <- new.env(parent = emptyenv())
  dc$adata <- adata
  dc$m <- m
  dc$kids <- split(m$id, factor(m$parent, levels = c(0L, m$id)))
  dc$cache <- vector("list", nrow(m))
  dc$envs <- list()
  dc$fill <- integer()
  root <- .mirror_decode(dc, .mirror_children(dc, 0L))
  while (length(dc$fill)) {
    id <- dc$fill[1]
    dc$fill <- dc$fill[-1]
    .mirror_fill_env(dc, id)
  }
  root
}

.mirror_children <- function(dc, id, kind = NULL) {
  k <- dc$kids[[as.character(id)]]
  if (!is.null(kind)) {
    k <- k[dc$m$kind[k] %in% kind]
  }
  k
}

.mirror_attrs_of <- function(dc, id) {
  a <- .mirror_children(dc, id, "attr")
  stats::setNames(lapply(a, function(i) .mirror_decode(dc, i)), dc$m$key[a])
}

.mirror_set_attrs <- function(x, at) {
  if (length(at)) {
    attributes(x) <- at
  }
  x
}

.mirror_decode <- function(dc, id) {
  hit <- dc$cache[[id]]
  if (!is.null(hit)) {
    return(hit[[1]])
  }
  r <- dc$m[id, ]
  val <- .mirror_decoder(r$how, r$path)(dc, id, r)
  dc$cache[[id]] <- list(val)
  val
}

.mirror_decoder <- function(how, path) {
  switch(
    how,
    ref = function(dc, id, r) .mirror_decode(dc, as.integer(r$location)),
    null = ,
    unsupported = function(dc, id, r) NULL,
    missing = function(dc, id, r) .MIRROR_MISSING,
    symbol = function(dc, id, r) as.name(r$detail),
    builtin = function(dc, id, r) .mirror_builtin(r$detail),
    envname = function(dc, id, r) .mirror_env_from_name(r$detail),
    envref = function(dc, id, r) dc$envs[[r$location]],
    env = .mirror_dec_env,
    language = function(dc, id, r) {
      .mirror_set_attrs(.mirror_parse(.mirror_fetch(dc$adata, r$location, ""), r$type), .mirror_attrs_of(dc, id))
    },
    call = .mirror_dec_call,
    closure = .mirror_dec_closure,
    long = .mirror_dec_long,
    list = .mirror_dec_list,
    s4 = .mirror_dec_s4,
    atomic = .mirror_dec_atomic,
    stop("se_from_h5ad: unknown manifest entry '", how, "' at ", path)
  )
}

.mirror_dec_env <- function(dc, id, r) {
  e <- new.env(parent = emptyenv())
  dc$envs[[as.character(id)]] <- e
  dc$fill <- c(dc$fill, id)
  e
}

.mirror_fill_env <- function(dc, id) {
  e <- dc$envs[[as.character(id)]]
  for (b in rev(.mirror_children(dc, id, "binding"))) {
    v <- .mirror_decode(dc, b)
    if (inherits(v, "se_mirror_missing")) {
      eval(call("<-", as.name(dc$m$key[b]), quote(expr = )), e)
    } else {
      assign(dc$m$key[b], v, envir = e)
    }
  }
  p <- .mirror_children(dc, id, "parent")
  if (length(p)) {
    parent.env(e) <- .mirror_decode(dc, p)
  }
  at <- .mirror_attrs_of(dc, id)
  if (length(at)) {
    attributes(e) <- at
  }
}

# Decoded elements, with empty arguments put back in place.
.mirror_decode_elems <- function(dc, ids) {
  x <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    v <- .mirror_decode(dc, ids[i])
    if (inherits(v, "se_mirror_missing")) {
      x[i] <- list(quote(expr = ))
    } else if (!is.null(v)) {
      x[[i]] <- v
    }
  }
  x
}

.mirror_dec_call <- function(dc, id, r) {
  parts <- .mirror_decode_elems(dc, .mirror_children(dc, id, "elem"))
  x <- if (r$type == "expression") as.expression(parts) else as.call(parts)
  .mirror_set_attrs(x, .mirror_attrs_of(dc, id))
}

.mirror_dec_closure <- function(dc, id, r) {
  fm <- .mirror_decode(dc, .mirror_children(dc, id, "formals"))
  bd <- .mirror_decode(dc, .mirror_children(dc, id, "body"))
  f <- eval(call("function", as.pairlist(fm), bd)) # builds the function; nothing in it runs
  environment(f) <- .mirror_decode(dc, .mirror_children(dc, id, "env"))
  .mirror_set_attrs(f, .mirror_attrs_of(dc, id))
}

.mirror_dec_long <- function(dc, id, r) {
  vals <- .mirror_decode(dc, .mirror_children(dc, id, "values"))
  lens <- .mirror_decode(dc, .mirror_children(dc, id, "lengths"))
  x <- unname(split(vals, factor(rep.int(seq_along(lens), lens), levels = seq_along(lens))))
  .mirror_set_attrs(x, .mirror_attrs_of(dc, id))
}

.mirror_dec_list <- function(dc, id, r) {
  x <- .mirror_decode_elems(dc, .mirror_children(dc, id, "elem"))
  if (r$detail == "pairlist") {
    x <- as.pairlist(x)
  }
  x <- .mirror_set_attrs(x, .mirror_attrs_of(dc, id))
  if (r$s4) asS4(x) else x
}

# The prototype of an S4 class, from its package; NULL when not available.
.mirror_s4_prototype <- function(cls) {
  tryCatch(
    {
      pkg <- attr(cls, "package")
      if (!is.null(pkg)) {
        requireNamespace(pkg, quietly = TRUE)
      }
      methods::getClassDef(cls[1], package = pkg %||% "")@prototype
    },
    error = function(e) NULL
  )
}

.mirror_dec_s4 <- function(dc, id, r) {
  at <- .mirror_attrs_of(dc, id)
  proto <- .mirror_s4_prototype(at$class)
  if (is.null(proto) || typeof(proto) != "S4") {
    warning(sprintf(
      "se_from_h5ad: class %s is not available; %s comes back as a list of its slots",
      at$class[1],
      r$path
    ))
    return(at)
  }
  attributes(proto) <- at
  proto
}

.mirror_dec_atomic <- function(dc, id, r) {
  v <- .mirror_fetch(dc$adata, r$location, r$detail)
  storage.mode(v) <- r$type
  x <- .mirror_set_attrs(v, .mirror_attrs_of(dc, id))
  if (r$s4) asS4(x) else x
}

# ---- comparison ---------------------------------------------------------------

# TRUE when `a` and `b` are identical(), except that environments other than
# namespaces and the global, base, and empty environments are compared by
# their contents (no copy, not even a qs2 round trip, keeps an environment's
# address). Otherwise the path of the first difference.
.se_mirror_compare <- function(a, b, path = "se", seen = new.env()) {
  if (identical(a, b)) {
    return(TRUE)
  }
  if (typeof(a) != typeof(b)) {
    return(sprintf("%s: typeof %s vs %s", path, typeof(a), typeof(b)))
  }
  cmp <- switch(
    typeof(a),
    environment = .se_mirror_cmp_env,
    closure = .se_mirror_cmp_closure,
    list = ,
    pairlist = ,
    expression = .se_mirror_cmp_list,
    S4 = .se_mirror_cmp_attrs,
    language = .se_mirror_cmp_bare,
    if (is.atomic(a)) .se_mirror_cmp_bare
  )
  if (is.null(cmp)) {
    return(sprintf("%s: %s differs", path, typeof(a)))
  }
  cmp(a, b, path, seen)
}

.se_mirror_cmp_attrs <- function(a, b, path, seen) {
  aa <- attributes(a)
  ab <- attributes(b)
  if (!setequal(names(aa), names(ab))) {
    return(sprintf("%s: attributes %s vs %s", path, paste(names(aa), collapse = ","), paste(names(ab), collapse = ",")))
  }
  for (n in names(aa)) {
    r <- .se_mirror_compare(aa[[n]], ab[[n]], paste0(path, "@", n), seen)
    if (!isTRUE(r)) {
      return(r)
    }
  }
  TRUE
}

.se_mirror_cmp_env <- function(a, b, path, seen) {
  sa <- .mirror_special_env(a)
  sb <- .mirror_special_env(b)
  if (!is.null(sa) || !is.null(sb)) {
    return(if (identical(sa, sb)) TRUE else sprintf("%s: environment %s vs %s", path, sa %||% "local", sb %||% "local"))
  }
  key <- paste(format(a), format(b))
  if (!is.null(seen[[key]])) {
    return(TRUE)
  }
  seen[[key]] <- TRUE
  r <- .se_mirror_cmp_bindings(a, b, path, seen)
  if (!isTRUE(r)) {
    return(r)
  }
  r <- .se_mirror_compare(parent.env(a), parent.env(b), paste0(path, " parent"), seen)
  if (!isTRUE(r)) {
    return(r)
  }
  .se_mirror_cmp_attrs(a, b, path, seen)
}

.se_mirror_cmp_bindings <- function(a, b, path, seen) {
  na <- sort(ls(a, all.names = TRUE))
  if (!identical(na, sort(ls(b, all.names = TRUE)))) {
    return(sprintf("%s: environment bindings differ", path))
  }
  for (n in na) {
    r <- .se_mirror_compare(.mirror_binding(a, n), .mirror_binding(b, n), paste0(path, "$", n), seen)
    if (!isTRUE(r)) {
      return(r)
    }
  }
  TRUE
}

.se_mirror_cmp_closure <- function(a, b, path, seen) {
  r <- .se_mirror_compare(as.list(formals(a)), as.list(formals(b)), paste0(path, " formals"), seen)
  if (!isTRUE(r)) {
    return(r)
  }
  if (!identical(body(a), body(b))) {
    return(sprintf("%s: function bodies differ", path))
  }
  r <- .se_mirror_compare(environment(a), environment(b), paste0(path, " env"), seen)
  if (!isTRUE(r)) {
    return(r)
  }
  .se_mirror_cmp_attrs(a, b, path, seen)
}

.se_mirror_cmp_list <- function(a, b, path, seen) {
  if (length(a) != length(b)) {
    return(sprintf("%s: length %d vs %d", path, length(a), length(b)))
  }
  ua <- if (typeof(a) == "list") unclass(a) else as.list(a)
  ub <- if (typeof(b) == "list") unclass(b) else as.list(b)
  for (i in seq_along(ua)) {
    r <- .se_mirror_compare(.mirror_elem(ua, i), .mirror_elem(ub, i), paste0(path, "[[", i, "]]"), seen)
    if (!isTRUE(r)) {
      return(r)
    }
  }
  .se_mirror_cmp_attrs(a, b, path, seen)
}

.se_mirror_cmp_bare <- function(a, b, path, seen) {
  ba <- a
  bb <- b
  attributes(ba) <- NULL
  attributes(bb) <- NULL
  if (!identical(ba, bb)) {
    return(sprintf("%s: values differ", path))
  }
  .se_mirror_cmp_attrs(a, b, path, seen)
}
