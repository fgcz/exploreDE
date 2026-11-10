# Copied by ui-tests.yaml into the installed app directory as .Rprofile.
# shinytest2 starts the app's R process in that directory, so this profiles
# it: a sample every 0.5 s to EXPLOREDE_UI_ARTIFACTS/Rprof-<pid>.out, of
# elapsed time on Linux, so that a process blocked outside R is sampled too,
# with samples taken during garbage collection marked <GC>.
utils::Rprof(
  file.path(Sys.getenv("EXPLOREDE_UI_ARTIFACTS"), paste0("Rprof-", Sys.getpid(), ".out")),
  interval = 0.5,
  gc.profiling = TRUE,
  event = if (Sys.info()[["sysname"]] == "Linux") "elapsed" else "cpu"
)
