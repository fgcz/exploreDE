# --- telemetry_module.R ---
library(shiny)
library(shiny.telemetry)
library(RSQLite)
library(httr)
library(jsonlite)

# --- 1. JAVASCRIPT CODE ---
JS_TELEMETRY_CODE <- "
$(document).ready(function(){
  $.getJSON('https://api.ipify.org?format=json', function(data){
    Shiny.setInputValue('public_ip_js', data.ip);
  }).fail(function() {
    Shiny.setInputValue('public_ip_js', 'Block/Fail');
  });
});
$(document).on('click', function(e) {
  Shiny.setInputValue('page_click', {
    x: e.pageX, y: e.pageY, elem: e.target.id || e.target.tagName 
  }, {priority: 'event'});
});
"

# --- 2. UI HELPER ---
use_telemetry_js <- function() {
  tagList(
    use_telemetry(), 
    tags$head(tags$script(HTML(JS_TELEMETRY_CODE)))
  )
}

# --- 3. GLOBAL SETUP (With Concurrency Fix) ---
setup_telemetry_global <- function(app_name, data_dir = "/srv/shiny-server/data") {
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  # Single Persistent File
  db_filename <- paste0(app_name, ".sqlite")
  db_path <- file.path(data_dir, db_filename)
  
  message(paste("Telemetry: Logging to persistent file:", db_path))

  # --- CONCURRENCY: Enable WAL Mode ---
  # We connect briefly to set the PRAGMAs that allow concurrent writing.
  # 'journal_mode=WAL' is persistent (stays with the file).
  # 'busy_timeout' makes it wait 5000ms for a lock instead of crashing instantly.
  tryCatch({
    con <- dbConnect(RSQLite::SQLite(), db_path)
    dbExecute(con, "PRAGMA journal_mode=WAL;")
    dbExecute(con, "PRAGMA busy_timeout=5000;") 
    dbExecute(con, "PRAGMA synchronous=NORMAL;") # Faster writes, slightly less safe on power loss
    dbDisconnect(con)
  }, error = function(e) {
    message("Telemetry Warning: Could not set WAL mode. Concurrency might be limited.")
  })
  
  # Initialize Storage
  data_storage <- DataStorageSQLite$new(db_path = db_path)
  telemetry_obj <- Telemetry$new(app_name = app_name, data_storage = data_storage)
  
  return(telemetry_obj)
}

# --- 4. SERVER LOGIC ---
start_telemetry_server <- function(input, session, telemetry_obj) {
  
  get_user_id <- reactive({
    sp_user <- Sys.getenv("SHINYPROXY_USERNAME")
    if (sp_user != "") sp_user else session$token
  })
  
  get_location_info <- reactive({
    req(input$public_ip_js)
    user_ip <- input$public_ip_js
    if (user_ip == "Block/Fail") return("AdBlock/Firewall prevented IP lookup")
    
    tryCatch({
      url <- paste0("https://ipinfo.io/", user_ip, "/geo")
      data <- GET(url)
      if (status_code(data) != 200) return(paste("API Error:", status_code(data)))
      jsondata <- fromJSON(rawToChar(data$content))
      paste("City:", jsondata$city, "| Region:", jsondata$region, "| Country:", jsondata$country, "| Coords:", jsondata$loc)
    }, error = function(e) paste("Location Error:", e$message))
  })
  
  get_project_info <- reactive({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query$data) && query$data != "") return(query$data) else return("Unknown")
  })

  observe({
    req(input$public_ip_js)
    user <- get_user_id()
    full_loc <- paste0("IP: ", input$public_ip_js, " -> ", get_location_info())
    
    telemetry_obj$start_session(track_inputs = TRUE, track_values = TRUE, username = user)
    
    telemetry_obj$log_custom_event("session_details", list(
      project_info = get_project_info(),
      user_location = full_loc,
      user_agent   = session$request$HTTP_USER_AGENT
    ))
  }) |> bindEvent(input$public_ip_js, once = TRUE)
  
  observeEvent(input$page_click, {
    telemetry_obj$log_custom_event("click_coordinates", list(
      x = input$page_click$x, y = input$page_click$y, element = input$page_click$elem
    ))
  })
}
