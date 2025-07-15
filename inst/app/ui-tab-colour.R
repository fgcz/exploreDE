tabItem(
  tabName = "colourTab",
  fluidPage(
    fluidRow(
      column(
        width = 9,
        box(
          title = "Group Colours",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          helpText(
            "Now you can more easily set custom colours without having to choose them each time you load the app.
            Download the colour template Excel file and change the colours for each factor level in the second column.
            You can use hex codes (e.g., #750912, #59a7cf) or ggplot colour names (e.g., dodgerblue, plum4, etc.).
            Then, upload that file and watch the colour palette and your figures automagically update to the new colours!
            Don't change the values in the first column, or it won't work!
            Hopefully this is a much easier way to set colours once."
          ),
          splitLayout(
            downloadButton(outputId = "dlColourTemplate", label = paste("Download Colour Template")),
            fileInput(inputId = "importColourFile", accept = ".xlsx", width = "85%", label = NULL, placeholder = "Upload Colour Template")
          ),
          uiOutput("colourPaletteUI"),
          uiOutput(outputId = "colourPickerUI", inline = T)
        )
      )
    )
  )
)