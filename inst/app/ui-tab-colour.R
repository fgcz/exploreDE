tabItem(
  tabName = "colourTab",
  fluidPage(
    fluidRow(
      column(
        width = 12,
        box(
          title = "Group Colours & Labels",
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
          helpText(
            HTML(
              "You can also give any group a <b>plot label</b> in the box under its colour picker. Leave it empty to keep the original name.
            Labels accept a few HTML tags: <code>&lt;i&gt;</code> italics, <code>&lt;b&gt;</code> bold, <code>&lt;sup&gt;</code> superscript,
            <code>&lt;sub&gt;</code> subscript and <code>&lt;br&gt;</code> line break, plus entities like <code>&amp;Delta;</code> or <code>&amp;alpha;</code>
            for Greek letters. For example, <code>&lt;i&gt;Cxcr4&lt;/i&gt;&lt;sup&gt;fl/fl&lt;/sup&gt; (Treated)</code> gives <i>Cxcr4</i><sup>fl/fl</sup> (Treated).
            Labels are saved in the colour template too.
            Note: the plots shown in the app may display italics as upright text, but the downloaded figures will have them."
            )
          ),
          splitLayout(
            downloadButton(outputId = "dlColourTemplate", label = paste("Download Colour Template")),
            fileInput(
              inputId = "importColourFile",
              accept = ".xlsx",
              width = "85%",
              label = NULL,
              placeholder = "Upload Colour Template"
            )
          ),
          uiOutput("colourPaletteUI"),
          uiOutput("groupLabelWarningsUI"),
          uiOutput(outputId = "colourPickerUI", inline = T)
        )
      )
    )
  )
)
