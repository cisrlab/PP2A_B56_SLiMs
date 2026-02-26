#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

#' Run itcpredictapp
#'
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
#' @import shiny
itcpredictApp <- function(...) {

  # Define UI for application that draws a histogram
  ui <- fluidPage(

    # Application title
    titlePanel("Predict PP2A-B56 SLiMs Kd values"),

    # Sidebar with a textArea input for sequences
    sidebarLayout(
        sidebarPanel(
            textAreaInput(inputId = "sequences", label="Sequences", value="<Sequences>"),
          actionButton(inputId = "submit", label="Submit"),
          actionButton(inputId = "example", label="Example"),
          actionButton(inputId = "clear", label="Clear"),
          textOutput("status"),
          textOutput("have_meme")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           downloadButton("download", "Download .tsv"),
           tableOutput("table")
           #plotOutput("distPlot")
        )
    )
)

exampleClicked<-function(...) {

  example_seqs <- c(
      "WTSFFSG.CSPIEEEAH",
      "WTSFFSG.CSPIEDEAH",
      "WTSFFSG.CSPIE(Sp)EAH",
      "WTSFFSG.CSPIEEAAH",
      "WTSFFSG.CSPLEEEAH",
      "WTSFFSG.CSPVEEEAH",
      "WTSFFSG.ASPIEEEAH",
      "WTSFFSG.CSPIEEDAH",
      "WTSFFSG.VSPIEEEAH",
      "WTSFFSG.ISPIEEEAH",
      "WTSFFSG.YSPIEEEAH",
      "WTSFFSG.WSPIEEEAH",
      "WTSFFSG.CSPIEEEEH",
      "WTSFFSG.CSPIEEEAD",
      "WTSFFSG.CSPIEEEAE",
      "WTSFFSG.CSPIEEEDH",
      "WTSFFSG.C(pS)PIEEEAH",
      "WTSFFSG.MSPIEEEAH",
      "WTSFFSG.FSPIEEEAH",
      "WTSFFSG.LSPIEEEAH",
      "WTSFFSG.LSPAEEEAH",
      "WTSFFSG.LSPIEAEAH",
      "WTSFFSG.LSPYEEEAH",
      "WTSFFSG.LSPIE(Sp)EAH",
      "WTSFFSG.LSPFEEEAH",
      "WTSFFSG.LSPMEEEAH",
      "WTSFFSG.LSPCEEEAH",
      "WTSFFSG.LSPWEEEAH",
      "WTSFFSG.LSPLEEEAH",
      "WTSFFSG.LSPIEEAAH",
      "WTSFFSG.LSPVEEEAH",
      "LSTLREQSSQS",
      "L(pS)TLREQSSQS",
      "LS(pT)LREQSSQS",
      "L(pS)(pT)LREQSSQS",
      "LSTLREQS(pS)QS",
      "LSTLREQSSQ(pS)",
      "LSTIDESGSIL",
      "L(pS)TIDESGSIL",
      "LS(pT)IDESGSIL",
      "LSTIDE(pS)GSIL",
      "L(pS)TIDE(pS)GSIL",
      "MQDIPEETES",
      "MQDIPEETE(pS)",
      "LRQSP.MQTIQENKPAT",
      "AQTAQENK",
      "LQTIQENK",
      "MLTPINEEA",
      "APPVQEDDE",
      "ASTAPEEEG",
      "LSIKK.LSPIIEDSREA",
      "GISGY.LPTLNEDEEWK",
      "VSTQE.LYSIPEDQEPE",
      "PGPEI.MRTIPEEELTD",
      "CSPIEEEAH",
      "LEPIEEEAH",
      "LEPIEEEPE",
      "MQAISE",
      "MPPIHE",
      "LESVAEEHE",
      "LSTIDESGS",
      "LSTIDEEGS",
      "LSTIDEEGE",
      "LSTIDDEGE",
      "LSTID(Tp)EGE",
      "L(pS)PIIEDDREADH",
      "LEPIIEDE",
      "K.L(pS)PIIEDE",
      "K.L(pS)PIIED(pS)",
      "K.L(pS)PIIED",
      "KKPL.L(pS)PIPELPE"
  )

  example_seqs_str <- paste(example_seqs, collapse="\n", sep="\n")
  updateTextAreaInput(inputId = "sequences", value = example_seqs_str)
}

clearClicked<-function(...) {
  updateTextAreaInput(inputId = "sequences", value ="<Sequences>")
  RV$data <- data.frame()
}

RV <- reactiveValues(data = data.frame(), status="Waiting...")

parseSequences<-function(sequences_str) {
  print(sequences_str)
  ans <- strsplit(sequences_str, '\n')[[1]]
  print(ans)
  return(ans)
}

submitClicked<-function(input,...) {
  print(names(input))
  message("submit clicked")
  RV$status <- "Processing..."
  sequences <- parseSequences(input$sequences)
  predictions <- PP2A.B56.SLiMs::getITC_2024_06_27_cv(sequences, enforce_6E = TRUE)

  RV$data <- data.frame(
      Sequence = sequences,
      ITC = predictions
  )

}


# Define server logic required to allow submission
server <- function(input, output) {
    observeEvent(input$example, {exampleClicked()})
    observeEvent(input$submit, {submitClicked(input)})
    observeEvent(input$clear, {clearClicked()})
    output$table <- renderTable(RV$data)
    output$status <- renderText(RV$status)
    output$download <- downloadHandler(
        filename = function() {
            "itc.tsv"
        },
        content = function(file) {
            vroom::vroom_write(RV$data, file)
        }
    )
    output$have_meme <- renderText(paste0("have fimo (meme):", itcpredictr:::have_ext_fimo()))

}

# Run the application
shinyApp(ui = ui, server = server)

}

itcpredictApp()


