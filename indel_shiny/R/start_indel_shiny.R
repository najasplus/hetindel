#' @export start_indel_shiny
#'
#' @return shiny application object
#'
#' @example \dontrun {start_indel_shiny()}
#'
#' @import shiny
#'
start_indel_shiny <-
function(){
	shinyApp(ui = shinyAppUi, server = shinyAppServer)
}