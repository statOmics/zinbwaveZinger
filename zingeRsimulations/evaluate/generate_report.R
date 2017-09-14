library(rmarkdown)

generate_report <- function(dds_list, show_code = FALSE, output_format = NULL,
                            output_file = NULL, output_dir = "./", ...){
  ## This function was written by Nicholas Hamilton and obtained from 
  ## http://stackoverflow.com/questions/37097535/generate-report-in-r
  
  ## Give the path to the template file
  theFile <- "simulation_comparison_template_list.Rmd"
  
  ## Process the arguments
  args <- list(...)
  args$input <- theFile
  args$output_dir <- output_dir
  args$output_format <- output_format
  args$output_file <- output_file
  
  ## Render the report
  outputFileName <- do.call('render', args = args)
  invisible(outputFileName)
}
