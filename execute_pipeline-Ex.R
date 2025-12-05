#!/usr/bin/env Rscript
library(rmarkdown)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (report number)", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = ""
}

# Call this script with run_templates.sh and pass it your arguments. The args below
# are just an example. 
user.report.num=args[1]
user.suffix=args[2]


render_report = function(suffix, template.filename){
        config.file=paste0("config",  suffix, ".yaml")
        # template.filename=paste0("scRNA_template_PART", report.num, ".Rmd")
        # output.filename=paste0("report_part", report.num, suffix, "_output.html")
	output.name <- gsub(".Rmd", ".pdf", template.filename)

        rmarkdown::render(template.filename,
                          params = list( config.args = config.file), 
			  envir = new.env(parent = globalenv()), 
                          output_file = output.name )

}


render_report(user.suffix, user.report.num)
