args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=2){
  stop("Two arguments (inputDir, sampleName) must be supplied.",call.=FALSE)
}

dir=as.character(args[1])
sample=as.character(args[2])

#Check to see if package(s) are installed, install if not and then load
  CheckInstallPackages <- function(pkgs) {                     #pkgs is a vector of strings with length >= 1
    x <- lapply(pkgs, function(pkg){                           #For each pkg in pkgs (attempt to load each package one at a time):
      if(!do.call("require", list(pkg))) {                     #  Load the package if available,
        try(install.packages(pkg, lib=.Library, repos="http://cran.us.r-project.org")) #    Silently attempt to install into the default library
        tryCatch(do.call("library", list(pkg)),                #    Now attempt to load the package, catch error if it wasn't installed
          error = function(err) {                              #    Catch if we're unable to install into the default library
            if(!interactive()) {                               #      If non-interactive, install into this user's personal library
              personalLibPath <- Sys.getenv("R_LIBS_USER")     #        Get the path to this user's personal library
              if(is.na(match(personalLibPath, .libPaths()))) { #        If the personal library is not in the list of libraries
                dir.create(personalLibPath, recursive = TRUE)  #          Then create the personal library
                .libPaths(personalLibPath)                     #          And add the personal library to the list of libraries
              }
              install.packages(pkg, lib=personalLibPath,       #        Attempt to install the package into the personal library
                              repos="http://cran.us.r-project.org") #          if this fails, raise the error back to the report
              do.call("library", list(pkg))                    #        Finally, attempt to load the package
            }
          }
        )
      }
    })
  }


##if(!require("fitdistrplus")){
  ##install.packages('fitdistrplus',repos='http://cran.us.r-project.org',lib=dir)
##}

CheckInstallPackages("fitdistrplus")
##library(fitdistrplus)


mut_type=c("A_C","A_G","A_T","C_A","C_G","C_T")
setwd(dir)

dupchk = list()
idx = 0

for(k in 1:6){
  errRate=c()
  
  inputFile=paste(sample,"_errorStatus_",mut_type[k],sep="")
  a=read.table(paste0(inputFile,'.txt'), header=T)
  ratio=a$ExpectedBE/a$TotalDepth
  x=sort(ratio)
  x_over_zero=x[x<quantile(x,0.99)]  
  y_over_zero=mledist(x_over_zero,"exp")
    
  mean=1/y_over_zero$estimate
  idx = idx + 1
  dupchk[[idx]] = paste(mut_type[k],mean,sep="\t")
}

header = paste("MutType","estimatedMean",sep="\t")
outputFileName=paste(sample,"_totalErrorStatus.txt",sep="")
write(header,outputFileName)
lapply(unique(dupchk), write, outputFileName, append=TRUE, ncolumns=2)
