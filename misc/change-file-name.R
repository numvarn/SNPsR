# Config values
setwd("/Volumes/Sirikanlaya/case-control-replicated/0.0")


filez <- list.files("./")
i <- 1
for (name in filez) {
     if (i < 10) {
          prefix <- "000"
     }else if (i < 100) {
          prefix <- "00"
     }else if (i <= 1000) {
          prefix <- ""
     } 
     file.rename(name, paste(prefix, i, ".csv", sep = ""))
     cat(sprintf("%s %s\n", i, name))
     i <- i + 1
}