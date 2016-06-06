options(echo=TRUE)

list.of.packages <- c("devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(devtools)
install_github("asifrim/mrmosaic",ref="development",subdir="mrmosaic")
print("Mr.Mosaic R Package Installed...")
