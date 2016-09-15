options(echo=TRUE)

list.of.packages <- c("devtools","argparse","yaml","zoo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos="https://mirrors.ebi.ac.uk/CRAN/")

library(devtools)
install_github("asifrim/mrmosaic",ref="master",subdir="mrmosaic")
print("Mr.Mosaic R Package Installed...")
