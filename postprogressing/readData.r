file <- read.csv("postprogressing/cell_data.csv",header = FALSE)
library(pracma)
file2 <- read.table("/mnt/e/mywork/programDesign/cppProject/riemann/postprogressing/data")
trapz(file2$V2,file2$V3)
