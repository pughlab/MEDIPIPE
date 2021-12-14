
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(MEDIPS)

IN = args[1];
OUT = args[2];

#print(args[1])
a <- read.table(IN, header = T);
write.table(a, file = OUT)
