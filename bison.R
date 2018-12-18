library(ggplot2)

setwd(dir="~/work/bison/")
bison_data<-read.table("results.txt", header=T)
bison_diff_demes<-read.table("results_diff_demes.txt", header=T)

ggplot(bison_data,aes(class,het))+geom_boxplot() 

ggplot(bison_diff_demes,aes(class,het))+geom_boxplot() 


#system("conda activate my-momi-env")
# haplo = system("python bison_sim.py", intern=T)
# cat(paste0(haplo, "\n"))
