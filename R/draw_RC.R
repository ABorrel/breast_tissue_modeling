#!/usr/bin/env Rscript
library(ggplot2)




data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
p_response = args[1]
hormone = args[2]

#p_response = "/mnt/c/Users/AlexandreBorrel/research/SSI/mixture_breast_cycle/RESULTS/H295R/Responce_curve/ESTRADIOL/50-28-2"
#hormone = "Estradiol"


d_resp = read.csv(p_response, sep = "\t")
d_resp = as.data.frame(d_resp)

chemical = d_resp$chnm[1]
eff_pot = d_resp$Efficacy.potency[1]
fold_change = d_resp$CR.results[1]
  
d_resp$Concentration = log10(d_resp$Concentration)
d_resp$Response = log10(d_resp$Response)
d_resp$sd = sd(d_resp$Response)


ggplot(data = d_resp, aes(x = Concentration, y = Response)) +
  geom_point() +
  stat_summary(fun=mean, geom="line", colour="red", size=0.5, shape=4, linetype = "dashed") +
  labs(x = "Log Concentration Chemical (µM)", y = paste("Log Concentration", hormone, "(µM)", sep = " "), 
       title=paste(chemical, " (", basename(p_response) ,")", "\nEfficacy/potency: ", eff_pot, "\nFold change (Kaumas 2012): ", fold_change, sep=""))
  

ggsave(paste(p_response, ".png", sep = ""),  width = 5, height = 5, dpi = 300, bg="transparent")





