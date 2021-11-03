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
pval = as.double(args[3])

#p_response = "/mnt/c/Users/AlexandreBorrel/research/SSI/mixture_breast_cycle/RESULTS/H295R/Responce_curve/Estradiol/borderline_active/629-76-5"
#hormone = "Estradiol"
#pval = 0.05

d_resp = read.csv(p_response, sep = "\t")
d_resp = as.data.frame(d_resp)

chemical = d_resp$chnm[1]
eff_pot = d_resp$Efficacy.potency[1]
fold_change = d_resp$CR.results[1]
  
d_resp$Concentration = log10(d_resp$Concentration)
d_resp$Response = log10(d_resp$Response)
d_resp$sd = sd(d_resp$Response)

Anova_pval = rep(paste("pval <", pval), length(d_resp$Response))
Anova_pval[which(d_resp$Anova <= pval)] = paste("pval >", pval)

ggplot(data = d_resp, aes(x = Concentration, y = Response)) +
  geom_point(aes(colour = Anova_pval)) +
  stat_summary(fun=mean, geom="line", colour="black", size=0.5, shape=4, linetype = "dashed") +
  labs(x = "Log Concentration Chemical (µM)", y = paste("Log Concentration", hormone, "(µM)", sep = " "), 
       title=paste(chemical, " (", basename(p_response) ,")", "\nEfficacy/potency: ", eff_pot, "\nFold change (Kaumas 2012): ", fold_change, sep=""))
  

ggsave(paste(p_response, ".png", sep = ""),  width = 5, height = 5, dpi = 300, bg="transparent")





