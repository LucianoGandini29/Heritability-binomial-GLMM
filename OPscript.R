# Quantitative genetics analysis of oviposition preference
# Mixed-effects binomial models (glmmTMB)
# Author: Luciano Gandini
# This script reproduces the analyses presented in: 
#"Fanara, J. J., Beti, M. I., Gandini, L., & Hasson, E. (2023). 
# Oviposition behaviour in Drosophila melanogaster: 
# Genetic and behavioural decoupling between oviposition 
# acceptance and preference for natural fruits. 
#Journal of Evolutionary Biology, 36(1), 251-263."


rm(list=ls())
library(glmmTMB)
library(DHARMa)
library(multcomp)
library(performance)
library(MuMIn)
library(ggplot2)
library(car)
library(insight)
library(ResourceSelection)
library(VCA)
library(GLMMmisc)

# Set working directory to project root (adjust if needed)
# setwd("path/to/project/root")
getwd()

#############################OP################
OP_data<-read.csv("data/OP.csv", sep=",")

#Descriptive plots 

ggplot(OP_data, aes(x=Fruit_Failure, y=Prop_Successes)) + geom_point(size=3)


varPlot(form=Prop_Successes~(Fruit_Failure+Line), Data=OP_data)


str(OP_data)
OP_data$Line<-as.factor(OP_data$Line)
OP_data$ID<-as.factor(OP_data$ID)
OP_data$Fruit_Failure<-as.factor(OP_data$Fruit_Failure)


# Fit a binomial glmm with the proportion of successes (oviposition in grape) in n events (total number of ovipositions) as dependent variables,
#The alternative fruit (grape or tomato) as a fixed effect explanatory variable and for random effect variables we use "DGRP Line" and 
#"Line by alternative fruit" interaction

m0<-glmmTMB(Prop_Successes~Fruit_Failure+(1|Line)+(1|Line:Fruit_Failure), family = binomial, weights = Sum_Total, data = OP_data)

#We use DHARMA package to check residuals

simulationOutput <- simulateResiduals(fittedModel = m0, refit = F, plot = T)

check_overdispersion(m0)
#There are multiple deviations from model assumptions likely caused by overdispersion in the data
#We add an observation level random effect (OLRE) variable to account for this

m1<-glmmTMB(Prop_Successes~Fruit_Failure+(1|Line)+(1|Line:Fruit_Failure)+(1|ID), family = binomial, weights = Sum_Total, data = OP_data)
summary(m1)

#We use DHARMA package to check residuals

simulationOutput <- simulateResiduals(fittedModel = m1, refit = F, plot = T)

check_overdispersion(m1)

#Actual vs. fitted values plot

plot(fitted(m1), OP_data$Prop_Successes, ylab="Actual", xlab="Fitted",main="Observed vs. Fitted values" )

#Goodness of fit
hoslem.test(OP_data$Prop_Successes, fitted(m1), g=10)

#We get a good fit and can now analize the results

#We determine significance of fixed effects by a type II Wald test
Anova(m1)

#We extract variance components and calculate Variance partition coefficients (VPCs)
var<-get_variance(m1)

Variance_Components<-cbind(var$var.intercept[[2]],var$var.intercept[[1]],var$var.intercept[[3]],(var$var.intercept[[2]]+var$var.intercept[[1]]+var$var.intercept[[3]]))
colnames(Variance_Components)<-c("Line_by_Resource","Line","Residual","Total_Variance")
Variance_Components<-rbind(Variance_Components, c(0,0,0,0))
row.names(Variance_Components)<-c("Variance", "%")

for(i in 1:ncol(Variance_Components)) {
  Variance_Components[2,i]<-(Variance_Components[1,i]/Variance_Components[1,4])*100
}

write.table(Variance_Components,file="results/VPC_OP.tsv",sep="\t")

#########################################################################################
# Because the previous model cannot differentiate the variance explained by the Line in each contrast ("Orange vs grape" 
# and "Tomato vs Grape"), we need to analize both contrasts separately in order to get the Variance components

datO<-subset(OP_data, Fruit_Failure=="Orange")
datT<-subset(OP_data, Fruit_Failure=="Tomato")

#Orange
m0O<-glmmTMB(Prop_Successes~1+(1|Line), family = binomial, weights = Sum_Total, data = datO)


simulationOutputN <- simulateResiduals(fittedModel = m0O, refit = F, plot = T)
check_overdispersion(m0O)

#Overdispersion, we add an OLRE

m1O<-glmmTMB(Prop_Successes~1+(1|Line)+ (1|ID), family = binomial, weights = Sum_Total, data = datO)
summary(m1O)

#We use DHARMA package to check residuals

simulationOutputO <- simulateResiduals(fittedModel = m1O, refit = F, plot = T)

check_overdispersion(m1O)
get_variance(m1O)

#Actual vs. fitted values plot

plot(fitted(m1O), datO$Prop_Successes, ylab="Actual", xlab="Fitted",main="Actual vs Fitted scatter-plot" )

#We make a null model (No Line effect)
m1O_null<-glmmTMB(Prop_Successes~1+ (1|ID), family = binomial, weights = Sum_Total, data = datO)

#We use a likelihood ratio test (LRT) implemented in the anova() function to compare both models
anova(m1O,m1O_null)

#Variance components
varO<-get_variance(m1O)
Variance_components_O<-cbind(varO$var.intercept[[1]],varO$var.intercept[[2]],(varO$var.intercept[[2]]+varO$var.intercept[[1]]))
colnames(Variance_components_O)<-c("Line","Residual","Total_Variance")
Variance_components_O<-rbind(Variance_components_O, c(0,0,0))
row.names(Variance_components_O)<-c("Variance", "%")

for(i in 1:ncol(Variance_components_O)) {
  Variance_components_O[2,i]<-(Variance_components_O[1,i]/Variance_components_O[1,3])*100
}

##################Tomato###############

m0T<-glmmTMB(Prop_Successes~1+(1|Line), family = binomial, weights = Sum_Total, data = datT)


simulationOutputT <- simulateResiduals(fittedModel = m0T, refit = F, plot = T)
check_overdispersion(m0T)
#Same as with Orange, we have overdispersion

m1T<-glmmTMB(Prop_Successes~1+(1|Line)+(1|ID), family = binomial, weights = Sum_Total, data = datT)

summary(m1T)


simulationOutputT <- simulateResiduals(fittedModel = m1T, refit = F, plot = T)

check_overdispersion(m1T)


p<-plot(fitted(m1T), datT$Prop_Successes, ylab="Actual", xlab="Fitted",main="Actual vs Fitted scatter-plot" )
ggsave("results/figures/Actual vs Fitted.png", plot = p, width = 8, height = 6, dpi = 300)

#null model with no Line effect
m1T_null<-glmmTMB(Prop_Successes~1+(1|ID), family = binomial, weights = Sum_Total, data = datT)

#LRT
anova(m1T,m1T_null)

#Variance components
varT<-get_variance(m1T)
Variance_components_T<-cbind(varT$var.intercept[[1]],varT$var.intercept[[2]],(varT$var.intercept[[2]]+varT$var.intercept[[1]]))
colnames(Variance_components_T)<-c("Line","Residual","Total_Variance")
Variance_components_T<-rbind(Variance_components_T, c(0,0,0))
row.names(Variance_components_T)<-c("Variance", "%")

for(i in 1:ncol(Variance_components_T)) {
  Variance_components_T[2,i]<-(Variance_components_T[1,i]/Variance_components_T[1,3])*100
}
