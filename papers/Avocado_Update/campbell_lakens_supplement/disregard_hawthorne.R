library(gee)
library(geepack)
library(BayesFactor)
library(RCurl)
library(httr)

#############################################################
# This loads custom code for ANOVA non-inferiority testing from github repo:
script <- getURL("https://raw.githubusercontent.com/harlanhappydog/noninfANOVAlm/master/noninfANOVA.R", ssl.verifypeer = FALSE)
eval(parse(text = script))


#############################################################
## This reads in the data, and formats it appropriately:
Hdata <-read.csv(text= getURL("https://raw.githubusercontent.com/harlanhappydog/noninfANOVAlm/master/analysisdataset2018.csv"))

Hdata$group<-as.factor(Hdata$group)
Hdata$t<-as.factor(Hdata$t)
Hdata$participant_ID <-as.factor(Hdata$participant_ID)

Hdata<- Hdata[!(Hdata$group)=="",]
Hdata$group<-as.factor(as.character(Hdata$group))

Hdata<-Hdata[order(Hdata$participant_ID),]
Hdata<-Hdata[Hdata$itt==1,]

#############################################################
## This reads reduces the data to subset of complete cases (i.e., participants that have both baseline and follow-up recorded for at least one outcome) :

Hdata_complete<-Hdata[Hdata$participant_ID%in%(c(names(table(Hdata$participant_ID)[table(Hdata$participant_ID)==2]))),]

Hdata_complete$participant_ID<-as.factor(as.character((Hdata_complete$participant_ID)))
Hdata_complete<-Hdata_complete[order(Hdata_complete$participant_ID),]

Hdata_complete <- Hdata_complete[order(Hdata_complete$participant_ID),]



base_data<-na.omit(Hdata[Hdata$t=="Baseline",][,c("participant_ID","totaldrinking", "group")])

followup_data<-na.omit(Hdata[Hdata$t=="Followup",][,c("participant_ID","totaldrinking")])

side_data<-merge(base_data , followup_data, by="participant_ID", all=TRUE)
side_data$totaldrinking.diff<-side_data$totaldrinking.y-side_data$totaldrinking.x


#############################################################
## This creates a "wide" version of the complete cases data:



#### sample sizes don't quite line up with what is published for Followup outcomes:
library(xtable)


### Group A
n_A <- length(na.omit(side_data[side_data$group=="A",]$totaldrinking.x))
m_A <- mean(na.omit(side_data[  side_data$group=="A",]$totaldrinking.x))
sd_A <- sd(na.omit(side_data[  side_data$group=="A",]$totaldrinking.x))


### Group B
n_B <- length(na.omit(side_data[  side_data$group=="B",]$totaldrinking.x))
m_B <- mean(na.omit(side_data[  side_data$group=="B",]$totaldrinking.x))
sd_B <- sd(na.omit(side_data[  side_data$group=="B",]$totaldrinking.x))

### Group C
n_C <- length(na.omit(side_data[  side_data$group=="C",]$totaldrinking.x))
m_C <- mean(na.omit(side_data[  side_data$group=="C",]$totaldrinking.x))
sd_C <- sd(na.omit(side_data[  side_data$group=="C",]$totaldrinking.x))

### Total
n_Total<- length(na.omit(side_data[,]$totaldrinking.x))
m_Total <- mean(na.omit(side_data[,]$totaldrinking.x))
sd_Total <- sd(na.omit(side_data[,]$totaldrinking.x))


baseline <-c(c(n_A,m_A,sd_A),c(n_B,m_B,sd_B),c(n_C,m_C,sd_C),c(n_Total,m_Total,sd_Total))
### Group A
n_A <- length(na.omit(side_data[side_data$group=="A",]$totaldrinking.y))
m_A <- mean(na.omit(side_data[  side_data$group=="A",]$totaldrinking.y))
sd_A <- sd(na.omit(side_data[  side_data$group=="A",]$totaldrinking.y))


### Group B
n_B <- length(na.omit(side_data[  side_data$group=="B",]$totaldrinking.y))
m_B <- mean(na.omit(side_data[  side_data$group=="B",]$totaldrinking.y))
sd_B <- sd(na.omit(side_data[  side_data$group=="B",]$totaldrinking.y))

### Group C
n_C <- length(na.omit(side_data[  side_data$group=="C",]$totaldrinking.y))
m_C <- mean(na.omit(side_data[  side_data$group=="C",]$totaldrinking.y))
sd_C <- sd(na.omit(side_data[  side_data$group=="C",]$totaldrinking.y))

### Total
n_Total<- length(na.omit(side_data[,]$totaldrinking.y))
m_Total <- mean(na.omit(side_data[,]$totaldrinking.y))
sd_Total <- sd(na.omit(side_data[,]$totaldrinking.y))


followup <-c(c(n_A,m_A,sd_A),c(n_B,m_B,sd_B),c(n_C,m_C,sd_C),c(n_Total,m_Total,sd_Total))


### Group A
n_A <- length(na.omit(side_data[side_data$group=="A",]$totaldrinking.diff))
m_A <- mean(na.omit(side_data[  side_data$group=="A",]$totaldrinking.diff))
sd_A <- sd(na.omit(side_data[  side_data$group=="A",]$totaldrinking.diff))


### Group B
n_B <- length(na.omit(side_data[  side_data$group=="B",]$totaldrinking.diff))
m_B <- mean(na.omit(side_data[  side_data$group=="B",]$totaldrinking.diff))
sd_B <- sd(na.omit(side_data[  side_data$group=="B",]$totaldrinking.diff))

### Group C
n_C <- length(na.omit(side_data[  side_data$group=="C",]$totaldrinking.diff))
m_C <- mean(na.omit(side_data[  side_data$group=="C",]$totaldrinking.diff))
sd_C <- sd(na.omit(side_data[  side_data$group=="C",]$totaldrinking.diff))

### Total
n_Total<- length(na.omit(side_data[,]$totaldrinking.diff))
m_Total <- mean(na.omit(side_data[,]$totaldrinking.diff))
sd_Total <- sd(na.omit(side_data[,]$totaldrinking.diff))


diff<-c(c(n_A,m_A,sd_A),c(n_B,m_B,sd_B),c(n_C,m_C,sd_C),c(n_Total,m_Total,sd_Total))

table_for_paper<-t(rbind(baseline, followup, diff))

rownames(table_for_paper)<-c(c("n_A","m_A","sd_A"),c("n_B","m_B","sd_B"),c("n_C","m_C","sd_C"),c("n_Total","m_Total","sd_Total"))

xtable(table_for_paper)




#### Analysis:


Xmatrix <- model.matrix(totaldrinking.diff  ~ group, data= side_data)
lmmodel <- lm(totaldrinking.diff  ~ group , data= side_data)

R2 <- summary(lmmodel)$r.squared
Fstat <- summary(lmmodel)$fstatistic[1]
K <- dim(Xmatrix)[2] - 1
N <- dim(Xmatrix)[1]
Delta <- 0.01
 
pf(Fstat,df1=K,df2=N-K-1,ncp=(N*Delta)/(1-Delta),lower.tail=TRUE)

linearReg.R2stat(N=N, p=K, R2= R2, simple=TRUE)




### The code bellow replicates the results published in McCambridge et al. (2019)
### Note: there appears to be a typo in McCambridge et al. (2019) Table 2: 
### p-values 0.89 and 0.86 are switched.

Hdata$group<-relevel(Hdata$group,"A")

mod0 <- geeglm(totaldrinking ~ + group+t, id= participant_ID, corstr="independence", data= Hdata, x=TRUE)
mod1 <- geeglm(totaldrinking ~ group*t + group+t, id= participant_ID, corstr="independence", data= Hdata)
(anova(mod1,mod0))
summary(mod1)$coefficients

Hdata$group<-relevel(Hdata$group,"C")
mod1a <- geeglm(totaldrinking ~ group*t + group+t, id= participant_ID, corstr="independence", data= Hdata)
summary(mod1a)


