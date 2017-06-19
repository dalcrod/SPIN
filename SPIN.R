
## ANALISIS OF THE SPIN cohort DATA #################

#Load gmodels package
#install.packages(c("gmodels", "rpart", "randomForest", "knitr", "psych", "pgirmess", "Hmisc", "car", "fBasics", "sm", "foreign", "xlsx"))
#install.packages(c("car", "xlsx", "ggplot2", "gridExtra"), repos='http://cran.us.r-project.org')

library(grid)
#library(gmodels)
#library(knitr)
#library(psych)
#library(survival)
#library(Hmisc)
library(car)
#library(stats)
#library(fBasics)
#library(sm)
library(xlsx)
#library(pgirmess)
library(ggplot2)
library(gridExtra)
#library(pROC)
#library(randomForest)
#library(rpart)
#library(rpart.plot)
#library(pwr)

options(scipen=8, warn=-1)

printWithNumber = function(gg_plot, top_right=counter, resetcounter=FALSE) 
{
  plot(gg_plot)
  if (resetcounter==TRUE){
    counter <<- 0
  }
  counter <<- counter+1
  label = textGrob(top_right,
                     x = 0.98,  # right side
                     y = 0.98,   # top
                     just="right", 
                     hjust = NULL,
                     vjust = 1,
                     gp=gpar(fontsize=10, col="gray"))
    grid.draw(label)
  }



#### RETRIEVE DATA and ARRANGE VARIABLES######

#DATA in Santpau
setwd("//dspau.santpau.es/U/CarpUMNeurologia/BaseUdM/Export-analisis/Resultados")

data <- data.frame(read.xlsx ("//dspau.santpau.es/U/CarpUMNeurologia/BaseUdM/Export-analisis/Tablas exportadas/CohBasalGlobal.xlsx", 1))
data1a <- data.frame(read.xlsx ("//dspau.santpau.es/U/CarpUMNeurologia/BaseUdM/Export-analisis/Tablas exportadas/CCoh1aGlobalExportable.xlsx", 1))
data2a <- data.frame(read.xlsx ("//dspau.santpau.es/U/CarpUMNeurologia/BaseUdM/Export-analisis/Tablas exportadas/CCoh2aGlobalExportable.xlsx", 1))
data3a <- data.frame(read.xlsx ("//dspau.santpau.es/U/CarpUMNeurologia/BaseUdM/Export-analisis/Tablas exportadas/CCoh3aGlobalExportable.xlsx", 1))
data4a <- data.frame(read.xlsx ("//dspau.santpau.es/U/CarpUMNeurologia/BaseUdM/Export-analisis/Tablas exportadas/CCoh4aGlobalExportable.xlsx", 1))

#DATA in MACBOOK
#data <- data.frame(read.xlsx ("/Users/Daniel/Google Drive/WORK/ExportacionBaseUdM/CohBasalGlobal.xlsx", 1))
#setwd("/Users/Daniel/Google Drive/WORK/Projectes/02-PROYECTOS ACTIVOS/SPIN cohort/Analisis R")


data$STDIAGCODE <- factor(data$STDIAGCODE)
data1a$STDIAGCODE <- factor(data1a$STDIAGCODE)
data2a$STDIAGCODE <- factor(data2a$STDIAGCODE)
data3a$STDIAGCODE <- factor(data3a$STDIAGCODE)
data4a$STDIAGCODE <- factor(data4a$STDIAGCODE)
data1a$STDIAGCODE1a <- factor(data1a$STDIAGCODE1a)
data2a$STDIAGCODE2a <- factor(data2a$STDIAGCODE2a)
data3a$STDIAGCODE3a <- factor(data3a$STDIAGCODE3a)
data4a$STDIAGCODE4a <- factor(data4a$STDIAGCODE4a)
levels(data$STDIAGCODE) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
levels(data1a$STDIAGCODE) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
levels(data1a$STDIAGCODE1a) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
levels(data2a$STDIAGCODE) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
levels(data2a$STDIAGCODE2a) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
levels(data3a$STDIAGCODE) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
levels(data3a$STDIAGCODE3a) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
levels(data4a$STDIAGCODE) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
levels(data4a$STDIAGCODE4a) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
dx <- levels(data$STDIAGCODE)
#Set up Cut-off point for ABETA42, TAU and P-TAU

for (i in 1:length(data$DBCODE)) {
  data$ABETA42status[i] <- ifelse(data$CSFABETA42[i] < 550, "ABETA42pos",
                                  ifelse(data$CSFABETA42[i] >550, "ABETA42neg", NA))
  data$TAUstatus[i] <- ifelse (data$CSFTAU[i]>350, "TAUpos",
                               ifelse(data$CSFTAU[i] <350, "TAUneg", NA))
  data$PTAUstatus[i] <- ifelse (data$CSFPTAU[i]>61, "PTAUpos", 
                                ifelse(data$CSFPTAU[i] <61, "PTAUneg", NA))
  }

data$RatioTauAbeta <- data$CSFTAU / data$CSFABETA42
data$RatioTauPtau <- data$CSFTAU / data$CSFPTAU
data$RatioTausappbeta <- data$CSFTAU / data$CSFSAPPBETA
data$RatioPtausappbeta <- data$CSFPTAU / data$CSFSAPPBETA
data$RatioYKLsappbeta <- data$CSFYKL40 / data$CSFSAPPBETA
data$RatioAbetasappbeta <- data$CSFABETA42 / data$CSFSAPPBETA

for (i in 1:length(data$RatioTauAbeta)) {
  data$RatioAD[i] <- ifelse(data$RatioTauAbeta[i] > 0.52, "RatioAD+", "RatioAD-")
}


for (i in 1:length(data$APOE)) {
      data$APOE4[i] <- ifelse(data$APOE[i]=="22"|
                                    data$APOE[i]=="23"|
                                    data$APOE[i]=="33", "APOE4-",
                              ifelse(data$APOE[i]=="24"|
                                           data$APOE[i]=="34"|
                                           data$APOE[i]=="44", "APOE4+", NA))
}

myvars <- c("DBCODE", "STDIAGCODE", "FUPDATEbasal", "NPSDATE", "CSFDATE", "MRIDATE")
myvars1 <- c("DBCODE", "STDIAGCODE", "FUPDATE1a", "NPSDATE")
myvars2 <- c("DBCODE", "STDIAGCODE", "FUPDATE2a", "NPSDATE", "CSFDATE", "MRIDATE")
myvars3 <- c("DBCODE", "STDIAGCODE", "FUPDATE3a", "NPSDATE")
myvars4 <- c("DBCODE", "STDIAGCODE", "FUPDATE4a", "NPSDATE", "CSFDATE", "MRIDATE")

subdata <- data[myvars]
subdata$VISIT <- "BASELINE"
subdata1a <- data1a[myvars1]
subdata1a$CSFDATE <- NA
subdata1a$MRIDATE <- NA
names(subdata1a) <- myvars
subdata1a$VISIT <- "YEAR 1"
subdata2a <- data2a[myvars2]
names(subdata2a) <- myvars
subdata2a$VISIT <- "YEAR 2"
subdata3a <- data3a[myvars3]
subdata3a$CSFDATE <- NA
subdata3a$MRIDATE <- NA
names(subdata3a) <- myvars
subdata3a$VISIT <- "YEAR 3"
subdata4a <- data4a[myvars4]
names(subdata4a) <- myvars
subdata4a$VISIT <- "YEAR 4"

longdata <- rbind(subdata, subdata1a, subdata2a, subdata3a, subdata4a)
rm(subdata)
rm(subdata1a)
rm(subdata2a)
rm(subdata3a)
rm(subdata4a)
longdata <- longdata[!is.na(longdata$FUPDATE),]
levels(longdata$STDIAGCODE) <- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
names(longdata)<- c("DBCODE", "STDIAGCODE", "FUPDATE", "NPSDATE", "CSFDATE", "MRIDATE", "VISIT")

##TABLE SUMMARY####
attach(data)
quantvariables <-data.frame(AGE, EDUC, MMSE, FCSRTTOTALFREE, FCSRTTOTAL, FCSRTDELTOT, CERADLISTRECALL, DIGDIR, DIGREV, TMTA, TMTB, BNT, REYCOPY, PHONFLU, SEMFLU, VOSP, NPI.Q.TOTAL, CSFABETA42, AB1.42, CSFTAU, CSFPTAU, CSFYKL40, CSFSAPPBETA, CSFNFL, CSFPROGRANULIN)
propvariables <- data.frame(SEX,  APOE4, RatioAD, ABETA42status, TAUstatus, PTAUstatus, FBPVISUAL, FDGVISUAL, STDIAGCODE)

DATATABLE <- matrix(nrow=length(quantvariables)*4, ncol=length(dx))
colnames(DATATABLE) <- dx
NAMES<-vector()
for(i in 1:length(quantvariables)) {
VAR <- quantvariables[i]
NAMES <- c(NAMES, paste("Summary-",names(VAR), sep=""), paste(names(VAR),"-N", sep=""), paste(names(VAR),"-MEAN", sep=""), paste(names(VAR),"-SD", sep=""))
for (k in 1:length(dx)) {
  DATATABLE[i*4-3,k] <- paste("")
  DATATABLE[i*4-2,k] <- length(VAR[data$STDIAGCODE==dx[k]&!is.na(VAR),])
  DATATABLE[i*4-1,k] <- round(mean(VAR[data$STDIAGCODE==dx[k],],na.rm=TRUE), digits=1)
  DATATABLE[i*4,k] <- round(sd(VAR[data$STDIAGCODE==dx[k],],na.rm=TRUE), digits=1)
}
}
row.names(DATATABLE) <- NAMES
colnames(DATATABLE)<- c("Control", "SCI", "MCI", "AD", "LBD", "FTD", "Down", "ALS")
date<-paste("Data updated on", format(Sys.Date(), "%Y-%m-%d"))
NAMESprop <- vector()
proptable <- table(dx)
for(i in 1:length((propvariables))) {
      VAR <- propvariables[,i]
      a<-(table(VAR, propvariables$STDIAGCODE))
      proptable <- rbind(proptable, a)
}
PROPTABLE <- proptable[2:(nrow(proptable)-length(dx)),]

DATATABLE <- rbind(DATATABLE, PROPTABLE)
write.xlsx2(DATATABLE, "DATATABLE.xlsx", sheetName=date, showNA=FALSE)


Sys.setlocale(locale = "en")
caption <- paste("Summary data from the SPIN cohort, Hospital de la Santa Creu i Sant Pau, Barcelona\nLast update on", format(Sys.time(), "%A, %B %d %Y at %H:%M"))

palettenumber<- "BuPu"



blank <- ggplot(data, aes()) + 
  geom_blank() +
  scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0))+
  theme_minimal()


### SUMMARY PDF ######
######################
pdf("SPINsummary.pdf", 8, 6)
#### >>pdf section 1####
grid.text(label="SPIN cohort",x=0.5, y=0.6,just="centre", gp=gpar(fontsize=24, col="#177db7", fontface="bold"))
grid.text(label="Sant Pau Initiative on Neurodegeneration",x=0.5, y=0.5, gp=gpar(fontsize=14, col="black", fontface="bold"))
grid.text(label=caption, x=0.5, y=0.22, vjust=0, hjust=0.5, gp=gpar(fontsize=9, col="gray", fontface="bold"))

#### >>pdf section 1####
printWithNumber(blank, resetcounter = TRUE)
grid.text(label="1. Summary Table",hjust="centre", vjust="bottom", gp=gpar(fontsize=16, col="#177db7", fontface="bold"))

printWithNumber(blank)
table1 <- (tableGrob(DATATABLE[1:20,], theme=ttheme_minimal(base_size=7)))
grid.draw(table1)

printWithNumber(blank)
table2 <- (tableGrob(DATATABLE[21:40,], theme=ttheme_minimal(base_size=7)))
grid.draw(table2)

printWithNumber(blank)
table3 <- (tableGrob(DATATABLE[41:60,], theme=ttheme_minimal(base_size=7)))
grid.draw(table3)

printWithNumber(blank)
table4 <- (tableGrob(DATATABLE[61:80,], theme=ttheme_minimal(base_size=7)))
grid.draw(table4)

printWithNumber(blank)
table5 <- (tableGrob(DATATABLE[81:100,], theme=ttheme_minimal(base_size=7)))
grid.draw(table5)

printWithNumber(blank)
table6 <- (tableGrob(DATATABLE[101:nrow(DATATABLE),], theme=ttheme_minimal(base_size=7)))
grid.draw(table6)

#### >>pdf section 2###
printWithNumber(blank)
grid.text(label="2. Cohort Summary Counts - BASELINE",hjust="centre", vjust="bottom", gp=gpar(fontsize=16, col="#177db7", fontface="bold"))

#### COHORT STATUS GRAPHICS ###
### Summary graphs#####
cohortstatus <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE),
             fill=SEX)) +
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - Count", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
      theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortstatus)


plotages <-
ggplot(data=data, 
       aes(x=factor(data$STDIAGCODE), y=AGE, fill=factor(data$STDIAGCODE))) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Age (years)", title="SPIN Cohort - Age", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_y_continuous(limits=c(0, 110), breaks=seq(0, 100, 20))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(plotages)


ploteduc <-
ggplot(data=data, 
       aes(x=factor(data$STDIAGCODE), y=EDUC, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="", y="Education (years)", title="SPIN Cohort - Education", caption=caption)+
  guides(fill=FALSE)+
  scale_y_continuous(limits=c(0, 25), breaks=seq(0, 20, 2))+
  theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(ploteduc)


plotmmse <-
      ggplot(data=data, 
             aes(x=factor(data$STDIAGCODE), y=MMSE, fill=factor(data$STDIAGCODE), legend)) +
      geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
      labs(x="",y="Mini-Mental State Examination", title="SPIN Cohort - MMSE", caption=caption)+
      guides(fill=FALSE)+
      theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
      scale_y_continuous(limits=c(0, 32), breaks=seq(0, 30, 2))+
      stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(plotmmse)

cohortapoe <- 
      ggplot(data=data, 
             aes(x=factor(data$STDIAGCODE),
                 fill=APOE4)) +
      geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
      geom_bar(colour="black", alpha=0.4, position="dodge")+
      labs(x="", y="Count", title="SPIN Cohort - APOE4", caption=caption)+
      scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
      scale_fill_brewer(palette=palettenumber)+
      theme_classic()+
      theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
      geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
      geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortapoe)



cohortabetastatus <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE),
             fill=ABETA42status)) +
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - ABETA42 status (cut-off=550 pg/ml)", fill="ABETA42", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortabetastatus)


cohorttaustatus <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE),
             fill=TAUstatus)) +
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - T-TAU status (cut-off=350 pg/ml)", fill="T-TAU", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohorttaustatus)


cohortptaustatus <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE),
             fill=PTAUstatus)) +
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - P-TAU status (cut-off=61 pg/ml)", fill="P-TAU", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortptaustatus)


cohortinclusion <- 
  ggplot(data, aes(y=STDIAGCODE,
                   x=FUPDATEbasal,
                   colour=(STDIAGCODE))) +
  geom_point(size=3, alpha=0.7)+
  labs(x="Baseline date",y="Diagnostic group", title="SPIN Cohort - Inclusion", caption=caption, colour="")+
  scale_x_date(date_breaks = "4 months", date_labels = "%b-%Y")+
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"), axis.text.x = element_text(face="bold", hjust=0, angle=330))
printWithNumber(cohortinclusion)



cohortstatusMRI <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE),
             fill=!is.na(data$MRIDATE))) +
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - MRI Count", fill="MRI is available", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
      theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortstatusMRI)

cohortstatusCSF <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE),
             fill=!is.na(data$CSFDATE))) +
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - CSF Count", fill="CSF is available", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
      theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.85)
printWithNumber(cohortstatusCSF)

cohortstatusFBP <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE),
             fill=!is.na(data$FBPDATE))) +
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - 18F-Florbetapir PET Count", fill="FBP is available", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
      theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.80)
printWithNumber(cohortstatusFBP)

cohortstatusFDG <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE),
             fill=!is.na(data$FDGDATE))) +
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - 18F-Fluorodeoxyglucose PET Count", fill="FDG is available", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
      theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.80)
printWithNumber(cohortstatusFDG)


cohortstatusRatiopos <- 
  ggplot(data=data, 
         aes(x=factor(STDIAGCODE),
             fill=factor(RatioAD)))+
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - RatioAD in CSF", fill="RatioAD<0.52\n(AD profile)", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
      theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.80)
printWithNumber(cohortstatusRatiopos)


cohortstatusFBPpos <- 
  ggplot(data=data, 
         aes(x=factor(STDIAGCODE),
             fill=factor(FBPVISUAL)))+
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - Visual FBP-PET", fill="Positive FBP", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
      theme( plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.80)
printWithNumber(cohortstatusFBPpos)

cohortstatusFDGpos <- 
  ggplot(data=data, 
         aes(x=factor(STDIAGCODE),
             fill=factor(FDGVISUAL)))+
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - Visual FDG-PET", fill="Positive FDG", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
      theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.80)
printWithNumber(cohortstatusFDGpos)





printWithNumber(blank)
grid.text(label="3. Longitudinal Follow-Up", hjust="centre", vjust="bottom", gp=gpar(fontsize=16, col="#177db7", fontface="bold"))



cohortstatusfup <- 
  ggplot(data=longdata, 
         aes(x=factor(STDIAGCODE),
             fill=factor(VISIT)))+
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - Follow-Up Count", fill="", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortstatusfup)


cohortstatusfupnps <- 
  ggplot(data=longdata[!is.na(longdata$NPSDATE),], 
         aes(x=factor(STDIAGCODE),
             fill=factor(VISIT)))+
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - Follow-Up NPS Count", fill="", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortstatusfupnps)

cohortstatusfupcsf <- 
  ggplot(data=longdata[!is.na(longdata$CSFDATE),], 
         aes(x=factor(STDIAGCODE),
             fill=factor(VISIT)))+
  geom_bar(data=data, aes(x=factor(data$STDIAGCODE), fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - Follow-Up CSF Count", fill="", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortstatusfupcsf)

cohortstatusfupmri <- 
  ggplot(data=longdata[!is.na(longdata$MRIDATE),], 
         aes(x=STDIAGCODE,
             fill=VISIT))+
  geom_bar(data=data, aes(x=data$STDIAGCODE, fill="", alpha=0.1), show.legend=FALSE)+
  geom_bar(colour="black", alpha=0.4, position="dodge")+
  labs(x="", y="Count", title="SPIN Cohort - Follow-Up MRI Count", fill="", caption=caption)+
  scale_y_continuous(limits=c(0, 315), breaks=seq(0, 300, 20))+
  scale_fill_brewer(palette=palettenumber)+
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  geom_text(stat="count", aes(label=..count..),cex=3,col="gray", position = position_dodge(width = 1), show.legend=FALSE, vjust = -0.6)+
  geom_text(data=data, aes(x=factor(data$STDIAGCODE), fill="",label=..count..),cex=6, show.legend=FALSE, stat="count", position = position_dodge(width = 0.8), vjust = -0.8)
printWithNumber(cohortstatusfupmri)





printWithNumber(blank)
grid.text(label="4. Neuropsychology - BASELINE", hjust="centre", vjust="bottom", gp=gpar(fontsize=16, col="#177db7", fontface="bold"))


cohortfcsrttotal <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=FCSRTTOTAL, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="FCSRT total", title="SPIN Cohort - FCSRT total", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortfcsrttotal)


cohortcerad <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=CERADLISTRECALL, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="CERAD list recall", title="SPIN Cohort - CERAD list recall", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortcerad)


cohortdigdir <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=DIGDIR, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Digits direct", title="SPIN Cohort - Digits direct", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortdigdir)


cohortdigrev <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=DIGREV, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Digits reverse", title="SPIN Cohort - Digits reverse", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortdigrev)


cohorttmta <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=TMTA, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Trail Making Test A", title="SPIN Cohort - Trail Making Test A", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohorttmta)

cohorttmtb <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=TMTB, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Trail Making Test B", title="SPIN Cohort - Trail Making Test B", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohorttmtb)


cohortbnt <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=BNT, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Boston Naming Test (60)", title="SPIN Cohort - Boston Naming Test", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortbnt)

cohortreycopy <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=REYCOPY, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Rey Complex Figure copy", title="SPIN Cohort - Rey Complex Figure copy", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortreycopy)


cohortPHONFLU <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=PHONFLU, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Phonetic Fluency", title="SPIN Cohort - Phonetic Fluency", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortPHONFLU)


cohortSEMFLU <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=SEMFLU, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Semantic Fluency", title="SPIN Cohort - Semantic Fluency", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortSEMFLU)


cohortVOSP <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=VOSP, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="VOSP", title="SPIN Cohort - VOSP", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortVOSP)


cohortNPI <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=NPITOTAL, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Neuropsychiatric Inventory", title="SPIN Cohort - Neuropsychiatric Inventory", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortNPI)

cohortNPIq <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=NPI.Q.TOTAL, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Neuropsychiatric Inventory Q-Total", title="SPIN Cohort - Neuropsychiatric Inventory Q-Total", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  #scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortNPIq)


printWithNumber(blank)
grid.text(label="5. CSF Biomarkers - BASELINE", hjust="centre", vjust="bottom", gp=gpar(fontsize=16, col="#177db7", fontface="bold"))


cohortabeta <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=CSFABETA42, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Abeta42 (pg/ml)", title="SPIN Cohort - CSF Abeta42", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortabeta)

cohortabetaage <- 
  ggplot(data, aes(y=CSFABETA42,
                   x=AGE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="ABeta42 (pg/ml)",y="Age (years)", title="SPIN Cohort - CSF Abeta42 - Age", caption=caption, colour="")+
  geom_hline(yintercept=550, linetype=2, colour="gray")+
  xlab("Age (years)") + xlim(0,100) +
  ylab("Abeta42 (pg/ml)") + ylim(0,2000) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortabetaage)

cohortabetammse <- 
  ggplot(data, aes(y=CSFABETA42,
                   x=MMSE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="ABeta42 (pg/ml)",y="Mini-Mental State Examination", title="SPIN Cohort - CSF Abeta42 - MMSE", caption=caption, colour="")+
  geom_hline(yintercept=550, linetype=2, colour="gray")+
  xlab("MMSE score") + xlim(0,30) +
  ylab("Abeta42 (pg/ml)") + ylim(0,2000) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortabetammse)


cohorttau <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=CSFTAU, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="Abeta42 (pg/ml)", title="SPIN Cohort - CSF T-Tau", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_y_continuous(limits=c(0, 2200), breaks=seq(0, 2000, 200))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,   vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohorttau)

cohorttauage <- 
  ggplot(data, aes(y=CSFTAU,
                   x=AGE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="T-Tau (pg/ml)",y="Age (years)", title="SPIN Cohort - CSF T-Tau - Age", caption=caption, colour="")+
  geom_hline(yintercept=350, linetype=2, colour="gray")+
    xlab("Age (years)") + xlim(0,100) +
  ylab("T-Tau (pg/ml)") + ylim(0,2500) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohorttauage)


cohorttaummse <- 
  ggplot(data, aes(y=CSFTAU,
                   x=MMSE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="T-Tau (pg/ml)",y="Mini-Mental State Examination", title="SPIN Cohort - CSF T-Tau - MMSE", caption=caption, colour="")+
  geom_hline(yintercept=350, linetype=2, colour="gray")+
  xlab("MMSE score") + xlim(0,30) +
  ylab("T-Tau (pg/ml)") + ylim(0,2300) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohorttaummse)


cohortptau <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=CSFPTAU, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="P-Tau (pg/ml)", title="SPIN Cohort - CSF P-Tau", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_y_continuous(limits=c(0, 300), breaks=seq(0, 300, 25))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE, vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortptau)

cohortptauage <- 
  ggplot(data, aes(y=CSFPTAU,
                   x=AGE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="P-Tau (pg/ml)",y="Age (years)", title="SPIN Cohort - CSF P-Tau - Age", caption=caption, colour="")+
  geom_hline(yintercept=61, linetype=2, colour="gray")+
  xlab("Age (years)") + xlim(0,100) +
  ylab("P-Tau (pg/ml)") + ylim(0,300) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortptauage)


cohortptaummse <- 
  ggplot(data, aes(y=CSFPTAU,
                   x=MMSE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="P-Tau (pg/ml)",y="Mini-Mental State Examination", title="SPIN Cohort - CSF P-Tau - MMSE", caption=caption, colour="")+
  geom_hline(yintercept=61, linetype=2, colour="gray")+
  xlab("MMSE score") + xlim(0,30) +
  ylab("P-Tau (pg/ml)") + ylim(0,300) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortptaummse)


cohortabetatau <- 
  ggplot(data, aes(y=CSFABETA42,
                   x=CSFTAU,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="ABeta42 (pg/ml)",y="T-Tau (pg/ml)", title="SPIN Cohort - CSF Abeta42 - T-Tau", caption=caption, colour="")+
  geom_hline(yintercept=550, linetype=2, colour="gray")+
  geom_vline(xintercept=350, linetype=2, colour="gray")+
  xlab("Total-Tau (pg/ml)") + xlim(0,2500) +
  ylab("Abeta42 (pg/ml)") + ylim(0,2000) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortabetatau)

cohortabetaptau <- 
  ggplot(data, aes(y=CSFABETA42,
                   x=CSFPTAU,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="ABeta42 (pg/ml)",y="P-Tau (pg/ml)", title="SPIN Cohort - CSF Abeta42 - P-Tau", caption=caption, colour="")+
  geom_hline(yintercept=550, linetype=2, colour="gray")+
  geom_vline(xintercept=61, linetype=2, colour="gray")+
  xlab("P-Tau (pg/ml)") + xlim(0,300) +
  ylab("Abeta42 (pg/ml)") + ylim(0,2000) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortabetaptau)


cohorttauptau <- 
  ggplot(data, aes(y=CSFTAU,
                   x=CSFPTAU,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="T-Tau (pg/ml)",y="P-Tau (pg/ml)", title="SPIN Cohort - CSF T-Tau - P-Tau", caption=caption, colour="")+
  geom_hline(yintercept=350, linetype=2, colour="gray")+
  geom_vline(xintercept=61, linetype=2, colour="gray")+
  xlab("P-Tau (pg/ml)") + xlim(0,300) +
  ylab("T-Tau (pg/ml)") + ylim(0,2300) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohorttauptau)


cohortykl <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=CSFYKL40, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="YKL-40 (ng/ml)", title="SPIN Cohort - CSF YKL-40", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_y_continuous(limits=c(0, 500), breaks=seq(0, 500, 50))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortykl)

cohoryklage <- 
  ggplot(data, aes(y=CSFYKL40,
                   x=AGE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="YKL-40 (ng/ml)",y="Age (years)", title="SPIN Cohort - CSF YKL-40 - Age", caption=caption, colour="")+
  xlab("Age (years)") + xlim(0,100) +
  ylab("YKL-40 (ng/ml)") + 
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohoryklage)


cohortyklmmse <- 
  ggplot(data, aes(y=CSFYKL40,
                   x=MMSE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="YKL-40 (ng/ml)",y="Mini-Mental State Examination", title="SPIN Cohort - CSF YKL-40  - MMSE", caption=caption, colour="")+
  xlab("MMSE score") + xlim(0,30) +
  ylab("YKL-40  (ng/ml)") +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortyklmmse)


cohortabetaykl <- 
  ggplot(data, aes(y=CSFABETA42,
                   x=CSFYKL40,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="Abeta42 (pg/ml)",y="YKL-40 (ng/ml)", title="SPIN Cohort - CSF Abeta42 - YKL-40", caption=caption, colour="")+
  geom_hline(yintercept=550, linetype=2, colour="gray")+
  xlab("YKL-40 (ng/ml)") +
  ylab("Abeta42 (pg/ml)") + ylim(0,2000) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortabetaykl)

cohorttauykl <- 
  ggplot(data, aes(y=CSFTAU,
                   x=CSFYKL40,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="T-Tau (pg/ml)",y="YKL-40 (ng/ml)", title="SPIN Cohort - CSF T-Tau - YKL-40", caption=caption, colour="")+
  geom_hline(yintercept=350, linetype=2, colour="gray")+
  xlab("YKL-40 (ng/ml)") +
  ylab("T-Tau (pg/ml)") + ylim(0,2300) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohorttauykl)




cohortnfl <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=CSFNFL, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="NFL (pg/ml)", title="SPIN Cohort - CSF Neurofilaments light", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
       scale_y_continuous(limits=c(0, 6000), breaks=seq(0, 6000, 500))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortnfl)

cohortnflage <- 
  ggplot(data, aes(y=CSFNFL,
                   x=AGE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="NFL (ng/ml)",y="Age (years)", title="SPIN Cohort - CSF NFL - Age", caption=caption, colour="")+
  xlab("Age (years)") + xlim(0,100) +
  ylab("NFL (ng/ml)") + 
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortnflage)


cohortnflmmse <- 
  ggplot(data, aes(y=CSFNFL,
                   x=MMSE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="NFL (ng/ml)",y="Mini-Mental State Examination", title="SPIN Cohort - CSF NFL  - MMSE", caption=caption, colour="")+
  xlab("MMSE score") + xlim(0,30) +
  ylab("NFL  (ng/ml)") +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortnflmmse)




cohortabetanfl <- 
  ggplot(data, aes(y=CSFABETA42,
                   x=CSFNFL,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="Abeta42 (pg/ml)",y="NFL (pg/ml)", title="SPIN Cohort - CSF Abeta42 - NFL", caption=caption, colour="")+
  geom_hline(yintercept=550, linetype=2, colour="gray")+
  xlab("NFL (ng/ml)") +
  ylab("Abeta42 (pg/ml)") + ylim(0,2000) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortabetanfl)

cohorttaunfl <- 
  ggplot(data, aes(y=CSFTAU,
                   x=CSFNFL,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="T-Tau (pg/ml)",y="NFL (pg/ml)", title="SPIN Cohort - CSF T-Tau - NFL", caption=caption, colour="")+
  geom_hline(yintercept=350, linetype=2, colour="gray")+
  xlab("NFL (pg/ml)") +
  ylab("T-Tau (pg/ml)") + ylim(0,2300) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohorttaunfl)



cohortsappb <- 
  ggplot(data=data, 
         aes(x=factor(data$STDIAGCODE), y=CSFSAPPBETA, fill=factor(data$STDIAGCODE), legend)) +
  geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
  labs(x="",y="sAPPbeta (ng/ml)", title="SPIN Cohort - CSF sAPPbeta", caption=caption)+
  guides(fill=FALSE)+
  theme_classic()+  
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_y_continuous(limits=c(0, 3000), breaks=seq(0, 3000, 500))+
  stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortsappb)


cohortsappbage <- 
  ggplot(data, aes(y=CSFSAPPBETA,
                   x=AGE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="sAPPbeta (ng/ml)",y="Age (years)", title="SPIN Cohort - CSF sAPPbeta - Age", caption=caption, colour="")+
  xlab("Age (years)") + xlim(0,100) +
  ylab("sAPPbeta (ng/ml)") + 
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortsappbage)


cohortsappbmmse <- 
  ggplot(data, aes(y=CSFSAPPBETA,
                   x=MMSE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="sAPPbeta (ng/ml)",y="Mini-Mental State Examination", title="SPIN Cohort - CSF sAPPbeta  - MMSE", caption=caption, colour="")+
  xlab("MMSE score") + xlim(0,30) +
  ylab("sAPPbeta  (ng/ml)") +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortsappbmmse)



cohortabetasappb <- 
  ggplot(data, aes(y=CSFABETA42,
                   x=CSFSAPPBETA,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="Abeta42 (pg/ml)",y="sAPPbeta (ng/ml)", title="SPIN Cohort - CSF Abeta42 - sAPPbeta", caption=caption, colour="")+
  geom_hline(yintercept=550, linetype=2, colour="gray")+
  xlab("sAPPbeta (ng/ml)") +
  ylab("Abeta42 (pg/ml)") + ylim(0,2000) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortabetasappb)

cohorttausappb <- 
  ggplot(data, aes(y=CSFTAU,
                   x=CSFSAPPBETA,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="T-Tau (pg/ml)",y="sAPPbeta (ng/ml)", title="SPIN Cohort - CSF T-Tau - sAPPbeta", caption=caption, colour="")+
  geom_hline(yintercept=350, linetype=2, colour="gray")+
  xlab("sAPPbeta (ng/ml)") +
  ylab("T-Tau (pg/ml)") + ylim(0,2300) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohorttausappb)


cohortprogranulin <- 
      ggplot(data=data, 
             aes(x=factor(data$STDIAGCODE), y=CSFPROGRANULIN, fill=factor(data$STDIAGCODE), legend)) +
      geom_boxplot(stat="boxplot", colour="black", alpha=0.4)+
      labs(x="",y="Progranulin (ng/ml)", title="SPIN Cohort - CSF Progranulin", caption=caption)+
      guides(fill=FALSE)+
      theme_classic()+
      theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
      scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, 2))+
      stat_summary(fun.y=mean, geom="text", show_guide = FALSE,                vjust=0, position="identity", aes(label=round(..y.., digits=1), fontface="bold"))
printWithNumber(cohortprogranulin)



cohortprogranulinage <- 
  ggplot(data, aes(y=CSFPROGRANULIN,
                   x=AGE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="Progranulin (ng/ml)",y="Age (years)", title="SPIN Cohort - CSF Progranulin - Age", caption=caption, colour="")+
  xlab("Age (years)") + xlim(0,100) +
  ylab("Progranulin (ng/ml)") + 
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortprogranulinage)


cohortprogranulinmmse <- 
  ggplot(data, aes(y=CSFPROGRANULIN,
                   x=MMSE,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="Progranulin (ng/ml)",y="Mini-Mental State Examination", title="SPIN Cohort - CSF Progranulin  - MMSE", caption=caption, colour="")+
  xlab("MMSE score") + xlim(0,30) +
  ylab("Progranulin  (ng/ml)") +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortprogranulinmmse)



cohortabetaprogranulin <- 
  ggplot(data, aes(y=CSFABETA42,
                   x=CSFPROGRANULIN,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="Abeta42 (pg/ml)",y="Progranulin (ng/ml)", title="SPIN Cohort - CSF Abeta42 - Progranulin", caption=caption, colour="")+
  geom_hline(yintercept=550, linetype=2, colour="gray")+
  xlab("Progranulin (ng/ml)") +
  ylab("Abeta42 (pg/ml)") + ylim(0,2000) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohortabetaprogranulin)

cohorttauprogranulin <- 
  ggplot(data, aes(y=CSFTAU,
                   x=CSFPROGRANULIN,
                   colour=(STDIAGCODE))) +
  geom_point(size=2, alpha=0.7)+
  labs(x="T-Tau (pg/ml)",y="Progranulin (ng/ml)", title="SPIN Cohort - CSF T-Tau - Progranulin", caption=caption, colour="")+
  geom_hline(yintercept=350, linetype=2, colour="gray")+
  xlab("Progranulin (ng/ml)") +
  ylab("T-Tau (pg/ml)") + ylim(0,2300) +
  theme_classic()+
  theme(plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))
printWithNumber(cohorttauprogranulin)


### Density plots ######
printWithNumber(blank)
grid.text(label="6. Biomarkers distribution - BASELINE", hjust="centre", vjust="bottom", gp=gpar(fontsize=16, col="#177db7", fontface="bold"))

for (i in 1:length(dx)) {
plot <- data[data$STDIAGCODE==dx[i],]

abeta <- round(quantile(plot$CSFABETA42, probs=c(0.05, 0.5, 0.95), na.rm=TRUE), digits=1)
abetamean <- round(mean(plot$CSFABETA42,na.rm=TRUE), digits=1)
abetasd <- round(sd(plot$CSFABETA42,na.rm=TRUE), digits=1)
labelabeta <- paste("ABETA42", "\n 5% Perc: ", abeta[1], "pg/ml\n50% Perc: ", abeta[2], "pg/ml\n95% Perc: ", abeta[3], "pg/ml\n\nMean: ",abetamean,"pg/ml\nSD: ",abetasd,"pg/ml", sep="")

abeta42 <- round(quantile(plot$AB1.42, probs=c(0.05, 0.5, 0.95), na.rm=TRUE), digits=1)
abeta42mean <- round(mean(plot$AB1.42,na.rm=TRUE), digits=1)
abeta42sd <- round(sd(plot$AB1.42,na.rm=TRUE), digits=1)
labelabeta42 <- paste("ABETA1-42 (Lumipulse)", "\n 5% Perc: ", abeta42[1], "pg/ml\n50% Perc: ", abeta42[2], "pg/ml\n95% Perc: ", abeta42[3], "pg/ml\n\nMean: ",abeta42mean,"pg/ml\nSD: ",abeta42sd,"pg/ml", sep="")


tau <- round(quantile(plot$CSFTAU, probs=c(0.05, 0.5, 0.95), na.rm=TRUE), digits=1)
taumean <- round(mean(plot$CSFTAU,na.rm=TRUE), digits=1)
tausd <- round(sd(plot$CSFTAU,na.rm=TRUE), digits=1)
labeltau <- paste("T-TAU", "\n 5% Perc: ", tau[1], "pg/ml\n50% Perc: ", tau[2], "pg/ml\n95% Perc: ", tau[3], "pg/ml\n\nMean: ", taumean,"pg/ml\nSD: ",tausd, "pg/ml",sep="")

ptau <- round(quantile(plot$CSFPTAU, probs=c(0.05, 0.5, 0.95), na.rm=TRUE), digits=1)
ptaumean <- round(mean(plot$CSFPTAU,na.rm=TRUE), digits=1)
ptausd <- round(sd(plot$CSFPTAU,na.rm=TRUE), digits=1)
labelptau <- paste("P-TAU", "\n 5% Perc: ", ptau[1], "pg/ml\n50% Perc: ", ptau[2], "pg/ml\n95% Perc: ", ptau[3], "pg/ml\n\nMean: ",ptaumean,"pg/ml\nSD: ",ptausd, "pg/ml",sep="")

age <- round(quantile(plot$AGE, probs=c(0.05, 0.5, 0.95), na.rm=TRUE), digits=1)
agemean <- round(mean(plot$AGE,na.rm=TRUE), digits=1)
agesd <- round(sd(plot$AGE,na.rm=TRUE), digits=1)
labelage <- paste("AGE", "\n 5% Perc: ", age[1],"years\n50% Perc: ", age[2],  "years\n95% Perc: ", age[3], "years\n\nMean: ",agemean,"years\nSD: ",agesd, "years",sep="")


plotabeta <- 
  ggplot(plot, aes(x=plot$CSFABETA42, fill=""))+
  #geom_histogram(alpha=0.2, bins=30, colour="black")+
  geom_density(alpha=0.2)+
  xlab("Abeta42 (pg/ml)")+
  ylab("Density estimate")+
  geom_vline(xintercept= abeta[2]) +
  geom_vline(xintercept= abeta[1], linetype=2) +
  geom_vline(xintercept= abeta[3], linetype=2) +
  xlim(0, 2000)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_fill_manual(values="green")+
  labs(title=paste("DIAGNOSTIC GROUP: ", dx[i]), subtitle="")+
  labs(caption=caption)+
  annotate(geom="text", x=Inf,y=Inf, label=labelabeta, size=3, hjust=1, vjust=1)
  

plotabeta42 <- 
  ggplot(plot, aes(x=plot$AB1.42, fill=""))+
  #geom_histogram(alpha=0.2, bins=30, colour="black")+
  geom_density(alpha=0.2)+
  xlab("Abeta42 - Lumipulse (pg/ml)")+
  ylab("Density estimate")+
  geom_vline(xintercept= abeta42[2]) +
  geom_vline(xintercept= abeta42[1], linetype=2) +
  geom_vline(xintercept= abeta42[3], linetype=2) +
  xlim(0, 3000)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_fill_manual(values="green")+
  labs(title=paste("DIAGNOSTIC GROUP: ", dx[i]), subtitle="")+
  labs(caption=caption)+
  annotate(geom="text", x=Inf,y=Inf, label=labelabeta42, size=3, hjust=1, vjust=1)

plottau <- 
  ggplot(plot, aes(x=plot$CSFTAU, fill=""))+
  #geom_histogram(alpha=0.2, bins=30, colour="black")+
  geom_density(alpha=0.2)+
  xlab("T-Tau (pg/ml)")+
  ylab("Density estimate")+
  geom_vline(xintercept= tau[2]) +
  geom_vline(xintercept= tau[1], linetype=2) +
  geom_vline(xintercept= tau[3], linetype=2) +
  xlim(0,2000)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_fill_manual(values="yellow")+
  labs(title=paste("DIAGNOSTIC GROUP: ", dx[i]), subtitle="")+
  labs(caption=caption)+
  annotate(geom="text", x = Inf, y=Inf, label=labeltau,size=3, hjust=1, vjust=1)

plotptau <- 
  ggplot(plot, aes(x=plot$CSFPTAU, fill=""))+
  #geom_histogram(alpha=0.2, bins=30, colour="black")+
  geom_density(alpha=0.2)+
  xlab("P-Tau (pg/ml)")+
  ylab("Density estimate")+
  geom_vline(xintercept= ptau[2]) +
  geom_vline(xintercept= ptau[1], linetype=2) +
  geom_vline(xintercept= ptau[3], linetype=2) +
  xlim(0,250)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_fill_manual(values="purple")+
  labs(title=paste("DIAGNOSTIC GROUP: ", dx[i]), subtitle="")+
  labs(caption=caption)+
  annotate(geom="text", x = Inf, y=Inf, label=labelptau,size=3, hjust=1, vjust=1)

plotage <- 
  ggplot(plot, aes(x=plot$AGE, fill=""))+
  #geom_histogram(alpha=0.2, bins=30, colour="black")+
  geom_density(alpha=0.2)+
  xlab("Age (years)")+
  ylab("Density estimate")+
  geom_vline(xintercept= age[2]) +
  geom_vline(xintercept= age[1], linetype=2) +
  geom_vline(xintercept= age[3], linetype=2) +
  xlim(0, 100)+
  theme_classic()+
  theme(legend.position="none", plot.title = element_text(face="bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5), plot.caption=element_text(size=6, face="italic"))+
  scale_fill_manual(values="blue")+
  labs(title=paste("DIAGNOSTIC GROUP: ", dx[i]), subtitle="")+
  labs(caption=caption)+
  annotate(geom="text", x = Inf, y=Inf, label=labelage, size=3, hjust=1, vjust=1)

printWithNumber(plotage)
printWithNumber(plotabeta)
printWithNumber(plotabeta42)
printWithNumber(plottau)
printWithNumber(plotptau)
}

dev.off()