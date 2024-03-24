library(dplyr)
library(verification)

##read the dataset in tsv format and store it in a dataframe called mydata
mydata<-read.table(file = 'data/CSF_corrected_data.tsv', sep = '\t', header = TRUE, row.names = 1)

#transpose the data,the dataset is arranged with each column being one variable and each row being one person
mydata<-t(mydata)
mydata<-data.frame(mydata)

#load the file with proteinnames in it and store it in a vector called proteinnames
load("data/proteinnames.RData")

#make a new dataset including only proteins that are Differentially expressed (DEPs=52)
mydata_proteins<-data.frame(mydata[,colnames(mydata)%in%proteinnames$proteinnames])

#make the protein data numeric
mydata_proteins<- data.frame(lapply(mydata_proteins, as.numeric), check.names = FALSE)

#store the additional data in a separate variable
mydata_additional <- data.frame(mydata[,c("sex","age","NEDA_EDA_2years","nARMSS","treatment_duration_index","ms_control")])

#make sex in the additional data a factor variable
mydata_additional[,c("sex")]<- as.numeric(mydata_additional[,c("sex")])
mydata_additional$sex<-factor(mydata_additional$sex,levels = c("0", "1"),labels = c("m", "f"))

#make NEDA_EDA_2years a numeic variable variable
mydata_additional$NEDA_EDA_2years<-factor(mydata_additional$NEDA_EDA_2years,levels = c("NEDA", "EDA"),labels = c("1", "0"))
mydata_additional[,c("NEDA_EDA_2years")]<- as.numeric(mydata_additional[,c("NEDA_EDA_2years")])

#store NEDA_EDA as a separate numeric variable
NEDA_EDA_2years_numeric<-mydata_additional$NEDA_EDA_2years#we store this numeric column for future use
NEDA_EDA_2years_numeric[which(NEDA_EDA_2years_numeric==1)]<-0
NEDA_EDA_2years_numeric[which(NEDA_EDA_2years_numeric==2)]<-1

#make NEDA_EDA-2years a factor variable
mydata_additional["NEDA_EDA_2years"]<-mydata["NEDA_EDA_2years"]
mydata_additional$NEDA_EDA_2years<-factor(mydata_additional$NEDA_EDA_2years,levels = c("NEDA", "EDA"),labels = c("NEDA", "EDA"))

#make age, treatment duration and nARMSS numeric
mydata_additional[,c("age","nARMSS","treatment_duration_index")]<- data.frame(lapply(mydata_additional[,c("age","nARMSS","treatment_duration_index")], as.numeric), check.names = FALSE)

#combine the two previous datasets to have a full dataset with all informations and proteins
mydata<-bind_cols(mydata_additional,mydata_proteins)

#devide the dataset into train (discovery cohort) and test (replication cohort) datasets
train_data<-mydata[c(1:92),]
train_data<-train_data[-c(35,73,75,76:86),]
test_data<-mydata[c(116:165),]

#storing variable "NEDA_EDA_2years" in two separate vectors for train dataset and test dataset
NEDA_EDA_2years_numeric_train<-NEDA_EDA_2years_numeric[1:92]
NEDA_EDA_2years_numeric_train<-NEDA_EDA_2years_numeric_train[-c(35,73,75,76:86)]
NEDA_EDA_2years_numeric_test<-NEDA_EDA_2years_numeric[116:165]

#storing the factor variable "NEDA_EDA_2years" in two separate vectors for train dataset and test dataset
NEDA_EDA_2years_factor_train<-mydata$NEDA_EDA_2years[1:92]
NEDA_EDA_2years_factor_train<-NEDA_EDA_2years_factor_train[-c(35,73,75,76:86)]
NEDA_EDA_2years_factor_test<-mydata$NEDA_EDA_2years[116:165]

#Building the logistic regression model for disease activity
model<-glm(NEDA_EDA_2years~NEFL+IL18+PDCD1+CD6,family=binomial,data = train_data)
trainresult<-roc.area(NEDA_EDA_2years_numeric_train,model$fitted.values)
AUC_train<-trainresult$A
pvalue_train<-trainresult$p.value

##testing the previously finalmodel in the test dataset and store the predictions in a dataframe called predicted
predicted<-predict(model,test_data,type=c("response"))
testresult<-roc.area(NEDA_EDA_2years_numeric_test,predicted)#AUC of the model performance in the test dataset
AUC_test<-testresult$A
pvalue_test<-testresult$p.value
print(model)#looking at the model
sprintf("%f is the AUC of the model for the Discovery cohort",AUC_train)
sprintf("%f is the pvalue of the model for the Discovery cohort",pvalue_train)
sprintf("%f is the AUC of the model for the Replication cohort",AUC_test)
sprintf("%f is the pvalue of the model for the Replication cohort",pvalue_test)