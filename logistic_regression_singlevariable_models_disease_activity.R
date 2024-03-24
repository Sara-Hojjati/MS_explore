library(dplyr)
library(verification)

#read the dataset in tsv format and store it in a dataframe called mydata
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

#building empty dataframes that will be filled later
auc_test<-NULL
auc_all<-NULL
auc_test_all<-NULL
pvalueofauc_all<-NULL
pvalueofauc_test_all<-NULL
allproteinstrain<-train_data[,-c(1:6)]#storing all the proteins in a dataframe for train dataset
allproteinstest<-test_data[,-c(1:6)]#storing all the proteins in a dataframe for test dataset
k<-ncol(allproteinstrain)
p_value<-NULL

#creat a loop that
for(i in 1:k){
  actuals<-NEDA_EDA_2years_factor_train #store the dependent variable from train dataset
  protein<-allproteinstrain[i]
  newdata<-bind_cols(actuals,protein)#create a new dataset that only holds the dependent variable and one independent variable at a time
  colnames(newdata)[1]<-"actuals"
  colnames(newdata)[2]  <- "protein" 
  model <- glm(actuals~protein,family=binomial,data=newdata)
  actuals<-NEDA_EDA_2years_factor_test #store the dependent variable from test dataset
  protein<-allproteinstest[i]
  newdata<-bind_cols(actuals,protein) #create a new dataset that has proteins from test dataset
  colnames(newdata)[1]<-"actuals"
  colnames(newdata)[2]  <- "protein" 
  predicted<-predict(model,newdata,type=c("response")) #store the predictions on the test dataset in a dataframe called predicted
  auc<-roc.area(NEDA_EDA_2years_numeric_train,model$fitted.values) #the AUC from training the model using one independent variable at a time in the loop 
  auc_test<-roc.area(NEDA_EDA_2years_numeric_test,predicted) #the AUC from testing the model using one independent variable at a time in the loop
  auc_all[i]<-auc$A #store all the AUCs from training the models
  auc_test_all[i]<-auc_test$A #store all the AUCs from testing the models
  pvalueofauc_all[i]<-auc$p.value #store all the pvalues of the AUCs from training the models
  pvalueofauc_test_all[i]<-auc_test$p.value #store all the pvalues of the AUCs from testing the models
}

auc_pvalue_all<-bind_cols(auc_all,pvalueofauc_all,auc_test_all,pvalueofauc_test_all)#store all train and test AUCs and pvalues of AUCs in a dataframe
auc_pvalue_all<-data.frame(auc_pvalue_all)
auc_pvalue_all<-data.frame(auc_pvalue_all)
colnames(auc_pvalue_all)[1]<-"auc"
colnames(auc_pvalue_all)[2]<-"pvalue"
colnames(auc_pvalue_all)[3]<-"auc_test"
colnames(auc_pvalue_all)[4]<-"pvalue_test"
rownames(auc_pvalue_all)<-colnames(allproteinstrain)
proteinnames<-colnames(allproteinstrain)
proteinnames<-data.frame(proteinnames)
completeaucs<-bind_cols(auc_pvalue_all,proteinnames)
orderedaucs<-completeaucs %>% arrange(desc(auc)) #arranging the dataframe from lowest (most significant) pvalue in test models to highest
write.table(orderedaucs, "results/Disease_activity_model_AUCs.tsv", sep='\t', quote=FALSE) #look at the ordered AUCs and pvalues of AUCs of all single protein models
