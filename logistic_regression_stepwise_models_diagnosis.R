library(dplyr)
library(verification)

#read the dataset in tsv format and store it in a dataframe called mydata
mydata<-read.table(file = 'data/CSF_corrected_data.tsv', sep = '\t', header = TRUE, row.names = 1)

#transpose the data,the dataset is arranged with each column being one variable and each row being one person
mydata<-t(mydata)

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

#make NEDA_EDA_2years a factor variable
mydata_additional$NEDA_EDA_2years<-factor(mydata_additional$NEDA_EDA_2years,levels = c("NEDA", "EDA"),labels = c("NEDA", "EDA"))

#store ms_control as a separate numeric variable
ms_control_numeric<-as.numeric(mydata_additional[,c("ms_control")])#we store this numeric column for future use

#make ms_control a factor variable
mydata_additional[,c("ms_control")]<- as.numeric(mydata_additional[,c("ms_control")])
mydata_additional$ms_control<-factor(mydata_additional$ms_control,levels = c("1", "0"),labels = c("MS", "HC"))

#make the numeric variable a correct 0 and 1 varibale
ms_control_numeric[which(ms_control_numeric==0)]<-3
ms_control_numeric[which(ms_control_numeric==1)]<-0
ms_control_numeric[which(ms_control_numeric==3)]<-1

#make age, treatment duration and nARMSS numeric
mydata_additional[,c("age","nARMSS","treatment_duration_index")]<- data.frame(lapply(mydata_additional[,c("age","nARMSS","treatment_duration_index")], as.numeric), check.names = FALSE)

#combine the two previous datasets to have a full dataset with all informations and proteins
mydata<-bind_cols(mydata_additional,mydata_proteins)

#devide the dataset into train (discovery cohort) and test (replication cohort) datasets
train_data<-mydata[c(1:115),]
test_data<-mydata[c(116:nrow(mydata)),]

#storing variable "ms_control_numeric" in two separate vectors for train dataset and test dataset
ms_control_numeric_train<-ms_control_numeric[1:115]
ms_control_numeric_test<-ms_control_numeric[116:length(ms_control_numeric)]

#storing the factor variable "ms_control" in two separate vectors for train dataset and test dataset
ms_control_factor_train<-train_data$ms_control
ms_control_factor_test<-test_data$ms_control

#forward selection followed by backward selection
nullmodel<-glm(ms_control~1,family=binomial,data = train_data)#null model is a model based on no variable
fullmodel<-glm(ms_control~.-NEDA_EDA_2years-nARMSS-treatment_duration_index,family=binomial,data = train_data)#full model is a model based on all variables
scope = list(lower=formula(nullmodel),upper=formula(fullmodel))#scop is the a range from nullmodel to fullmodel
forwardmodel<-step(nullmodel,scope,direction = "forward")#performing the initial forwardmodel
pvalues<-summary(forwardmodel)$coefficients[,4]#extracting the p-values of each variable from the model

#excluding the variable with highest insignificant pvalue
if (max(pvalues)>0.05)
{
  pvalues<-pvalues[pvalues != max(pvalues)]
}

#storing the forwardmodel variables in a new dataset
newdata <- forwardmodel$model[,names(forwardmodel$model) %in% names(pvalues)]

#constructing a backward model based on the forwardmodel variables
backwardmodel<-glm(ms_control_factor_train~.,family=binomial,data = newdata)

#a loop that checks for the most insignificant variable each time and excludes that until all variables become significant
p<-10
for (j in 1:p) {
  
  pvalues<-summary(backwardmodel)$coefficients[,4]
  pvaluesnoint<-pvalues[-1];
  if (max(pvaluesnoint)>0.05){
    pvaluesnoint<-pvaluesnoint[pvaluesnoint != max(pvaluesnoint)]
    newdata <- backwardmodel$model[,names(backwardmodel$model) %in% names(pvaluesnoint)]
    backwardmodel<-glm(ms_control_factor_train~.,family=binomial,data = newdata)
  }
}
finalmodel<-backwardmodel #this is the final model
print(summary(finalmodel))#looking at a summary of the final model

#model evaluation: AUC and p-value of the AUC
trainresult<-roc.area(ms_control_numeric_train,finalmodel$fitted.values)
AUC_train<-trainresult$A
pvalue_train<-trainresult$p.value

#testing the previously finalmodel in the test dataset and store the predictions in a dataframe called predicted
predicted<-predict(finalmodel,test_data,type=c("response"))
testresult<-roc.area(ms_control_numeric_test,predicted)#AUC of the model performance in the test dataset
AUC_test<-testresult$A
pvalue_test<-testresult$p.value
finalmodel#looking at the model
sprintf("%f is the AUC of the model for the Discovery cohort",AUC_train)
sprintf("%f is the pvalue of the model for the Discovery cohort",pvalue_train)
sprintf("%f is the AUC of the model for the Replication cohort",AUC_test)
sprintf("%f is the pvalue of the model for the Replication cohort",pvalue_test)
