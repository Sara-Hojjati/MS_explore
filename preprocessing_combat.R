library(ChAMP)
library(ggplot2)

source('champ_svd_homemade.R')

process_data <- function(input_data, clinical_file) {
  exp_data <- read.table(input_data, row.names=1, check.names = FALSE)
  
  clinical_data <- read.table(clinical_file, sep='\t', header=TRUE)
  
  rownames(clinical_data) <- clinical_data$sample_ID
  clinical_data <- clinical_data[, c('sex', 'age')]
  clinical_data[, "sex"][clinical_data[, "sex"] == "m"] = 0
  clinical_data[, "sex"][clinical_data[, "sex"] == "f"] = 1
  clinical_data[, "sex"] <- as.numeric(clinical_data[, "sex"])
  clinical_data[, "group"] <- 0
  clinical_data[, "group"][grep("MS", row.names(clinical_data))] <- 1
  clinical_data[, "site"] <- 0
  clinical_data[, "site"][grep("Dis", row.names(clinical_data))] <- 1
  
  sample_key <- clinical_data[colnames(exp_data),]
  
  data <- list()
  data$data <- data.frame(lapply(exp_data, as.numeric), check.names = FALSE)
  rownames(data$data) <- rownames(exp_data)
  data$sample_key <- sample_key
  
  return(data)
}

run_combat <- function(data){
  expdata_combat <- champ.runCombat(beta = data$data, 
                                    pd = data$sample_key, 
                                    variablename = "group",
                                    logitTrans = FALSE,
                                    batchname = c("site"))
  
  data_combat <- list()
  data_combat$data <- expdata_combat
  data_combat$sample_key <- data$sample_key
  
  return(data_combat)
}

# ComBat correction on CSF data
data <- process_data(input_data = "data/CSF_control_corrected.tsv",
                     clinical_file = "data/clinical_data.csv")

data_combat <- run_combat(data)

write.table(data_combat$data, "data/CSF_batch_corrected.tsv", sep="\t", quote=FALSE)

# ComBat correction on plasma data
data <- process_data(input_data="data/Plasma_control_corrected.tsv",
                     clinical_file = "data/clinical_data.csv")

data_combat <- run_combat(data)

write.table(data_combat$data, "data/Plasma_batch_corrected.tsv", sep="\t", quote=FALSE)
