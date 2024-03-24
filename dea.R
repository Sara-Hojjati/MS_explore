library('limma')

process_input_data <- function(input_file){
  data <- read.table(input_file, header=TRUE, check.names = FALSE)
  data <- as.matrix(data)
  proteins <- rownames(data)
  data <- apply(data, 2, as.numeric)
  rownames(data) <- proteins
  return(data)
}

make_design_matrix <- function(data, dataset="csf"){
  early_ms1a <- grep("Rep_MS_([1-9]$|[1][0-9]|2[0-1])", colnames(data)) #Replication cohort group A
  early_ms1b <- grep("Rep_MS_(2[2-9]|[3-4][0-9]|50|51)", colnames(data)) #Replication cohort group B
  early_ms5a <- grep("Dis_MS_([1-9]$|[1-3][0-9])", colnames(data)) #Discovery cohort group A
  early_ms5b <- grep("Dis_MS_([4-8][0-9])|9[0-2]", colnames(data)) #Discovery cohort group B
  control4 <- grep("Rep_HC", colnames(data))
  control6 <- grep("Dis_HC", colnames(data))
  
  factor_design <- c(rep(0, ncol(data)))
  factor_design[control6] <- 1
  factor_design[early_ms1a] <- 2
  factor_design[early_ms1b] <- 3
  factor_design[early_ms5a] <- 4
  factor_design[early_ms5b] <- 5
  
  design <- model.matrix(~0+factor(factor_design))
  if (dataset != "csf"){
    colnames(design) <- c("C4", "C6", "MS1b", "MS5a", "MS5b")
  } else {
    colnames(design) <- c("C4", "C6", "MS1a", "MS1b", "MS5a", "MS5b")
  }
  return(design) 
}

limma_analysis <- function(data, dataset="csf") {
  
  design <- make_design_matrix(data=data, dataset=dataset)
  
  if (dataset == "plasma") {
    contrast.matrix <- makeContrasts(MS1b-C4, (MS5a + MS5b)/2-C6, levels=design)
  } else if (dataset == "csf") {
    contrast.matrix <- makeContrasts((MS5a + MS5b)/2-C6, (MS1a + MS1b)/2-C4, levels=design)
  }
  
  fit <- lmFit(data, design)
  fit2 <- contrasts.fit(fit, contrasts = contrast.matrix, coefficients = NULL)
  fit2 <- eBayes(fit2, trend = FALSE)
  
  results <- decideTests(fit2, method="separate", adjust.method="BH", p.value=0.05, lfc=0)
  deps <- list()
  for (i in 1:ncol(contrast.matrix)){
    
    dep <- topTable(fit2, number=length(rownames(data)), coef=i, adjust="BH")
    compare_coef <- colnames(contrast.matrix)[i]
    print(paste(compare_coef, ': DEP(fdr) =', nrow(dep[dep$adj.P.Val < 0.05,]), ', DEP(nominal) =', nrow(dep[dep$P.Value < 0.05,])))

    dep <- dep[order(dep$P.Value),]
    
    deps[[i]] <- row.names(dep[dep$adj.P.Val < 0.05,])
  }
  return(deps)
}

#Read data
csf_data <- process_input_data(input_file="data/CSF_data.tsv")
plasma_data <- process_input_data(input_file="data/Plasma_data.tsv")

# Perform differential expression analysis
csf_deps <- limma_analysis(data=csf_data, dataset="csf")
plasma_deps <- limma_analysis(data=plasma_data, dataset="plasma")

# Write CSF Discovery DEPs to file 
write.table(csf_deps[[1]], "results/Discovery_DEP.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)