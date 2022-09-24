### Alvaro Benitez Mateo, 2021. ###

#Loading required libraries 
library(edgeR)
library(dplyr)

#Reading the gene count matrix that is going to be analyzed 
read_gene_count <- function(gene_count_matrix, isColumnName = TRUE, colName = "transcript_id", isFirstColumn = TRUE, colNum) {
  gene_matrix <- read.csv(gene_count_matrix)
  if (isColumnName == TRUE) {
  rownames(gene_matrix) <- gene_matrix$colName
  }
  if (isFirstColumn == TRUE) { 
  gene_matrix <- gene_matrix[, -c(1)]
  } else {
    gene_matrix <- gene_matrix[, -c(integer(colNum))]
  }
}

########################
gene_matrix <- read_gene_count(gene_count_matrix = "path/to/rnaseq_aba_transcript_count_matrix.csv", isColumnName = TRUE, colName = "transcript_id", isFirstColumn = TRUE)

#Calculate normalization factors and preprocess data
preprocessing <- function(gene_matrix) {
  d0 <- DGEList(gene_matrix)
  d0 <- calcNormFactors(d0)
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  processed_data <- d0[-drop,]
  if (dim(d0)[1] > dim(processed_data)[1] & dim(d0)[2] == dim(processed_data)[2]) {
  message("Dim ok :)")
  }
  return(processed_data)
}

########################
processed_data <- preprocessing(gene_matrix = gene_matrix)

########################
#Assigning sample groups
snames <- colnames(gene_matrix)
snames
group <- c("fh", "21", "21", "21", "22", "22", "22", "31", "31", "31", "32", "fh", "32", "32", "fh", "11", "11", "11", "12", "12", "12")

#Plot MDS to check sample clustering (It is not the same than a PCA)
plotMDS(processed_data, col =as.numeric(group), gene.selection = "pairwise") #MDS or PCoA

#Define model matrix to perform the differential expression analysis itself
mm <- model.matrix(~0 + group)

#Transform RNA-seq data ready for Linear Modelling
y <- voom(processed_data, mm, plot = T)

#Fit linear model for each gene
fit <- lmFit(y, mm)
head(coef(fit))
message("Setup completed")

#################################################################################################################################################
#It is time to perform the differential expression analysis between the selected groups named before (line 43). Select groups to be compared with: groupXX - groupYY.
#Data 1 can be use all in its own, ignore data 2 if no time zero or absolute reference available. 
#In case a global reference is available in the experimetal desing, it is recommended to compare two different conditions (e.g.: Treated / Non-Treated samples) vs time zero in data 1 and data 2.
#e.g.: DATA1 -> groupTreated - groupT0 // DATA2 -> groupNonTreated - groupT0

#DATA1
compare_data1 <- function(test, reference) {
  x <- c(paste(test, "-", reference, collapse = ""))
  contr <- makeContrasts(contrasts = x, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  return(tmp)
}

########################
contrast_data <- compare_data1(test = "group31", reference = "groupfh")


top.table1 <- topTable(tmp, sort.by = "P", n = Inf)
top.table1.p <- subset(top.table1, top.table1$adj.P.Val < 0.05)
top.table1.p.id <- top.table1.p
top.table1.p.id$id <- rownames(top.table1.p)

#Annotation
top.table1.p.id.annot <- top.table1.p.id
top.table1.p.id.annot$id <- substr(top.table1.p.id.annot$id, 1,nchar(top.table1.p.id.annot$id)-2)
top.table1.p.id <- top.table1.p.id.annot
cucurbita_pepo_gene_description <- read.delim("~/path/to/cucurbita_pepo_gene_description.txt", header=FALSE)
library(dplyr)
colnames(cucurbita_pepo_gene_description) <- c("id", "Annotation")
top.table1.p.id.annot.full <- dplyr::left_join(top.table1.p.id.annot, cucurbita_pepo_gene_description, by=c('id'='id'))
top.table1.p.id <- top.table1.p.id.annot.full
top.table1.p.id.over <- subset(top.table1.p.id, top.table1.p.id$logFC > 1.5)
top.table1.p.id.sub <- subset(top.table1.p.id, top.table1.p.id$logFC < -1.5)
top.table1.p.id.fc <- subset(top.table1.p.id, abs(top.table1.p.id$logFC) > 1.5)