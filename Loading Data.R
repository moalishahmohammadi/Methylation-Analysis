setwd("Your Working Directory")

#### Loading required packages
library(SummarizedExperiment)
library(TCGAbiolinks)
library(readxl)
library(DT)

#### Characterizing and downloading the desired data

query.exp <- GDCquery(project = "TCGA-CESC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts",
                      legacy = FALSE)

query.met <- GDCquery(project = "TCGA-CESC",
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450",
                      data.format = "TXT",
                      legacy = FALSE)

GDCdownload(query = query.exp , method = "api" , files.per.chunk = 5 )
GDCdownload(query = query.met , method = "api" , files.per.chunk = 5 )

exp <- GDCprepare(query.exp , summarizedExperiment = T, )

met <- GDCprepare(query.met , summarizedExperiment = T, )

