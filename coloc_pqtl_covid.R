#!/usr/bin/Rscript

# Install packages if not installed already (will need to install all from source)

# MR > 'plotly' ( > httr ('curl', openssl)), 'iterpc' (> gmp)

list.of.packages <- c("googleCloudStorageR", "data.table", "remotes", "MendelianRandomization")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

# Install coloc separately

# library(remotes, lib.loc = "Rlib")
# withr::with_libpaths(new = "Rlib", install_github("chr1swallace/coloc", dependencies = TRUE))

if(!require("remotes"))
  install.packages("remotes") # if necessary
library(remotes)
list.of.packages <- "coloc"
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#
if(length(new.packages)) install_github("chr1swallace/coloc")

# Import one sun protein coloc file and one covid coloc file from google storage

setwd("/home/mk31/")

library(googleCloudStorageR)

gcs_auth("/Users/mk31/open-targets-genetics-fc5b6cda58e5.json")

proj <- "open-targets-genetics"

gcs_global_bucket("genetics-portal-analysis")

buckets <- gcs_list_buckets(proj)

objects <- gcs_list_objects()

myobjects <- objects[grep("mohd", objects$name),]

colocds <- myobjects[grep("coloc-", myobjects$name),]

# remove TYRO for now (it does not appear to have any SNPs beyond chr 9, need to let James know) + there appears to be only 76 cis-SNPs for TXNIP but only 1 of THEM HAVE AFs from gnomad v2.1.1 
# TYRO included back as James corrected this

# colocds2 <- colocds[-grep("TYRO3|TXNIP", colocds$name),]

colocds2 <- colocds[-grep("TXNIP", colocds$name),]

# creating separate Sun and COVID coloc dataset URLs

cis_sun_coloc <- colocds2[grep("cis-sun-coloc-formatted", colocds2$name),]

# ensure sun proteins are repeated thrice to run coloc with each A2/B2/C2 covid datasets

cis_sun_coloc2 <- do.call("rbind", replicate(7, cis_sun_coloc, simplify = FALSE))

cis_sun_coloc3 <- cis_sun_coloc2[-grep("^mohd/cis-sun-coloc-formatted/$", cis_sun_coloc2$name),]

# ensure covid genes match with sun proteins

cis_covid_coloc <- colocds2[grep("covid4-oct2020-coloc-formatted", colocds2$name),]

cis_covid_coloc2 <- cis_covid_coloc[-grep("^mohd/covid4-oct2020-coloc-formatted/$", cis_covid_coloc$name),]

# sort coloc sun and covid

cis_covid_coloc2 <- cis_covid_coloc2[order(cis_covid_coloc2$name),]

cis_sun_coloc3 <- cis_sun_coloc3[order(cis_sun_coloc3$name),]

# both traits should have the same number of rows

colocdf <- lapply(1:length(cis_sun_coloc3$name), function (x) {
  
  # load libraries
  
  require(googleCloudStorageR)
  require(data.table)
  require(coloc)
  require(MendelianRandomization)
  
  # Read datasets
  
  n1 <- strsplit(cis_sun_coloc3$name[x], "/")[[1]][3] # pQTL dataset
  
  n2 <- strsplit(cis_covid_coloc2$name[x], "/")[[1]][3] # trait dataset
  
  print(n1)
  print(n2)
  
  parsed_download <- gcs_get_object(cis_sun_coloc3$name[x],saveToDisk = n1, overwrite = T)
  
  parsed_download <- gcs_get_object(cis_covid_coloc2$name[x],saveToDisk = n2, overwrite = T)
  
  df1 <- data.table::fread(n1, stringsAsFactors = F) 
  
  df2 <- data.table::fread(n2, stringsAsFactors = F) 
  
  # remove SNPs that have beta with missing values
  
  df1 <- df1[which(!is.na(df1$b)),]
  
  df2 <- df2[which(!is.na(df2$b)),]
  
  # merge datasets
  
  df <- merge(df1, df2, by = "SNP")
  
  # coloc analysis
  
  ## create function to convert eaf to maf and apply
  
  eaf_to_maf = function(eaf) {
    min(eaf, 1-eaf)
  }
  
  df$maf.x <- sapply(df$freq.x, eaf_to_maf)
  
  df$maf.y <- sapply(df$freq.y, eaf_to_maf)
  
  # remove SNPs with NA mafs
  
  df <- df[which(!is.na(df$maf.x)),]
  
  ## Make coloc dataset (left: pqtl)
  
  left_n = df$N[1]
  left_type = 'quant'
  left_data = with(df, list(
    pvalues=p.x,
    N=left_n,
    beta = b.x,
    varbeta = se.x^2,
    MAF=maf.x,
    type=left_type,
    snp = SNP
  ))
  
  # Make coloc dataset (right: trait)
  
  right_n = df$N_total[1]
  # right_ncases = df$N_cases[1]
  # right_prop = right_ncases / right_n
  right_type = 'cc'
  right_data = with(df, list(
    pvalues=p.y,
    N=right_n,
    beta = b.y,
    varbeta = se.y^2,
    MAF=maf.y,
    type=right_type,
    # s=right_prop,
    s = 0.5, # this is required for cc dataset 2 but it does not seem to affected coloc estimates, so having this 0.5 for now
    snp = SNP
  ))
  
  # Run coloc
  
  coloc_res <- coloc.abf(left_data, right_data)
  
  print("coloc complete")
  
  # Extract coloc summary results
  
  s <- data.frame(coloc_res$summary)
  s2 <- data.frame(t(s), stringsAsFactors = F)
  row.names(s2) <- NULL
  SNP <- coloc_res$results$snp[which.max(coloc_res$results$SNP.PP.H4)]
  s2$coloc_snp <- SNP # annotate with SNP with highest PP4
  s2$left_trait <- gsub(".csv", "_sun_pQTL", n1)
  s2$right_trait <- gsub(".csv", "", n2)
  s2$rsid <- df$rsid.x[which(df$SNP %in% s2$coloc_snp)]
  print(s2$rsid)
  
  # perform single-SNP MR with single coloc_snp
  
  mr1 <- df[which(df$SNP %in% s2$coloc_snp),]
  print(SNP)
  print(mr1$SNP)
  mr2 <- mr_ivw(mr_input(bx = mr1$b.x, bxse = mr1$se.x, by = mr1$b.y, byse = mr1$se.y))
  s2$mr_beta <- mr2@Estimate
  s2$mr_se <- mr2@StdError
  s2$mr_lci <- mr2@CILower
  s2$mr_uci <- mr2@CIUpper
  s2$mr_p <- mr2@Pvalue
  
  
  # upload s2 to a google bucket
  
  # f <- function(input, output) write.csv(input, row.names = FALSE, file = output)
  # 
  # gcs_upload(s2,
  #            object_function = f,
  #            name = "eg.csv")
  
  unlink(c(n1, n2))
  
  s2
  
})

colocdf2 <- do.call(rbind, colocdf)

write.csv(colocdf2, "coloc_sun_covid.csv", row.names = F, quote = F)