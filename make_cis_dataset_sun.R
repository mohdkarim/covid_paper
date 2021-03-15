# Install packages

list.of.packages <- c("googleCloudStorageR", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

# Extract dataset

proj <- "open-targets-genetics"

library(googleCloudStorageR)

gcs_auth("open-targets-genetics-fc5b6cda58e5.json")

buckets <- gcs_list_buckets(proj)

gcs_global_bucket("genetics-portal-analysis")

objects <- gcs_list_objects()

myobjects <- objects[grep("mohd", objects$name),]

colocds <- myobjects[grep("coloc-", myobjects$name),]

colocds2 <- colocds[-grep("TXNIP", colocds$name),]

cis_sun_coloc <- colocds[grep("cis-sun-coloc-formatted", colocds$name),]

cis_sun_coloc2 <- cis_sun_coloc[-grep("^mohd/cis-sun-coloc-formatted/$", cis_sun_coloc$name),]

# Loop to produce cis-Sun files for cis-GSMR analysis

sapply(seq_along(cis_sun_coloc2$name), function (x) {
  
  n1 <- strsplit(cis_sun_coloc2$name[x], "/")[[1]][3] # pQTL dataset
  
  require(googleCloudStorageR)
  require(dplyr)
  
  parsed_download <- gcs_get_object(cis_sun_coloc2$name[x],saveToDisk = n1, overwrite = T)
  
  df1 <- data.table::fread(n1, stringsAsFactors = F)
  
  unlink(n1)
  
  df2 <- df1 %>% select(SNP, A1, A2, freq, b, se, p, N) # tried to upload this to google bucket using gcs_upload but it does not work, so writing to vm instead.
  
  write.table(df2, paste0("sun_cis_gsmr_harmonised/", n1), row.names = F, quote = F)
  
  print('table written')
  
})

