# this R script requires bigrquery installed

list.of.packages <- c("bigrquery", "data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

sun2 <- data.table::fread("/Users/mk31/sun_proteins_ensid_coloc.csv", stringsAsFactors = F)

sapply(seq_along(sun2$ensembl_id), function (y) {
  
  print(sun2$hgnc[y])
  
  chr <- sun2$chr[y]
  r1 <- sun2$r1[y]
  r2 <- sun2$r2[y]
  
  # studies <- c("COVID19_HGI_A1_ALL_20201020_N", "COVID19_HGI_A2_ALL_leave_23andme_20201020_N", "COVID19_HGI_B1_ALL_20201020_N", "COVID19_HGI_B2_ALL_leave_23andme_20201020_N", "COVID19_HGI_C1_ALL_leave_23andme_20201020_N", "COVID19_HGI_C2_ALL_leave_23andme_20201020_N", "COVID19_HGI_D1_ALL_20201020_N")
  
  # studies <- "COVID19_HGI_B2_ALL_leave_23andme_20201020_N"
  
  # studies <- c("COVID19_HGI_A2_ALL_leave_23andme_20201020_N", "COVID19_HGI_C2_ALL_leave_23andme_20201020_N")
  
  # studies <- c("COVID19_HGI_A1_ALL_20201020_N", "COVID19_HGI_B1_ALL_20201020_N", "COVID19_HGI_C1_ALL_leave_23andme_20201020_N", "COVID19_HGI_D1_ALL_20201020_N")
  
  # n_cases <- c(269, 4336, 2430, 7885, 11085, 17965, 3204)
  
  # n_cases <- 7885
  
  # n_cases <- c(4336, 17965)
  
  # n_cases <- c(269, 2430, 11085, 3204)
  
  studies <- c("COVID_A1_h", "COVID_A2_h", "COVID_B1_h", "COVID_B2_h", "COVID_C1_h", "COVID_C2_h", "COVID_D1_h")
  
  n_cases <- c(269, 4336, 2430, 7885, 11085, 17965, 3204)
  
  # FORMAT OUTCOME TABLES FOR GSMR AND SAVE IN MY GOOGLE BUCKET
  
  sapply(seq_along(studies), function (x) {
    
    # build query
    
    require("bigrquery")
    
    project <- "open-targets-genetics"
    
    study <- studies[x]
    n_case <- n_cases[x]
    # chr <- sun2$chr[x]
    # r1 <- sun2$r1[x]
    # r2 <- sun2$r2[x]
    
    n_study <- paste0(sun2$ensembl_id[y], "_", sun2$SOMAMER_ID[y], "_", study)
    
    print(study)
    
    # THIS QUERY WILL REMOVE ROWS WITH DUPLICATE SNP IDs
    
    qry <- paste0("WITH 
  BIG_QUERY AS (
WITH
  BIG_SUB_QUERY AS (
  SELECT
    CONCAT(string_field_2, ':', string_field_3, ':', string_field_4, ':', string_field_5) AS SNP,
    string_field_1 as rsid,
    safe_cast(string_field_2 as INT64) AS chr,
    safe_cast(string_field_3 as INT64) AS pos,
    string_field_5 AS A1,
    string_field_4 AS A2,
    string_field_10 AS freq,
    string_field_6 AS b,
    string_field_17 AS se,
    string_field_18 AS p,
    string_field_19 AS N_total,",
                  n_case, " AS N_cases
    
  FROM
    `mohd_hypothesis.", study, "`", ")",
                  "SELECT
  *
FROM (
  SELECT
    *,
    ROW_NUMBER() OVER (PARTITION BY SNP) row_number
  FROM
    BIG_SUB_QUERY)
WHERE
  row_number = 1)
  SELECT
    SNP,
    rsid,
    chr,
    pos,
    A1,
    A2,
    freq,
    b,
    se,
    p,
    N_total,
    N_cases
   FROM
    BIG_QUERY
  WHERE
  chr =", chr, " AND pos BETWEEN ", r1, " AND ", r2)
    
    # run bigquery
    
    tb <- bigrquery::bq_project_query(project, qry)
    
    print("done")
    
    # save in my gcloud bucket
    
    bigrquery::bq_perform_extract(tb, paste0("gs://genetics-portal-analysis/mohd/covid4-oct2020-coloc-formatted/", n_study, ".csv"), destination_format = "CSV", print_header = TRUE)
    
    print(paste0(study, " coloc formatted and uploaded to GC bucket"))
    
  })
  
})