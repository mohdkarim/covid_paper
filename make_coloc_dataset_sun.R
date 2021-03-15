
sun2 <- data.table::fread("sun_proteins_ensid_coloc.csv", stringsAsFactors = F)

sapply(seq_along(sun2$sun_id), function (y) {
  
  print(sun2$sun_id[y])
  
  chr <- sun2$chr[y]
  r1 <- sun2$r1[y]
  r2 <- sun2$r2[y]
  study <- sun2$sun_id[y]
    
    # build query
    
    require("bigrquery")
    
    project <- "open-targets-genetics"
    
    n_study <- paste0(sun2$ensembl_id[y], "_", sun2$SOMAMER_ID[y])
    
    print(study)
    
    # THIS QUERY WILL ADD IN AF AND REMOVE ROWS WITH DUPLICATE SNP IDs
    
    qry <- paste0("WITH 
  BIG_QUERY AS (
 WITH
  BIG_SUB_QUERY AS (
  SELECT
  CONCAT(string_field_3, ':', string_field_4, ':', string_field_5, ':', string_field_6) AS SNP,
  string_field_2 AS rsid,
  safe_cast(string_field_3 as INT64) AS chr,
  safe_cast(string_field_4 as INT64) AS pos,
  string_field_6 AS A1,
  string_field_5 AS A2,
  g.af.gnomad_nfe AS freq,
  string_field_7 AS b,
  string_field_19 AS se,
  string_field_21 AS p,
  '3300' AS N
FROM
  `mohd_hypothesis.",study, "` m
INNER JOIN
  `200201.variants` g
ON
  CONCAT(m.string_field_3, ':', m.string_field_4, ':', m.string_field_5, ':', m.string_field_6) = CONCAT(g.chr_id, ':', g.position, ':', g.ref_allele, ':', g.alt_allele))
SELECT
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
    N
   FROM
    BIG_QUERY
  WHERE
  chr =", chr, " AND pos BETWEEN ", r1, " AND ", r2)
    
    # run bigquery
    
    tb <- bigrquery::bq_project_query(project, qry)
    
    print("done")
    
    # save in my gcloud bucket
    
    bigrquery::bq_perform_extract(tb, paste0("gs://genetics-portal-analysis/mohd/cis-sun-coloc-formatted/", n_study, ".csv"), destination_format = "CSV", print_header = TRUE)
    
    print(paste0(study, " coloc formatted and uploaded to GC bucket"))
  
})