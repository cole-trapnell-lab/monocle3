context("test-estimate_dispersions")

test_that("estimate_dispersions() properly validates its input",{
  lung <- load_lung()
  lung <- detect_genes(lung)
  lung <- estimate_size_factors(lung)
  lung <- estimate_dispersions(lung)
  lung_pData <- pData(lung)
  expect_equal(colnames(lung_pData), c("file", "total_mass", "internal_scale",
                                       "external_scale",
                                       "median_transcript_frags", "BioSample",
                                       "age", "genotype", "Sample.Name",
                                       "SRA.Sample", "MBases", "MBytes",
                                       "SRA.Study", "BioProject",
                                       "source_name", "strain", "tissue",
                                       "Assay.Type", "Center.Name",
                                       "Platform", "Consent", "Time",
                                       "Size_Factor", "Total_mRNAs",
                                       "endogenous_RNA", "Pseudotime", "State",
                                       "Parent", "num_genes_expressed"))
  expect_true(all(lung_pData$total_mass >= 392))
  expect_true(all(lung_pData$total_mass <= 4086580))
  expect_true(all(lung_pData$internal_scale >= 0.00153031))
  expect_true(all(lung_pData$internal_scale <= 11.0584))
  expect_true(all(lung_pData$external_scale == 1))
  expect_true(all(lung_pData$median_transcript_frags >= 1.03129e-06))
  expect_true(all(lung_pData$median_transcript_frags <= 3.48372))
  expect_true(all(substring(lung_pData$BioSample, 0, 5) == "SAMN0"))
  expect_equal(levels(lung_pData$age),
               c("Embryonic day 14.5", "Embryonic day 16.5",
                 "Embryonic day 18.5", "post natal day 107"))
  expect_equal(levels(lung_pData$genotype),
               c("Sftpc-Cre-ERT2-rtta -/- tetO-HIST1H2BJ-GFP+/-)",
                 "wild type"))
  expect_true(all(substring(lung_pData$Sample.Name, 0, 6) == "GSM127"))
  expect_equal(substring(lung_pData$Sample.Name[1], 7, 10), "1863")
  expect_equal(substring(lung_pData$Sample.Name[185], 7, 10), "2062")
  expect_equal(toString(lung_pData$SRA.Sample[1]), "SRS504890")
  expect_equal(toString(lung_pData$SRA.Sample[185]), "SRS505090")
  expect_true(all(lung_pData$MBases >= 0))
  expect_true(all(lung_pData$MBases <= 961))
  expect_true(all(lung_pData$MBytes >= 0))
  expect_true(all(lung_pData$MBases <= 1000))
  expect_equal(levels(lung_pData$SRA.Study), "SRP033209")
  expect_equal(levels(lung_pData$BioProject), "GSE52583")
  expect_equal(levels(lung_pData$source_name), "distal lung epithelium")
  expect_equal(levels(lung_pData$strain), "C57BL/6J")
  expect_equal(levels(lung_pData$tissue), "lung")
  expect_equal(levels(lung_pData$Assay.Type), "RNA-Seq")
  expect_equal(levels(lung_pData$Center.Name), "GEO")
  expect_equal(levels(lung_pData$Platform), "ILLUMINA")
  expect_equal(levels(lung_pData$Consent), "public")
  expect_true(all(lung_pData$Size_Factor >= 0.01582144)) #0.2896533
  expect_true(all(lung_pData$Size_Factor <= 5.020053)) #4.08098
  expect_true(all(lung_pData$Total_mRNAs >= 13177.56)) #50.43858
  expect_true(all(lung_pData$Total_mRNAs <= 127126.5)) #16184.72
  expect_true(all(lung_pData$endogenous_RNA >= 251.1154))
  expect_true(all(lung_pData$endogenous_RNA <= 194030))
  expect_true(all(lung_pData$Pseudotime >= 0))
  expect_true(all(lung_pData$Pseudotime <= 28.45905)) #16.67064
  expect_equal(levels(lung_pData$State), c("1", "2", "3"))
  expect_true((is.na(substring(lung_pData$Parent, 0, 6)) ||
                 substring(lung_pData$Parent, 0, 6) == "SRR103"))
  expect_true(all(lung_pData$num_genes_expressed >= 10))
  expect_true(all(lung_pData$num_genes_expressed <= 196))
  expect_equal(nrow(lung@disp_fit_info[["blind"]]$disp_table), 213)
  expect_equal(ncol(lung@disp_fit_info[["blind"]]$disp_table), 3)
  expect_equal(as.numeric(attr(lung@disp_fit_info[["blind"]]$disp_func,
                               "coefficients")[1]), 2.03474, tolerance = 1e-5)
  expect_equal(as.numeric(attr(lung@disp_fit_info[["blind"]]$disp_func,
                               "coefficients")[2]), 8.888378, tolerance = 1e-5)
})
