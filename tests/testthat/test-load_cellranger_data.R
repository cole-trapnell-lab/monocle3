# Contributed by Nigel Delaney with adaptations

v3dname = "../testdata/cr3.0"
v2dname = "../testdata/cr2.0"



test_that("Cell Ranger 3.0 Filtered Data Loading", {
  cds = load_cellranger_data(v3dname)
  expect_equal(as.character(fData(cds)$gene_short_name[1]), "RP11-34P13.3")
  expect_equal(as.character(pData(cds)$barcode[1]), "AAAGTAGCACAGTCGC-1")
  expect_equal(dim(fData(cds)), c(3,2))
  expect_is(cds, "cell_data_set")
})

test_that("Cell Ranger 3.0 No Data Fails", {
  expect_error(load_cellranger_data(v3dname, barcode_filtered = FALSE))
})

test_that("Cell Ranger 2.0 Filtered Data Loading", {
  cds = load_cellranger_data(v2dname)
  expect_equal(as.character(fData(cds)$gene_short_name[1]), "MS4A1")
  expect_equal(as.character(pData(cds)$barcode[1]), "ATGCCAGAACGACT-1")
  expect_equal(dim(fData(cds)), c(240,2))
  expect_is(cds, "cell_data_set")
})

test_that("Cell Ranger 2.0 No Data Fails", {
  expect_error(load_cellranger_data(v3dname, barcode_filtered = FALSE))
})

# Now make sure raw data works
test_that("load_cellranger_matrix on 2.0 raw matrices", {
  # Copy the filtered data over and pretend it is unfiltered
  tmpdir = tempdir()
  odir = file.path(tmpdir, "outs", "raw_gene_bc_matrices", "hg19")
  odir2 = file.path(tmpdir, "outs", "filtered_gene_bc_matrices")
  dir.create(odir, recursive = TRUE)
  dir.create(odir2, recursive = TRUE)
  src = "../testdata/cr2.0/outs/filtered_gene_bc_matrices/hg19/"
  lapply(dir(src), function(p) file.copy(file.path(src, p), odir, recursive = TRUE))
  cds = load_cellranger_data(tmpdir, barcode_filtered=FALSE)
  unlink(odir, recursive = TRUE)
  unlink(odir2, recursive = TRUE)
  expect_equal(as.character(fData(cds)$gene_short_name[1]), "MS4A1")
  expect_equal(as.character(pData(cds)$barcode[1]), "ATGCCAGAACGACT-1")
  expect_equal(dim(fData(cds)), c(240,2))
  expect_is(cds, "cell_data_set")
})

test_that("load_cellranger_matrix on 3.0 raw matrices", {
  # Copy the filtered data over and pretend it is unfiltered
  tmpdir = tempdir()
  odir = file.path(tmpdir, "outs", "raw_feature_bc_matrix")
  odir2 = file.path(tmpdir, "outs", "filtered_feature_bc_matrix")
  dir.create(odir, recursive = TRUE)
  dir.create(odir2, recursive = TRUE)
  src = "../testdata/cr3.0/outs/filtered_feature_bc_matrix"
  lapply(dir(src), function(p) file.copy(file.path(src, p), odir, recursive = TRUE))
  cds = load_cellranger_data(tmpdir, barcode_filtered=FALSE)
  unlink(odir, recursive = TRUE)
  unlink(odir2, recursive = TRUE)
  expect_equal(as.character(fData(cds)$gene_short_name[1]), "RP11-34P13.3")
  expect_equal(as.character(pData(cds)$barcode[1]), "AAAGTAGCACAGTCGC-1")
  expect_equal(dim(fData(cds)), c(3,2))
  expect_is(cds, "cell_data_set")
})

test_that("load_cellranger_matrix on 3.0 with genome", {
  cds = load_cellranger_data(v3dname, genome="hg19")
  expect_equal(as.character(fData(cds)$gene_short_name[1]), "RP11-34P13.3")
  expect_equal(as.character(pData(cds)$barcode[1]), "AAAGTAGCACAGTCGC-1")
  expect_equal(dim(fData(cds)), c(3,2))
  expect_is(cds, "cell_data_set")
})

test_that("Cell Ranger 2.0 Filtered Data Loading with genome", {
  cds = load_cellranger_data(v2dname, genome="hg19")
  expect_equal(as.character(fData(cds)$gene_short_name[1]), "MS4A1")
  expect_equal(as.character(pData(cds)$barcode[1]), "ATGCCAGAACGACT-1")
  expect_equal(dim(fData(cds)), c(240,2))
  expect_is(cds, "cell_data_set")
})

test_that("Cell Ranger 2.0 Filtered Data Loading fake genome", {
  expect_error(load_cellranger_data(v2dname, genome="donkey"))
})
