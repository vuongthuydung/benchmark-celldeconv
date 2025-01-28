# QC
## First: cells with library size, mitochondrial or ribosomal content further than three MAD away were discarded
filterCells <- function(filterParam){
	cellsToRemove <- which(filterParam > median(filterParam) + 3 * mad(filterParam) | filterParam < median(filterParam) - 3 * mad(filterParam) )
	cellsToRemove
}
libSizes <- colSums(data)
gene_names <- rownames(data)

mtID <- grepl("^MT-|_MT-", gene_names, ignore.case = TRUE)
rbID <- grepl("^RPL|^RPS|_RPL|_RPS", gene_names, ignore.case = TRUE)

mtPercent <- colSums(data[mtID, ])/libSizes
rbPercent <- colSums(data[rbID, ])/libSizes

lapply(list(libSizes = libSizes, mtPercent = mtPercent, rbPercent = rbPercent), filterCells) %>% 
	unlist() %>% 
	unique() -> cellsToRemove

if(length(cellsToRemove) != 0){
	data <- data[,-cellsToRemove]
	full_phenoData <- full_phenoData[-cellsToRemove,]
}

## Keep only "detectable" genes: at least 5% of cells (regardless of the group) have a read/UMI count different from 0
keep <- which(Matrix::rowSums(data > 0) >= round(0.05 * ncol(data)))
data = data[keep,]

# Data split into training/test
# set.seed(24)

original_cell_names = colnames(data) # cellID
# updates the column names of the data frame to be the corresponding cell types from the full_phenoData frame based on matching cell IDs
colnames(data) <- as.character(full_phenoData$cellType[match(colnames(data),full_phenoData$cellID)])
rnames <- rownames(data)

# Keep CTs with >= 50 cells after QC
cell_counts = table(colnames(data)) # the frequency of each unique cell type
to_keep = names(cell_counts)[cell_counts >= 50]
pData <- full_phenoData[full_phenoData$cellType %in% to_keep,]
to_keep = which(colnames(data) %in% to_keep)   
# data <- data[,to_keep]
data <- data.frame(as.list(data)[to_keep], check.names=FALSE)
rownames(data) <- rnames
original_cell_names <- original_cell_names[to_keep]

# Define training and test samples
train_samples <- c("human1","human2")
test_samples <- c("human3")

# Filter cells based on sample IDs for training and testing
train_indices <- which(pData$sampleID %in% train_samples)
test_indices <- which(pData$sampleID %in% test_samples)

# Create training and test data
train_data <- data[, train_indices]
test_data <- data[, test_indices]

# Update phenodata
train_phenoData <- pData[train_indices, ]
test_phenoData <- pData[test_indices, ]

# Update column names of train and test data
colnames(train_data) <- original_cell_names[train_indices]
colnames(test_data) <- original_cell_names[test_indices]

# Generate pseudo-bulk mixtures (T) on test data
generator <- Generator(sce = test_data, phenoData = test_phenoData, Num.mixtures = number_mixtures, pool.size = number_cells)
T <- generator[["T"]]
P <- generator[["P"]]

# Reference matrix C
C <- train_data

# Update phenodata for reference matrix C
phenoDataC <- train_phenoData
 if(length(grep("[N-n]ame",colnames(phenoDataC))) > 0){
        sample_column = grep("[N-n]ame",colnames(phenoDataC))
} else {
        sample_column = grep("[S-s]ample|[S-s]ubject",colnames(phenoDataC))
}

colnames(phenoDataC)[sample_column] = "SubjectName"
rownames(phenoDataC) = phenoDataC$cellID

# find common genes in bulk and scdata
keep = intersect(rownames(C),rownames(T)) 
C = C[keep,]
T = T[keep,]

require(xbioc)
C.eset <- Biobase::ExpressionSet(assayData = as.matrix(C),phenoData = Biobase::AnnotatedDataFrame(phenoDataC)) # single-cell data
T.eset <- Biobase::ExpressionSet(assayData = as.matrix(T)) # pseudo-bulk data

# Start deconvolving
# DWLS
ptm <- proc.time()
STRING=NULL
path=paste(getwd(),"/results_",STRING,sep="")

if(! dir.exists(path)){ #to avoid repeating marker_selection step when removing cell types; Sig.RData automatically created

    dir.create(path)
    Signature <- DWLS::buildSignatureMatrixMAST(scdata = C, id = as.character(phenoDataC$cellType), path = path, diff.cutoff = 0.5, pval.cutoff = 0.01)
} else {
	#re-load signature

    load(paste(path,"Sig.RData",sep="/"))
    Signature <- Sig
}
        
RESULTS_DWLS <- apply(T,2, function(x){
        b = setNames(x, rownames(T))
        tr <- DWLS::trimData(Signature, b)
        RES <- t(solveDampenedWLS(tr$sig, tr$bulk))
    })

rownames(RESULTS_DWLS) <- as.character(unique(phenoDataC$cellType))
RESULTS_DWLS = apply(RESULTS_DWLS,2,function(x) ifelse(x < 0, 0, x)) #explicit non-negativity constraint
RESULTS_DWLS = apply(RESULTS_DWLS,2,function(x) x/sum(x)) #explicit STO constraint
time_DWLS = proc.time() - ptm
RESULTS_DWLS = RESULTS_DWLS[gtools::mixedsort(rownames(RESULTS_DWLS)),]
RESULTS_DWLS = data.table::melt(RESULTS_DWLS)
colnames(RESULTS_DWLS) <-c("CT","tissue","observed_values")

# Music
sce <- as(as(C.eset, "SummarizedExperiment"), "SingleCellExperiment")
names(assays(sce))=c("counts")
ptm <- proc.time()
RESULTS_Music = t(MuSiC::music_prop(bulk.mtx = exprs(T.eset), sc.sce = sce, clusters = 'cellType',
                            markers = NULL, normalize = FALSE, samples = 'SubjectName', 
                            verbose = FALSE)$Est.prop.weighted)
time_Music = proc.time() - ptm
RESULTS_Music = RESULTS_Music[gtools::mixedsort(rownames(RESULTS_Music)),]
RESULTS_Music = data.table::melt(RESULTS_Music)
colnames(RESULTS_Music) <-c("CT","tissue","observed_values")

# SCDC
ptm <- proc.time()
RESULTS_scdc <- t(SCDC::SCDC_prop(bulk.eset = T.eset, sc.eset = C.eset, ct.varname = "cellType", sample = "SubjectName", ct.sub = unique(as.character(phenoDataC$cellType)), iter.max = 1000, verbose=FALSE)$prop.est.mvw)
time_scdc = proc.time() - ptm
RESULTS_scdc = RESULTS_scdc[gtools::mixedsort(rownames(RESULTS_scdc)),]
RESULTS_scdc = data.table::melt(RESULTS_scdc)
colnames(RESULTS_scdc) <-c("CT","tissue","observed_values")

# Combined1
ptm <- proc.time()
RESULTS_combined1 <- ReferenceBasedDecompositionNew1(T.eset, C.eset, markers=rownames(Signature), use.overlap=FALSE, verbose=FALSE)$bulk.props
time_combined1 = proc.time() - ptm
RESULTS_combined1 = RESULTS_combined1[gtools::mixedsort(rownames(RESULTS_combined1)),]
RESULTS_combined1 = data.table::melt(RESULTS_combined1)
colnames(RESULTS_combined1) <-c("CT","tissue","observed_values")

# Combined2
ptm <- proc.time()
RESULTS_combined2 <- ReferenceBasedDecompositionNew2(T.eset, C.eset, markers=NULL, use.overlap=FALSE, verbose=FALSE)$bulk.props
time_combined2 = proc.time() - ptm
RESULTS_combined2 = RESULTS_combined2[gtools::mixedsort(rownames(RESULTS_combined2)),]
RESULTS_combined2 = data.table::melt(RESULTS_combined2)
colnames(RESULTS_combined2) <-c("CT","tissue","observed_values")

#Bisque
ptm <- proc.time()
RESULTS_Bisque <- BisqueRNA::ReferenceBasedDecomposition(T.eset, C.eset, markers=NULL, use.overlap=FALSE, verbose=FALSE)$bulk.props
time_Bisque = proc.time() - ptm
RESULTS_Bisque = RESULTS_Bisque[gtools::mixedsort(rownames(RESULTS_Bisque)),]
RESULTS_Bisque = data.table::melt(RESULTS_Bisque)
colnames(RESULTS_Bisque) <-c("CT","tissue","observed_values")

# deconvSeq
singlecelldata = C.eset 
celltypes.sc = as.character(phenoDataC$cellType) #To avoid "Design matrix not of full rank" when removing 1 CT 
tissuedata = T.eset 
ptm <- proc.time()
design.singlecell = model.matrix(~ -1 + as.factor(phenoDataC$cellType))
colnames(design.singlecell) = levels(as.factor(celltypes.sc))
rownames(design.singlecell) = colnames(singlecelldata)

dge.singlecell = deconvSeq::getdge(singlecelldata,design.singlecell, ncpm.min = 1, nsamp.min = 3, method = "bin.loess")
saveRDS(dge.singlecell, file = "dge_singlecell.rds")
b0.singlecell = deconvSeq::getb0.rnaseq(dge.singlecell, design.singlecell, ncpm.min =1, nsamp.min = 3)
saveRDS(b0.singlecell, file = "b0_singlecell.rds")
# dge.singlecell <- readRDS("dge_singlecell.rds")
# b0.singlecell <- readRDS("b0_singlecell.rds")
dge.tissue = deconvSeq::getdge(tissuedata, NULL, ncpm.min = 1, nsamp.min = 3, method = "bin.loess")

RESULTS_deconv = t(deconvSeq::getx1.rnaseq(NB0 = "top_fdr",b0.singlecell, dge.tissue)$x1) #genes with adjusted p-values <0.05 after FDR correction
time_deconv = proc.time() - ptm
RESULTS_deconv = RESULTS_deconv[gtools::mixedsort(rownames(RESULTS_deconv)),]
RESULTS_deconv = data.table::melt(RESULTS_deconv)
colnames(RESULTS_deconv) <-c("CT","tissue","observed_values")

# BayesPrism
ptm <- proc.time()
myPrism <- BayesPrism::new.prism(reference=t(C),mixture=,t(T),input.type="count.matrix",cell.type.labels=phenoDataC$cellType,cell.state.labels=phenoDataC$cellType,key=NULL,outlier.cut=0.01,outlier.fraction=0.1)
bp.res <- BayesPrism::run.prism(prism=myPrism,n.cores=2)
RESULTS_Prism = t(BayesPrism::get.fraction(bp=bp.res,which.theta="final",state.or.type="type"))
time_Prism <- proc.time() - ptm
RESULTS_Prism = RESULTS_Prism[gtools::mixedsort(rownames(RESULTS_Prism)),]
RESULTS_Prism = data.table::melt(RESULTS_Prism)
colnames(RESULTS_Prism) <-c("CT","tissue","observed_values")

process_results <- function(P, RESULTS, method_name) {  
	RESULTS = merge(RESULTS, P)
	RESULTS$expected_values <- round(RESULTS$expected_values, 3)
	RESULTS$observed_values <- round(RESULTS$observed_values, 3)
	RESULTS$method <- rep(method_name, nrow(RESULTS))
  
	return(RESULTS)
}

if (!is.null(P)) {
	P = P[gtools::mixedsort(rownames(P)), ]
	P$CT = rownames(P)
	P = data.table::melt(P, id.vars = "CT")
	colnames(P) <- c("CT", "tissue", "expected_values")
	
	RESULTS_Bisque <- process_results(P, RESULTS_Bisque, "Bisque")
	RESULTS_DWLS <- process_results(P, RESULTS_DWLS, "DWLS")
	RESULTS_combined1 <- process_results(P, RESULTS_combined1, "DWLS_Bisque")
	RESULTS_combined2 <- process_results(P, RESULTS_combined2, "Music_Bisque")
	RESULTS_Music <- process_results(P, RESULTS_Music, "MuSic")
	RESULTS_scdc <- process_results(P, RESULTS_scdc, "SCDC")
	RESULTS_deconv <- process_results(P, RESULTS_deconv, "deconv")
	RESULTS_Prism <- process_results(P, RESULTS_Prism, "BayesPrism")
}

summarize_results <- function(RESULTS, time_data) {
	RESULTS_summary <- RESULTS %>%
		dplyr::summarise(
			RMSE = sqrt(mean((observed_values - expected_values)^2)) %>% round(., 4),
			Pearson = cor(observed_values, expected_values) %>% round(., 4)
		)
	RESULTS_summary$Running_time = round(as.numeric(time_data["elapsed"]), 3)
	return(RESULTS_summary)
}

RESULTS_Bisque <- summarize_results(RESULTS_Bisque, time_Bisque)
RESULTS_DWLS <- summarize_results(RESULTS_DWLS, time_DWLS)
RESULTS_combined1 <- summarize_results(RESULTS_combined1, time_combined1)
RESULTS_combined2 <- summarize_results(RESULTS_combined2, time_combined2)
RESULTS_Music <- summarize_results(RESULTS_Music, time_Music)
RESULTS_scdc <- summarize_results(RESULTS_scdc, time_scdc)
RESULTS_deconv <- summarize_results(RESULTS_deconv, time_deconv)
RESULTS_Prism <- summarize_results(RESULTS_Prism, time_Prism)

combine.results <- rbind(RESULTS_combined1, RESULTS_combined2, RESULTS_Bisque, RESULTS_DWLS, RESULTS_Music, RESULTS_scdc, RESULTS_deconv, RESULTS_Prism)
rownames(combine.results) <- c("DWLS_Bisque", "Music_Bisque", "Bisque", "DWLS", "MuSic", "SCDC", "deconv", "BayesPrism")
write.csv(combine.results, "results.csv")