# modified from original repos: 
# 1. https://github.com/cozygene/bisque
# 2. https://github.com/dtsoucas/DWLS

CountsToCPM <- function(eset) {
  Biobase::exprs(eset) <- base::sweep(Biobase::exprs(eset),
                                      2, base::colSums(Biobase::exprs(eset)),
                                      `/`) * 1000000
  indices <- base::apply(Biobase::exprs(eset), MARGIN=2,
                         FUN=function(column) {base::anyNA(column)})
  if (base::any(indices)) {
    n.cells <- base::sum(indices)
    base::stop(base::sprintf("Zero expression in selected genes for %i cells",
                             n.cells))
  }
  return(eset)
}

FilterZeroVarianceGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, stats::var) != 0)
  indices <- indices & (! base::is.na(indices))
  if (base::sum(indices) < base::length(indices)) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::message(base::sprintf("Filtered %i zero variance genes.",
                                genes.filtered))
  }
  return(eset)
}

FilterUnexpressedGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, base::sum) != 0)
  indices <- indices & (! base::is.na(indices))
  if (base::sum(indices) < base::length(indices)) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::message(base::sprintf("Filtered %i unexpressed genes.",
                                genes.filtered))
  }
  return(eset)
}

GetOverlappingSamples <- function(sc.eset, bulk.eset, subject.names, verbose) {
  bulk.samples <- Biobase::sampleNames(bulk.eset)
  sc.samples <- base::levels(base::factor(sc.eset[[subject.names]]))
  overlapping.samples <- base::intersect(bulk.samples, sc.samples)
  if (base::length(overlapping.samples) == 0) {
    base::stop("No overlapping samples in bulk and single-cell expression.")
  }
  remaining.samples <- base::setdiff(Biobase::sampleNames(bulk.eset),
                                     overlapping.samples)
  if (base::length(remaining.samples) == 0) {
    base::stop("All samples have single-cell data, nothing to process.")
  }
  samples <- base::list("overlapping"=overlapping.samples,
                        "remaining"=remaining.samples)
  if (verbose) {
    n.overlapping <- base::length(samples$overlapping)
    n.remaining <- base::length(samples$remaining)
    base::message(base::sprintf("Found %i samples ", n.overlapping),
                                "with bulk and single-cell expression.")
    base::message(base::sprintf("Remaining %i ", n.remaining),
                                "bulk samples will be decomposed.")
  }
  return(samples)
}

GetOverlappingGenes <- function(sc.eset, bulk.eset, markers, verbose) {
  bulk.genes <- Biobase::featureNames(bulk.eset)
  sc.genes <- Biobase::featureNames(sc.eset)
  overlapping.genes <- base::intersect(bulk.genes, sc.genes)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste0("No overlapping genes found between bulk and ",
                           "single-cell expression."))
  }
  overlapping.genes <- base::intersect(overlapping.genes, markers)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste0("No marker genes found in both bulk and ",
                           "single-cell expression."))
  }
  if (verbose) {
    n.genes <- base::length(overlapping.genes)
    base::message(base::sprintf("Using %i genes in both", n.genes),
                                " bulk and single-cell expression.")
  }
  return(overlapping.genes)
}

GenerateSCReference <- function(sc.eset, cell.types) {
  cell.labels <- base::factor(sc.eset[[cell.types]])
  all.cell.types <- base::levels(cell.labels)
  aggr.fn <- function(cell.type) {
    base::rowMeans(Biobase::exprs(sc.eset)[,cell.labels == cell.type, drop=F])
  }
  template <- base::numeric(base::nrow(sc.eset))
  sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
  return(sc.ref)
}

SupervisedTransformBulk <- function(gene, Y.train, X.train, X.pred) {
  Y.train.scaled <- base::scale(Y.train[gene,,drop=T])
  Y.center <- base::attr(Y.train.scaled, "scaled:center")
  Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
  X.train.scaled <- base::scale(X.train[gene,,drop=T])
  X.center <- base::attr(X.train.scaled, "scaled:center")
  X.scale <- base::attr(X.train.scaled, "scaled:scale")

  if (base::anyNA(X.train.scaled) & base::anyNA(Y.train.scaled)) {
    coeff <- Y.train[gene,,drop=T][1]/X.train[gene,,drop=T][1]
    if (coeff == 0 || ! is.finite(coeff)) {
      coeff = NaN
    }
    Y.pred <- base::matrix(X.pred[gene,,drop=T] * coeff,
                           dimnames=base::list(base::colnames(X.pred), gene))
  }
  
  else if (anyNA(X.train.scaled) || anyNA(Y.train.scaled)) {
    Y.pred <- base::matrix(X.pred[gene,,drop=T] * NaN,
                           dimnames=base::list(base::colnames(X.pred), gene))
  }

  else {
    X.pred.scaled <- base::scale(X.pred[gene,,drop=T],
                                 center=X.center,
                                 scale=X.scale)
    model <- stats::lm(Y.train.scaled ~ X.train.scaled +0)
    coeff <- base::as.numeric(stats::coefficients(model))
    Y.pred.scaled <- X.pred.scaled * coeff
    Y.pred <- base::matrix((Y.pred.scaled * Y.scale) + Y.center,
                           dimnames=base::list(base::colnames(X.pred), gene))
  }
  return(Y.pred)
}

SemisupervisedTransformBulk <- function(gene, Y.train, X.pred) {

  Y.train.scaled <- base::scale(Y.train[gene,,drop=T])
  Y.center <- base::attr(Y.train.scaled, "scaled:center")
  Y.scale <- base::attr(Y.train.scaled, "scaled:scale")
  n <- base::length(Y.train.scaled)
  
  shrink.scale <- base::sqrt(base::sum((Y.train[gene,,drop=T]-Y.center)^2)/(n+1))
  X.pred.scaled <- base::scale(X.pred[gene,,drop=T])
  Y.pred <- base::matrix((X.pred.scaled * shrink.scale) + Y.center,
                         dimnames=base::list(base::colnames(X.pred), gene))
  return(Y.pred)
}

ReferenceBasedDecomposition <- function(bulk.eset,
                                        sc.eset,
                                        markers=NULL,
                                        cell.types="cellType",
                                        subject.names="SubjectName",
                                        use.overlap=TRUE, 
                                        verbose=TRUE,
                                        old.cpm=TRUE,
                                        decom_type) {
  if ((! methods::is(sc.eset, "ExpressionSet")) || 
      (! methods::is(bulk.eset, "ExpressionSet"))) {
    base::stop("Expression data should be in ExpressionSet")
  }
  else if (! cell.types %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::sprintf("Cell type label \"%s\" ", cell.types),
               "not found in single-cell ExpressionSet varLabels.")
  }
  else if (! subject.names %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::sprintf("Individual label \"%s\"", subject.names),
               " not found in single-cell ExpressionSet varLabels.")
  }
  n.sc.individuals <-
    base::length(base::levels(base::factor(sc.eset[[subject.names]])))
  if (n.sc.individuals == 1) {
    base::stop("Only one individual detected in single-cell data. At least ",
               "two subjects are needed (three or more recommended).")
  }
  else if (n.sc.individuals == 2) {
    base::warning("Only two individuals detected in single-cell data. While ",
                  "Bisque will run, we recommend at least three subjects for",
                  " reliable performance.")
  }
  n.cell.types <-
    base::length(base::levels(base::factor(sc.eset[[cell.types]])))
  if (n.cell.types == 1) {
    base::stop("Single-cell pheno data indicates only one cell type",
               " present. No need for decomposition.")
  }
  if (verbose) {
    base::message(base::sprintf("Decomposing into %i cell types.",
                                n.cell.types))
  }
  if (use.overlap) {
    samples <- GetOverlappingSamples(sc.eset, bulk.eset, subject.names, verbose)
  }
  if (base::is.null(markers)) {
    markers <- Biobase::featureNames(sc.eset)
  }
  else {
    markers <- base::unique(base::unlist(markers))
  }
  genes <- GetOverlappingGenes(sc.eset, bulk.eset, markers, verbose)
  if (old.cpm) {
    sc.eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(sc.eset)[genes,],
                             phenoData=sc.eset@phenoData)
    bulk.eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(bulk.eset)[genes,],
                             phenoData=bulk.eset@phenoData)
  }
  if (verbose) {
    base::message("Converting single-cell counts to CPM and ",
              "filtering zero variance genes.")
  }
  sc.eset <- CountsToCPM(sc.eset)
  if (!old.cpm) {
    sc.eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(sc.eset)[genes,],
                             phenoData=sc.eset@phenoData)
  }
  sc.eset <- FilterZeroVarianceGenes(sc.eset, verbose)
  if (verbose) {
    base::message("Converting bulk counts to CPM and filtering",
                  " unexpressed genes.")
  }
  bulk.eset <- CountsToCPM(bulk.eset)
  if (!old.cpm) {
    bulk.eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(bulk.eset)[genes,],
                             phenoData=bulk.eset@phenoData)
  }
  bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose)
  genes <- base::intersect(Biobase::featureNames(sc.eset),
                           Biobase::featureNames(bulk.eset))
  if (base::length(genes) == 0) {
    base::stop("Zero genes remaining after filtering and ",
               "intersecting bulk, single-cell, and marker genes.")
  }
  if (verbose) {
    n.cells <- base::ncol(sc.eset)
    base::message("Generating single-cell based reference from ",
                  sprintf("%i cells.\n", n.cells))
  }
  sc.ref <- GenerateSCReference(sc.eset, cell.types)[genes,,drop=F] # sig
  sc.props <- BisqueRNA::CalculateSCCellProportions(sc.eset, subject.names, cell.types)
  sc.props <- sc.props[base::colnames(sc.ref),,drop=F]
  if (use.overlap) {
    if (verbose) {
      base::message("Learning bulk transformation from overlapping samples.")
    }
    # Y.train is pseudo-bulk expression based on reference profile weighted by
    #   cell type proportions estimated for single-cell samples.
    Y.train <- sc.ref %*% sc.props[,samples$overlapping,drop=F]
    # X.train is the actual bulk for the single-cell samples.
    X.train <- Biobase::exprs(bulk.eset)[genes,samples$overlapping,drop=F]
    # X.pred is the bulk for the remaining samples to be decomposed.
    X.pred <- Biobase::exprs(bulk.eset)[genes,samples$remaining,drop=F]
    template <- base::numeric(base::length(samples$remaining))
    base::names(template) <- samples$remaining
    if (verbose) {
      base::message("Applying transformation to bulk samples and decomposing.")
    }
    # Y.pred is the transformed bulk for samples to be decomposed.
    Y.pred <- base::matrix(base::vapply(X=genes, FUN=SupervisedTransformBulk,
                                        FUN.VALUE=template,
                                        Y.train, X.train, X.pred,
                                        USE.NAMES=TRUE),
                           nrow=base::length(samples$remaining))
    sample.names <- samples$remaining
  }
  else {
    if (verbose) {
      base::message("Inferring bulk transformation from single-cell alone.")
    }
    # Y.train is pseudo-bulk expression based on reference profile weighted by
    #   cell type proportions estimated for single-cell samples.
    Y.train <- sc.ref %*% sc.props
    # X.pred is the bulk for the remaining samples to be decomposed.
    X.pred <- Biobase::exprs(bulk.eset)[genes,,drop=F]
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    template <- base::numeric(base::length(sample.names))
    base::names(template) <- sample.names
    if (verbose) {
      base::message("Applying transformation to bulk samples and decomposing.")
    }
    # Y.pred is the transformed bulk for samples to be decomposed.
    Y.pred <- base::matrix(base::vapply(X=genes,
                                        FUN=SemisupervisedTransformBulk,
                                        FUN.VALUE=template,
                                        Y.train, X.pred,
                                        USE.NAMES=TRUE),
                           nrow=base::length(sample.names))
  }
  # Columns in Y.pred with NaN indicate transformation could not be learned
  #   for that gene. 
  indices <- base::apply(Y.pred, MARGIN=2,
                         FUN=function(column) {base::anyNA(column) || base::any(column<=0)})
  if (base::any(indices)) {
    if (verbose) {
      n.dropped <- base::sum(indices)
      base::message(base::sprintf("Dropped an additional %i genes", n.dropped),
                    " for which a transformation could not be learned.")
    }
    if (sum(!indices) == 0) {
      base::stop("Zero genes left for decomposition.")
    }
    Y.pred <- Y.pred[,!indices,drop=F]
    sc.ref <- sc.ref[!indices,,drop=F]
  }
  if (decom_type == "new1") {
    results <- NULL
    # signature <- signature[!indices,,drop=F]
    for(id in seq_len(nrow(Y.pred))){
      subresult <- solveDampenedWLS(sc.ref,Y.pred[id,])
      results<-cbind(results,subresult)
    }
    results <- as.matrix(results)

    Y.pred <- base::t(Y.pred) 
    base::rownames(Y.pred) <- base::rownames(sc.ref)
    base::colnames(Y.pred) <- sample.names
    base::rownames(results) <- base::colnames(sc.ref)
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    base::colnames(results) <- sample.names
    results <- base::list(bulk.props=results[base::colnames(sc.ref),,drop=F],
                        sc.props=sc.props,
                        genes.used=base::rownames(sc.ref),
                        transformed.bulk=Y.pred)
  }

  if (decom_type == "new2") {
    sc.eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(sc.eset)[genes,],
                             phenoData=sc.eset@phenoData)

    sce <- as(as(sc.eset[!indices,,drop=F], "SummarizedExperiment"), "SingleCellExperiment")
    names(assays(sce))=c("counts")
    Y.pred <- base::t(Y.pred) 
    base::rownames(Y.pred) <- base::rownames(sc.ref)
    base::colnames(Y.pred) <- sample.names
    T.eset <- Biobase::ExpressionSet(assayData=Y.pred)
    results = t(MuSiC::music_prop(bulk.mtx = exprs(T.eset), sc.sce = sce, clusters = 'cellType',
                              markers = NULL, normalize = FALSE, samples = 'SubjectName', 
                              verbose = FALSE)$Est.prop.weighted)
    base::rownames(results) <- base::colnames(sc.ref)
    sample.names <- base::colnames(Biobase::exprs(bulk.eset))
    base::colnames(results) <- sample.names
    results <- base::list(bulk.props=results)
  }
  
  return(results)
}

ReferenceBasedDecompositionNew1 <- function(bulk.eset,
                                            sc.eset,
                                            markers=NULL,
                                            cell.types="cellType",
                                            subject.names="SubjectName",
                                            use.overlap=TRUE, 
                                            verbose=TRUE,
                                            old.cpm=TRUE){
  return(ReferenceBasedDecomposition(bulk.eset,sc.eset,markers,cell.types,subject.names,use.overlap,verbose,old.cpm,decom_type="new1"))
}

ReferenceBasedDecompositionNew2 <- function(bulk.eset,
                                            sc.eset,
                                            markers=NULL,
                                            cell.types="cellType",
                                            subject.names="SubjectName",
                                            use.overlap=TRUE, 
                                            verbose=TRUE,
                                            old.cpm=TRUE){
  return(ReferenceBasedDecomposition(bulk.eset,sc.eset,markers,cell.types,subject.names,use.overlap,verbose,old.cpm,decom_type="new2"))
}

#----------------
## Modify DWLS
#return cell number, not proportion
solveOLSInternal<-function(S,B,nnls=FALSE){
  if(!nnls){
    D<-t(S)%*%S
    d<-t(S)%*%B
    A<-cbind(diag(dim(S)[2]))
    bzero<-c(rep(0,dim(S)[2]))
    solution<-solve.QP(D/norm(D,"2"),d,A,bzero)$solution
  }
  else{
    solution<-nnls(S,B)
    solution<-solution$x
  }
  names(solution)<-colnames(S)
  return(solution)
}

#solve using WLS with weights dampened by a certain dampening constant
solveDampenedWLS<-function(S,B){
  solution<-solveOLSInternal(S,B)
  iterations<-0
  changes<-c()
  #find dampening constant for weights using cross-validation
  j<-findDampeningConstant(S,B,solution)
  change<-1
  while(change>.01 & iterations<1000){
    newsolution<-solveDampenedWLSj(S,B,solution,j)
    #decrease step size for convergence
    solutionAverage<-rowMeans(cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
    change<-norm(as.matrix(solutionAverage-solution))
    solution<-solutionAverage
    iterations<-iterations+1
    changes<-c(changes,change)
  }
  return(solution/sum(solution))
}

#find a dampening constant for the weights using cross-validation
findDampeningConstant<-function(S,B,goldStandard){
  solutionsSd<-NULL
  #goldStandard is used to define the weights
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsScaledMinusInf<-wsScaled
  #ignore infinite weights
  if(max(wsScaled)=="Inf"){
    wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
  }
  #try multiple values of the dampening constant (multiplier)
  #for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:ceiling(log2(max(wsScaledMinusInf)))){
    multiplier<-1*2^(j-1)
    wsDampened<-wsScaled
    wsDampened[which(wsScaled>multiplier)]<-multiplier
    solutions<-NULL
    seeds<-c(1:100)
    for (i in 1:100){
      set.seed(seeds[i]) #make nondeterministic
      subset<-sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
      #solve dampened weighted least squares for subset
      fit = lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
      sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
      solutions<-cbind(solutions,sol)
    }
    solutionsSd<-cbind(solutionsSd,apply(solutions,1,sd))
  }
  #choose dampening constant that results in least cross-validation variance
  j<-which.min(colMeans(solutionsSd^2))
  return(j)
}

#solve using WLS with weights dampened by a certain dampening constant
solveDampenedWLS_nonopt<-function(S,B){
  solution<-solveOLSInternal(S,B)
  iterations<-0
  changes<-c()
  #find dampening constant for weights using cross-validation
  j<-findDampeningConstant_nonopt(S,B,solution)
  change<-1
  while(change>.01 & iterations<1000){
    newsolution<-solveDampenedWLSj(S,B,solution,j)
    solutionAverage<-rowMeans(cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
    change<-norm(as.matrix(solutionAverage-solution))
    solution<-solutionAverage
    iterations<-iterations+1
    changes<-c(changes,change)
  }
  return(solution/sum(solution))
}

findDampeningConstant_nonopt<-function(S,B,goldStandard){
  solutionsSd<-NULL
  #goldStandard is used to define the weights
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsScaledMinusInf<-wsScaled
  #ignore infinite weights
  if(max(wsScaled)=="Inf"){
    wsScaledMinusInf<-wsScaled[-which(wsScaled=="Inf")]
  }

  j <- quantile(wsScaledMinusInf, probs = 0.9)
  return(j)
}

solveDampenedWLSj<-function(S,B,goldStandard,j,nnls=FALSE){
  multiplier<-1*2^(j-1)
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsDampened<-wsScaled
  wsDampened[which(wsScaled>multiplier)]<-multiplier
  if(!nnls){
    W<-diag(wsDampened)
    D<-t(S)%*%W%*%S
    d<- t(S)%*%W%*%B
    A<-cbind(diag(dim(S)[2]))
    bzero<-c(rep(0,dim(S)[2]))
    sc <- norm(D,"2")
    solution<-solve.QP(D/sc,d/sc,A,bzero)$solution
  }
  else{
    B_weight<-B*sqrt(wsDampened)
    S_weight<-sweep(S, 1, sqrt(wsDampened), '*')
    solution<-nnls(S_weight,B_weight)
    solution<-solution$x
  }

  names(solution)<-colnames(S)
  return(solution)
}