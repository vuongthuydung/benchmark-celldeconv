# modified from original repo: https://github.com/favilaco/deconv_benchmark

Generator <- function(sce, phenoData, Num.mixtures=1000, pool.size=100, min.percentage = 1, max.percentage = 99, seed=24){
    CT = unique(phenoData$cellType) # name of cell types
    ?stopifnot(length(CT) >= 2)
    # set.seed(10)
    b <- abs(rnorm(nrow(sce), mean = 0, sd = 0.5)) + 1
    a <- abs(rnorm(nrow(sce), mean = 0, sd = 0.5)) 
    set.seed(seed)
    require(dplyr)
    require(gtools)  
    
    # create a frame of 2 columns: CT and its frequency
    cell.distribution = data.frame(table(phenoData$cellType),stringsAsFactors = FALSE) 
    colnames(cell.distribution) = c("CT","max.n")
    
    Tissues = list() 
    Proportions = list()

    for(y in 1:Num.mixtures){
        while(!exists("P")){
    
           num.CT.mixture = sample(x = 2:length(CT),1) # choose the number of cell types to be present in the pseudo bulk 
           selected.CT = sample(CT, num.CT.mixture, replace=FALSE) # vectors of cell types present

           P = runif(num.CT.mixture, min.percentage, max.percentage) # choose proportion for each cell type present in the above vector
           P = round(P/sum(P), digits = log10(pool.size))  #sum to 1
           P = data.frame(CT = selected.CT, expected = P, stringsAsFactors = FALSE)
           missing.CT = CT[!CT %in% selected.CT]
           missing.CT = data.frame(CT = missing.CT, expected = rep(0, length(missing.CT)), stringsAsFactors = FALSE)
          
           P = rbind.data.frame(P, missing.CT)
           potential.mix = merge(P, cell.distribution)
           # Ensure the total number of original cells from one cell type are no less than the number of cells of the same type chosen for a pseudo bulk
           potential.mix$size = potential.mix$expected * pool.size
           if( !all(potential.mix$max.n >= potential.mix$size) | sum(P$expected) != 1){
                rm(list="P")
           }
        }
        # Using info in P to build T simultaneously
        chosen_cells <- sapply(which(P$expected != 0), function(x){
            n.cells = P$expected[x] * pool.size
            chosen = sample(phenoData$cellID[phenoData$cellType == P$CT[x]],
                      n.cells)
            chosen
        }) %>% unlist()

        T <- Matrix::rowSums(sce[,colnames(sce) %in% chosen_cells]) %>% as.data.frame()
        colnames(T) = paste("mix",y,sep="")
        P = P[,c("CT","expected")]
        P$mix = paste("mix",y,sep="")

        for(i in 1:nrow(T)){
            T[i,] <- b[i] * T[i,] + a[i]
        }

        Tissues[[y]] <- T
        Proportions[[y]] <- P
        rm(list=c("T","P","chosen_cells","missing.CT"))

    }
    P = do.call(rbind.data.frame, Proportions)

    T = do.call(cbind.data.frame, Tissues)
    P = data.table::dcast(P, CT ~ mix, 
                        value.var = "expected",
                        fun.aggregate = sum) %>% data.frame(.,row.names = 1)
    P = P[,gtools::mixedsort(colnames(P))]
    return(list(T = T, P = P))
}