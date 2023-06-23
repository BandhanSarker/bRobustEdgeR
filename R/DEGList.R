#' @title Sample RNA-seq count data with outliers value
#' @description sample RNAseq count \code{data-matrix} with outlier’s value that is used for testing
#' @format RNAseq count \code{data_matrix},whose row contains 1000 genes and column contains 6 samples, where 3 control groups and 3 case groups
#' @examples
#' data(sampleDataOut)
"sampleDataOut"

#' @title Sample RNA-seq count data with missing value
#' @description sample RNAseq count \code{data-matrix} with missing value that is used for testing
#' @format RNAseq count \code{data_matrix},whose row contains 1000 genes and column contains 6 samples, where 3 control groups and 3 case groups
#' @examples
#' data(sampleDataMiss)
"sampleDataMiss"
##################################################################################
######                           Boosting Robust edgeR                 ######
##################################################################################
#' Identify differentially expressed gene list
#' @description This function is used to identify differentially expressed genes from RNA-seq count data.
#' @param data  RNA-seq count data matrix, whose row contain genes and column contain samples
#' @param n1    #Control group
#' @param n2    #Case group
#' @param p.threshold adjusted pvalue by default=0.05
#' @return Differentially expressed genes list
#' @import edgeR
#' @import limma
#' @import missForest
#' @import stats
#' @examples
#' Brobustedger::DEGList(sampleDataOut,3,3,0.05)
#' Brobustedger::DEGList(sampleDataMiss,3,3,0.05)
#' @export
DEGList<-function(data,n1,n2,p.threshold=0.05){

  ##############################
  message("Brobustedger::INFO:Computation takes few times... ...")
  #check_missing_values or not
  set.seed(100)
  if (sum(is.na(data)) > 0) {
    mm<-missForest::missForest(data, maxiter = 5, ntree = 100)
    nmdata<-mm$ximp
      } else {
    nmdata<-data
  }
  #############
  sub<-nmdata
  #estimate sequencing depth and compute cutoff
  sd <- mean(colSums(sub,na.rm=T))
  sdcut <- 1/sd
  suppressWarnings(
    #implement iterative scheme
    outl <- apply(sub, 1, function(z) {
      x <- z
      nbp <- 0
      out <- rep(0,0)
      track <- c(1:length(x))

      #iterative scheme
      while((min(nbp,na.rm=T) < sdcut) & (length(x)>2)) {

        #build matrix with rows representing leave-one out observation
        mat <- matrix(rep(x,length(x)),ncol=length(x),byrow=T)
        diag(mat) <- NA
        tmp <- t(mat)
        mat <- t(matrix(tmp[!is.na(tmp)],nrow=(length(x)-1),ncol=(length(x))))

        #fit negative binomial or Poisson distribution
        nbfit <- apply(mat,1,function(y) {
          if(length(y)>1) {
            v <- var(y)
            m <- mean(y)
            if(all(y==0)=="TRUE") {
              output <- NA
            } else if (v>m) {
              p <- mean(y)/var(y)
              r <- mean(y)^2/(var(y)-mean(y))
              output <- c(p,r)
            } else {
              lamb <- mean(y)
              output <- c(lamb)
            }
          } else output <- NA

          list(output)
        })

        nbfit <- lapply(nbfit, "[[", 1)
        #compute probabilities for leave-one out observation
        nbp <- rep(0,0)
        for (i in 1:length(nbfit)) {
          if(length(nbfit[[i]])==2) {
            nbp <- c(nbp,dnbinom(x[i],prob=nbfit[[i]][1],size=nbfit[[i]][2]))
          } else {
            nbp <- c(nbp,dpois(x[i],lambda=nbfit[[i]][1]))
          }
        }

        #compare probabilities to cutoff
        sel <- which(nbp < sdcut)
        if(length(sel)>0) x <- x[-sel]
        if(length(out)==0) {
          out <- c(out,track[sel])
        } else {
          out <- c(out,track[-out][sel])
        }

        fout <- rep(NA,length(z))
        if(length(out)>0) {
          fout[out] <- z[out]
        }
      }

      list(fout)
    })
  )
  #new data matrix with outliers (all other data is NA’ed)
  identout <- matrix(unlist(outl),nrow=nrow(sub),ncol=ncol(sub),byrow=T)
  colnames(identout) <- colnames(sub)
  rownames(identout) <- rownames(sub)

  outData<-identout
  imputeData<-nmdata
  ######Test each gene whether it is an outlying gene or not using iLOO######
  location<-which(!is.na(outData))
  if(length(location)==0){
    dataM<-nmdata
  }else{
    ###Outlier’s value in each row and column should be considered as NULL#########

    imputeData[!is.na(outData)]<-NA
    dataM<-imputeData
  }


  ########Imputation of nullify values using RF############
  m<-missForest::missForest(dataM, maxiter = 5, ntree = 100)
  dataC<-m$ximp

  ########
  dataR<-round(dataC,0)
  # Assign condition
  control<-"control";case<-"case"
  condition <- factor(c(rep(control, n1), rep(case, n2)))
  coldata <- data.frame(row.names=colnames(dataR), condition)
  ################

  dge <- edgeR::DGEList(counts=dataR, group=coldata$condition)
  # Create the contrast matrix
  design.mat <-model.matrix(~ 0 + dge$samples$group)
  colnames(design.mat) <- levels(dge$samples$group)

  # Estimate dispersion parameter for GLM
  dge<-edgeR::estimateGLMRobustDisp(dge,  design.mat,maxit = 5, residual.type = "pearson")

  # Design matrix
  design.mat <-model.matrix(~ 0 + dge$samples$group)
  colnames(design.mat) <- c(control, case)
  # Model fitting

  fit.edgeR <- edgeR::glmFit(dge, design.mat)
  # Differential expression
  contrasts.edgeR <- limma::makeContrasts(control - case, levels=design.mat)
  lrt.edgeR <- edgeR::glmLRT(fit.edgeR, contrast=contrasts.edgeR)
  #logfold
  edgeR_results <- lrt.edgeR$table
  sig.edgeR <- edgeR::decideTestsDGE(lrt.edgeR, adjust.method="BH", p.value = p.threshold)
  genes.edgeR <- row.names(edgeR_results)[which(sig.edgeR != 0)]
  return(genes.edgeR)
}
#####################################################################################
####                                END                                   ###########
#####################################################################################

