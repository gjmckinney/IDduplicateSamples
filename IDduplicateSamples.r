library(vcfR)
library(pbapply)
library(ggplot2)
library(stringr)


#function to compare genotypes between all pairs of samples
#input format: rows are loci, columns are samples, alleles are 0 and 1, missing genotypes are NA
#heterozygous calls must be consistent within a locus, ie all 0/1, no 1/0
IDduplicateSamples<-function(genotypes,MAF=NULL){
  #function to calculate MAF
  calcMAF<-function(locusGenos){
    allele1Counts<-sum(str_count(locusGenos,"0"),na.rm=TRUE)
    allele2Counts<-sum(str_count(locusGenos,"1"),na.rm=TRUE)
    allele1Freq<-allele1Counts/sum(allele1Counts,allele2Counts)
    if(allele1Freq>0.5){
      MAF<-1-allele1Freq
    }else{
      MAF<-allele1Freq
    }
    return(MAF)
  }
  
  #filter loci using MAF if threshold is specified
  if(!is.null(MAF)){
    #calculate MAF
    message(paste("MAF threshold applied:",MAF,"MAF",sep=" "))
    message("calculating MAF")
    locusMAF<-pbapply(genotypes,1,calcMAF)
    #convert to dataframe
    locusMAF<-data.frame(locus_ID=names(locusMAF),MAF=locusMAF,row.names=NULL)
    locusMAF$locus_ID<-as.character(locusMAF$locus_ID)
    #filter loc based on MAF threshold
    genotypes<-genotypes[rownames(genotypes)%in%locusMAF$locus_ID[locusMAF$MAF>=MAF],]
  }else{
    message("No MAF threshold applied, using all loci")
  }
  
  #make matrix of called vs NA genotypes for faster counting of missing data
  genotypes_NAmatrix<-genotypes
  genotypes_NAmatrix[!is.na(genotypes_NAmatrix)]<-0
  genotypes_NAmatrix[is.na(genotypes_NAmatrix)]<-1
  class(genotypes_NAmatrix)<-"numeric"
  
  #identify all unique pairs of samples
  allPairs<-combn(dim(genotypes)[2], 2)
  ncombo<-dim(allPairs)[2]
  nloci<-dim(genotypes)[1]
  message(paste("number of samples:",dim(genotypes)[2],sep=" "))
  message(paste("number of loci:",nloci,sep=" "))
  message(paste("number of sample pairs:",ncombo,sep=" "))
  #reshape into 2xn matrix
  allPairs<-matrix(allPairs,nrow=2)
  
  #function to compare genotypes
  compareGenos<-function(samplePair){
    #count number of loci that are genotyped in both samples
    NAcounts<-genotypes_NAmatrix[,samplePair[1]]+genotypes_NAmatrix[,samplePair[2]]
    sharedCounts<-nloci-(sum(NAcounts)-sum(NAcounts[NAcounts==2])/2)
    #count number of loci with genotypes that match in both samples
    genotypeMatches<-sum(genotypes[,samplePair[1]]==genotypes[,samplePair[2]],na.rm=TRUE)
    #return counts of genotype matches and shared loci
    return(c(genotypeMatches,sharedCounts))
  }
  
  #do all pairwise sample comparisons
  message("comparing genotypes")
  matches<-pbapply(allPairs,2,compareGenos)
  
  #make dataframe of results
  comparisonResults<-data.frame(matrix(NA,nrow=dim(allPairs)[2],ncol=7))
  colnames(comparisonResults)<-c("Sample1","Sample2","matchedGenotypes","commonGenotypes","proportionMatch","proportionCommon","totalLoci")
  comparisonResults$Sample1<-colnames(genotypes)[allPairs[1,]]
  comparisonResults$Sample2<-colnames(genotypes)[allPairs[2,]]
  comparisonResults$matchedGenotypes<-matches[1,]
  comparisonResults$commonGenotypes<-matches[2,]
  comparisonResults$proportionMatch<-comparisonResults$matchedGenotypes/comparisonResults$commonGenotypes
  comparisonResults$proportionCommon<-comparisonResults$commonGenotypes/nloci
  comparisonResults$totalLoci<-nloci
  return(comparisonResults)
}


