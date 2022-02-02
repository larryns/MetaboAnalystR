# Scripts to process metabolomics data using MetaboanalystR. Most of the code in here is replicated from 
# running an analysis on the webpage version.

library(MetaboAnalystR)
library(magrittr)
library(ggplot2)
library(stats)

RootDir <- "/data/singhln/Projects/Kilbaugh/Metabolomics/"
DataDir <- paste(RootDir, "Data/", sep="")
MBDir <- paste(RootDir, "MetaboAnalyst/", sep="")
ResultsDir <- paste0(MBDir, "Results/")

# Mappings to collapse the different amino acids:
AAMap <- list()
AAMap[["Essential"]] <- c("AA_His", "AA_IL", "AA_Leu", "AA_Lys", "AA_METH", "AA_Phe", "AA_Thr", "AA_Trp", "AA_Val")
AAMap[["NonEssential"]] <- c("AA_Ala", "AA_Arg", "AA_Asn", "AA_Asp", "AA_Cys", "AA_Gln", "AA_Glu", "AA_Gly", "AA_Pro", "AA_Ser", "AA_Tyr")
AAMap[["Glucogenic"]] <- c("AA_Ala", "AA_Arg", "AA_Asn", "AA_Asp", "AA_Cys", "AA_Gln", "AA_Glu", "AA_Gly","AA_His", "AA_IL", "AA_METH", "AA_Phe", "AA_Pro", "AA_Ser", "AA_Thr", "AA_Trp", "AA_Tyr", "AA_Val"
)
AAMap[["Ketogenic"]] <- c("AA_IL", "AA_Leu", "AA_Lys", "AA_Phe", "AA_Thr", "AA_Trp", "AA_Tyr")
AAMap[["BCAA"]] <- c("AA_IL", "AA_Leu", "AA_Val")
AAMap[["Pyr"]] <- c("AA_Ala", "AA_Cys", "AA_Gly", "AA_Ser", "AA_Thr", "AA_Trp")
AAMap[["AcetylCoA"]] <- c("AA_IL", "AA_Leu", "AA_Lys", "AA_Trp")

# We can't do Fumarate, because it only has two AA's.
#AAMap[["Fumarate"]] <- c("AA_Phe", "AA_Tyr")
AAMap[["SuccinylCoA"]] <- c("AA_IL", "AA_METH", "AA_Val")

# All AA
AAMap[["AA"]] <- c("AA_ASP","AA_GLU","AA_ASN","AA_SER","AA_GLN","AA_GLY","AA_THR","AA_CIT","AA_ARG","AA_TAU","AA_ALA","AA_GABA","AA_TYR","AA_METH","AA_VAL","AA_PHE","AA_IL","AA_LEU","AA_ORN","AA_LYS")

# Organic Acids
AAMap[["OrganicAcids"]] <- c("OA_Lactic", "OA_Pyruvic", "OA_Succinic", "OA_Fumaric", "OA_a-Ketoglutaric", "OA_Citric", "OA_Malic")

inFile <- paste(DataDir, "CardiacArrest1023MetabolomicsInitialToLarry.csv", sep="")
Triple <- F # Are we using 3 samples: Triple = T, or just two, Triple = F

plotPCACustom <- function(mSet, fileName) {
  cls <- mSet$dataSet$cls
  pcaDF <- data.frame(Treatment=cls, mSet$analSet$pca$x[,c("PC1", "PC2")], label=rownames(mSet$analSet$pca$x))
  
  gp <- ggplot(pcaDF, aes(x=PC1, y=PC2, colour=Treatment, fill=Treatment)) + 
    geom_label(aes(label=label), colour="white", fontface="bold", family="Helvetica", size=2.0, show.legend=F) +
    stat_ellipse() +
    theme_bw() + 
    scale_colour_manual(values=c("red", "blue", "orange"), aesthetics=c("fill", "colour"))
  
  ggsave(gp, file=fileName)
}

# Initialize the mSet objects and type of data
initMSet <- function(infile, plotNorm=FALSE) {
  # Initialize the data object, we are dealing with concentration data, and want to do statistics on it
  mSet <- InitDataObjects("conc", "stat", FALSE) %>%
    Read.TextData(inFile, "rowu", "disc") %>%   # Read in the data from file
    SanityCheckData() %>%
    ReplaceMin() %>%
    PreparePrenormData() %>%
    
    # No row normalization needed, since the samples are normalized to weight
    Normalization("NULL", "CrNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
  
  # Should we plot the normalization results
  if(plotNorm) {
    mSet <- PlotNormSummary(mSet, "norm", "pdf", width=NA)
    mSet <- PlotSampleNormSummary(mSet, "snorm", "pdf", width=NA)
  }
  
  mSet
}

# Takes an mSet and performs the various statistical computations on it that we need for analysis.
analyzeMSet <- function(mSet, suffix="", DPI="auto", combined) {
  addSuffix <- function(...) paste0(..., suffix)
  
  # First perform PCA, and create the plots
  mSet <- PCA.Anal(mSet) 
  
  save(mSet, file=paste0("mSet", combined, ".Rdata"))
  
  # Save the pca plot
  plotPCACustom(mSet, paste0("pcaScores_", combined, ".pdf"))
  
  # Continue analysis
  mSet <- mSet %>%
    PlotPCAPairSummary(imgName=addSuffix("pca_pair"), format="pdf", dpi=DPI, width=NA, pc.num=2) %>%
    PlotPCAScree(imgName=addSuffix("pca_scree"), format="pdf", width=NA, dpi=DPI, scree.num=2) %>%
    PlotPCALoading(imgName=addSuffix("pca_loading"), format="pdf", width=NA, dpi=DPI, inx1=1, inx2=2) %>%
    PlotPCABiplot(imgName=addSuffix("pca_biplot"), format="pdf", width=NA, inx1=1, inx2=2, dpi=DPI) %>%
    PlotPCA3DScoreImg(imgName=addSuffix("pca_score3d"), format="pdf", dpi=DPI, width=NA, inx1=1, inx2=2, inx3=3, angl=40)
  
  # Now perform anova analysis
  mSet <- ANOVA.Anal(mSetObj=mSet, nonpar=T, thresh=0.05, post.hoc="tukey") %>%
    plotAnova(addSuffix("aov"), "pdf", dpi=400, width=NA) 
  
  # Random forest analysis
  set.seed(20190901)
  mSet <- RF.Anal(mSet, treeNum=500, tryNum=7, randomOn=1) %>%
    PlotRF.Classify(addSuffix("rf_cls"), "pdf", dpi=DPI, width=NA) %>%
    PlotRF.VIP(addSuffix("rf_imp"), "pdf", dpi=DPI, width=NA) %>%
    PlotRF.Outlier(addSuffix("rf_outlier"), "pdf", dpi=DPI, width=NA)
  
  # Write the RF analyis data to file
  sink(paste0(addSuffix("rf_res"), ".txt"))
  print(mSet$anal$rf)
  sink()
  
  mSet
}

# Replacement for PlotPCAPairSummary to plot the PCA components
plotPCA <- function(mSetObj, imgName, format, dpi, width) {
  # See, https://github.com/xia-lab/MetaboAnalystR/blob/dbbe0b1ebc34c5c20e2407d8a7c132fc648ec403/R/stats_chemometrics.R
  
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="")
  w <- ifelse(is.na(width), 10, ifelse(width == 0, 8, width))
  
  mSetObj$imgSet$pca.pair <- imgName
  
  h <- w
  if(mSet$dataSet$cls.type != "disc") stop("Do not know how to handle non disc cls type")
  
  return(mSetObj)
}

# Perform analysis with HD-CPR and SHAM both combined and separate.
analyzeDuo <- function(mSet, suffix="", computePairs) {
  # First perform the analysis with all the classes
  analyzeMSet(mSet, suffix, combined="1")
  
  # Now combine HD-CPR and SHAM
  mSetCombined <- mSet
  
  if(computePairs) {
    # adjust the classes
    newClasses <- as.character(mSet$dataSet$cls)
    newClasses[newClasses == "HD-CPR" | newClasses == "SHAM"] <- "HD-CPRSHAM"
    newClasses <- factor(newClasses)
    
    # Now change the appropriate variables
    mSetCombined$dataSet$cls <- newClasses
    mSetCombined$dataSet$orig.cls <- newClasses
    mSetCombined$dataSet$prenorm.cls <- newClasses
    mSetCombined$dataSet$proc.cls <- newClasses
    
    # Now re-run the analyses
    analyzeMSet(mSetCombined, suffix=paste0("_HD-CPRsham_", suffix), combined="2")
  
    # Now combine Std-CPR and SHAM
    mSetCombined <- mSet
    
    # adjust the classes
    newClasses <- as.character(mSet$dataSet$cls)
    newClasses[newClasses == "Std-CPR" | newClasses == "SHAM"] <- "Std-CPRSHAM"
    newClasses <- factor(newClasses)
    
    # Now change the appropriate variables
    mSetCombined$dataSet$cls <- newClasses
    mSetCombined$dataSet$orig.cls <- newClasses
    mSetCombined$dataSet$prenorm.cls <- newClasses
    mSetCombined$dataSet$proc.cls <- newClasses
    
    # Now re-run the analyses
    analyzeMSet(mSetCombined, suffix=paste0("_Std-CPRsham_", suffix), combined="3")
  }
}

# Do all the combinations, in essence do all the analyses
analyze <- function() {
  setwd(ResultsDir)
  
  # Load and initialize the data
  mSet <- initMSet(inFile, TRUE)
  
  # First compute for everything
  newDir <- paste0(ResultsDir, "All")
  unlink(newDir, recursive=T)
  dir.create(newDir)
  setwd(newDir)
  
  analyzeDuo(mSet, suffix="_all", Triple)
  
  # Now go through all the combinations
  for(aaType in names(AAMap)) {
    # Determine which columns we need
    colNames <- toupper(colnames(mSet$dataSet$orig)) %in% toupper(AAMap[[aaType]])
    colNames[colnames(mSet$dataSet$orig) %in% c("Samples", "Class")] <- TRUE
    
    mSetComb <- mSet
    mSetComb$dataSet <- within(mSetComb$dataSet, {
                               norm <- norm[,colNames]
                               row.norm <- row.norm[,colNames]
                               proc <- proc[,colNames]
                               preproc <- preproc[,colNames]
                               orig <- orig[,colNames]
                               cmpd <- colNames
    })
    
    # Change output directory
    newDir <- paste0(ResultsDir, aaType)
    unlink(newDir, recursive=T)
    dir.create(newDir)
    setwd(newDir)
    analyzeDuo(mSetComb, suffix=paste0("_", aaType), Triple)
  }
}

# My version of PlotANOVA which adds the compund names.
plotAnova <- function(mSetObj=NA, imgName, format="png", dpi=72, width=NA){
  
  mSetObj <- mSetObj
  
  lod <- mSetObj$analSet$aov$p.log;
  imgName = paste(imgName, "dpi", dpi, ".", format, sep="");
  
  if(is.na(width)){
    w <- 9;
  }else if(width == 0){
    w <- 7;
  }else{
    w <- width;
  }
  h <- w*6/9;
  
  mSetObj$imgSet$anova <- imgName;
  
  # convert the lod to a dataframe
  lodDF <- data.frame(x=sub("AA_", "", names(lod)), y=lod, col=ifelse(mSetObj$analSet$aov$inx.imp,"red","green"))
  
  gp <- ggplot(data=lodDF, aes(x=x, y=y, col=col)) + 
    geom_point(size=2.5) + 
    theme_bw() + 
    geom_hline(yintercept=mSetObj$analSet$aov$thresh, linetype=2) +
    scale_colour_manual(values=c("black", "orange")) +
    theme(axis.text.x=element_text(angle=90, hjust=1)) +
    theme(legend.position="none") +
    xlab("Metabolites") + ylab("-log10(p-value)")
  
  ggsave(filename=imgName, plot=gp, units="in", dpi=dpi, width=w, height=h)
  
  return(mSetObj)
}

# Run everything
analyze()

