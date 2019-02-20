# [ADMIXR] #####################################################################

# Package: admixr
# Title: Various Functions For Admixture-related Analyses of SNP data
# Version: 0.0.0.9000
# Author: Martin Sikora (martin.sikora@snm.ku.dk)
# Description: Various functions for admixture-related analyses of SNP data
# Depends: R (>= 3.1.2)
# License: GPL (>= 3)
# LazyData: true
# Suggests: rARPACK
# Imports: ggplot2,dplyr,grid
# Source Code: https://github.com/martinsikora/admixr

# Adapted by: Pavel Salazar-Fernandez
# Developed at: PRBB - Spain / LANGEBIO - Mexico
# Last Edit: May 26 2016

# Requirements:
# - Genotypes: RAW file (use --recode A in PLINK to generate this file)
# - SNP info: BIM file
# - Populations info: FAM file (Use FAMID to tag population membership)
# - Admixture: Q file. Optional: P file.
# - R packages: admixr, ggplot2, data.table, ape. Optional: foreach, doMC.

# Pipeline:
# 1. Reads input data and transforms it
# 2. Generates PCA plots
# 3. Gets allele frequencies
# 4. Calculates Fst for population pairs
# 5. Performs outgroup f3 statistics calculations
# 6. Obtains D statistics
# 7. Makes Admixture plots and projections

# Procedure: 
# 1. Set the working and output directory.
# 2. Select the analyses to be performed.
# 3. Declare the names of the input files.
# 4. Adjust the parameters for each analysis.
# 5. Source the script. Messages will inform the progress.
# 6. Once finished, all the output files will be in the working directory.

# Warning:
# This script has been adapted and automated. In case of errors, please check 
# your input files and parameters, and run lines and sections manually.

#<START> #######################################################################

#[1. PREPARATIONS] =============================================================

#<INPUT> -----------------------------------------------------------------------

# Working and output directory

setwd("PATH/")

# Select modules to be run:
# 1. Data Input (execute in the first run, skip for re-runs)
# 2. PCA plot (2.1 = Full PCA; 2.2 = Fast PCA)
# 3. Allele frequencies (Required for [6])
# 4. Fst 
# 5. f3 statistics (5.1 = Outgroup; 5.2 = Paired; 5.3 = Grouped)
# 6. D statistics (Requires Module [3])
# 7. Admixture (7.1 = Plot; 7.2 = Projections)

run.modules <- c(1,2.2,3,4,5.1,6)


# 1. Genotyping Data (Required for all modules)

raw.file <- ".raw"
bim.file <- ".bim"
fam.file <- ".fam"

# 2,3. PCA & Allele Frequencies

## All required data is obtained from Module 1

# 4. Fst: Outgroup for the NJ tree
nj.outgroup <- "XYZ"

# 5. f3 Statistics

## 5.1 Outgroup f3

# Outgroup
idxO <-"XYZ"

# Test
idxX <- "XYZ"

## 5.2 Paired
popd1 <- "ABC"
popd2 <- "UVW"

## 5.3 Grouped
popd <- "ABC"
groups <- c(rep("outgroup",1),rep("derived",1))



# 6. D Statistic (Outgroup,Test)(X,POP2)
ds.outgroup <- "ABC"
ds.test <- "XYZ"
ds.pop2 <- "UVW"

# 7. ADMIXTURE

## 7.1 Admixture Plot
qfile <- ".Q"
fam <- read.table(".fam", stringsAsFactors = FALSE)

## 7.2 Admixture Projection

### Reference panel with omitted populations
proj.pops <- "XYZ"
proj.qfile <- ".Q"
proj.fam <-  ".fam"
proj.bim <- ".bim"
proj.pfile <- ".P"

#</INPUT>

#<PREPARATIONS>

message(">>> INITIALIZING: ADMIXR+ ")

#<LIBRARIES>
message("> Loading libraries...")
# Load required packages
library("admixr")
library("ggplot2")
library("data.table")

# Optional: for parallel computations
library("foreach")
library("doMC")

#</LIBRARIES>

if(1 %in% run.modules){ 
  message("> STARTING: Module [1. PREPARATIONS]")
  #<DATA>
  message("> Loading data...")
  # Load RAW file
  gt <- data.frame(fread(raw.file, header=T, 
                         drop = c("FID","PAT","MAT","SEX","PHENOTYPE"))
                   ,row.names = 1)
  colnames(gt) <- gsub('^{1}.', '',colnames(gt))
  colnames(gt) <- gsub('.{2}$', '',colnames(gt))
  gt <- t(gt)
  
  # Load BIM file
  snpInfo  <- read.table(bim.file, header=F,
                         colClasses = c(NA,NA,"NULL",NA,"NULL","NULL"), 
                         col.names = c("chromosome","posId","","position","",""),
                         stringsAsFactors = FALSE)
  
  # Load FAM file
  sampleInfo <- read.table(fam.file, colClasses = c(NA,NA,rep("NULL",4)), 
                           header=F, col.names =  c("popId","sampleId",rep("",4)),
                           stringsAsFactors = FALSE)
  popInfo <- data.frame(unique(sampleInfo$popId), stringsAsFactors = FALSE,
                        rainbow(length(unique(sampleInfo$popId))))
  colnames(popInfo) = c("popId","color")
  popInfo$color <- gsub('.{2}$', '',popInfo$color)
  message("> Data loaded correctly")
  #</DATA>
} else {
  message("> SKIPPED: Module [1. PREPARATIONS]")
}


message(">> STARTING ANALYSIS MODULES")

#[/PREPARATIONS]

#[2. PCA] ======================================================================

if (any(c(2.1,2.2) %in% run.modules)){
  message("> STARTING: Module [2. PCA]")
  # Do PCA on the genotype matrix
  res <- list()
  
  message("> PROCESSING: Module [2. PCA]")
  # Full PCA
  if(2.1 %in% run.modules) { res$full <- getPcaGT(gt) }
  
  # Approximate PCA (faster)
  if(2.2 %in% run.modules) { res$apx <- getPcaGTFast(gt) }
  
  #<PLOT>

  
  idx <- matrix(paste("PC", 1:10, sep = ""), nrow = 2)
  idxC <- popInfo$color
  names(idxC) <- popInfo$popId
  
  for(k in names(res)){
      d <- as.data.frame(res[[k]]$summary$pca)
      d$sampleId <- rownames(d)
      d$popId <- sampleInfo$popId[match(d$sampleId, sampleInfo$sampleId)]
      pdf(paste(k, "pca.plot.pdf", sep = "."), width = 4, height = 4)
      for(i in 1:ncol(idx)){
          p <- ggplot(d, aes_string(x = idx[1, i], y = idx[2, i], 
                                    colour = "popId", fill = "popId"))
          print(p + geom_text(aes(label = gsub("pop", "", popId)), size = 2, 
                              alpha = 0.5) + 
                  scale_colour_manual(name = "Population", values = idxC) + 
                  scale_fill_manual(name = "Population", values = idxC) + 
                  labs(x = idx[1, i], y = idx[2, i]) + theme_bw() + 
                  theme(panel.grid.minor=element_blank(), 
                        panel.grid.major=element_blank(), 
                        plot.title = element_text(size = 10), 
                        legend.position = "none"))
      }
      dev.off()
  }
  message("> OUTPUT: Module [2. PCA] > pca.plot.pdf")
  message("> FINISHED: Module [2. PCA]")
  #</PLOT>
} else {
  message("> SKIPPED: Module [2. PCA]")
}

#[/PCA]

#[3. ALLELE FREQUENCIES] =======================================================

if (3 %in% run.modules){
  message("> STARTING: MODULE [3. ALLELE FREQUENCIES]")
  # Get allele frequencies
  ## Index vector for population labels of individuals
  idxP <- sampleInfo$popId[match(colnames(gt), sampleInfo$sampleId)] 
  ## Allele counts for both alleles in each population
  counts <- getAlleleCountsGT(gt, idxP)
  ## Frequency of allele 2 in each population
  freq <- counts[,,2] / (counts[,,1] + counts[,,2]) 
  message("> FINISHED: MODULE [3. ALLELE FREQUENCIES]")
} else {
  message("> SKIPPED: Module [3. ALLELE FREQUENCIES]")
  message("! WARNING: Module [3] required for Module [4] and [6]")
}

#[/ALLELEFREQS]

#[4. FST] ======================================================================

if (4 %in% run.modules){
  message("> STARTING: Module [4. FST]")  
  # Calculate FST for all pairs of populations
  ## Index vector for population labels of individuals
  idxP <- sampleInfo$popId[match(colnames(gt), sampleInfo$sampleId)]
  ## Index matrix of all pairwise comparisons
  idxM <- t(combn(popInfo$popId, 2)) 
  fst <- as.data.frame(idxM, stringsAsFactors = FALSE)
  colnames(fst) <- c("p1", "p2")
  fst$fst <- NA
  ## Second data frame for calculation from allele counts
  fst1 <- fst 
  message("> PROCESSING: Module [4. FST]")
  # Weir & Cockerham 1984, for diploid genotypes
  for(i in 1:nrow(idxM)){
      cat(i, "\r")
      idx <- idxP %in% idxM[i,]
      fst$fst[i] <- getFstGT(gt[, idx], idxP[idx], region = TRUE)
  }
  
  # Weir & Hill 2002, from allele frequencies assuming no inbreeding
  for(i in 1:nrow(idxM)){
      cat(i, "\r")
      fst1$fst[i] <- getFstAlleleCounts(counts[, idxM[i,],], region = TRUE) #
  }
  
  # Compare the estimators
  pdf("fst.comparison.plot.pdf", width = 5, height = 5)
  plot(fst$fst, fst1$fst, xlab = "WC1984", ylab = "WH2002", asp = 1, type = "n")
  abline(0, 1, col = "grey")
  points(fst$fst, fst1$fst, pch = 21, col = "#00000055", bg = "#00000055")
  dev.off()
  message("> OUTPUT: Module [4. FST] > fst.comparison.pdf")
  
  # Convert to matrix and make NJ tree from WC1984 estimator results
  fst1 <- fst[, c(2:1, 3)]
  colnames(fst1) <- colnames(fst)
  fst1 <- rbind(fst, fst1)
  fst1$p1 <- factor(fst1$p1, levels = popInfo$popId)
  fst1$p2 <- factor(fst1$p2, levels = popInfo$popId)
  d <- reshape2::acast(fst1, p1 ~ p2, value.var = "fst")
  d1 <- ape::bionj(d) ## bionj algorithm from ape package
  
  #<PLOT>
  pdf("fst.nj.plot.pdf", width = 8, height = 7)
  plot(ape::root(d1, outgroup = nj.outgroup), font = 1, tip.color = idxC[d1$tip.label], 
       cex = 1)
  plot(d1, type = "unrooted", lab4ut = "axial", font = 1, 
       tip.color = idxC[d1$tip.label], cex = 1)
  dev.off()
  message("> OUTPUT: Module [4. FST] > fst.nj.plot.pdf")
  
  #</PLOT>
  
  #<PLOT>
  pal <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
  pdf("fst.matrix.plot.pdf", width = 6, height = 5)
  p <- ggplot(fst1, aes(x = p1, y = p2, fill = fst))
  print(p + geom_tile() + geom_text(aes(label = formatC(fst, digits = 2, 
                                                        format = "f")), size = 1) 
        + scale_fill_gradientn(name = expression(F[ST]), 
                               colours = pal, limits = c(0, max(fst$fst))) + 
          coord_equal() + xlab("Population 1") + ylab("Population 2") + theme_bw()
        + theme(panel.grid.minor=element_blank(), 
                panel.grid.major = element_blank(), 
                strip.background = element_rect(fill = "white", colour = NA), 
                strip.text = element_text(size = 2), 
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                                           color = idxC), 
                axis.text.y = element_text(color = idxC), 
                axis.ticks.length = grid::unit(0, "cm")))
  dev.off()
  message("> OUTPUT: Module [4.FST] > fst.matrix.plot.pdf")
  #</PLOT>
  message("> COMPLETED: Module [4. FST]")
} else {
  message("> SKIPPED: Module [4. FST]")
}

#[/FST]

#[5. F3 STATISTICS] ============================================================

if (any(c(5.1,5.2,5.3) %in% run.modules)){
  message("> STARTING: Module [5. F3 STATISTICS]")
  # [5.1 Outgroup f3 Statistics] -----------------------------------------------
  if(5.1 %in% run.modules){
    message("> PROCESSING: Module [5.1 OUTGROUP F3 STATISTICS]")
    # Outgroup f3 Statistic
    blockIdx <- getChromosomeBlocks(snpInfo$chromosome, snpInfo$position, 5e6)
    
    f3 <- list()
    ## Run 4 threads in parallel for f3 computation
    registerDoMC(4) 
    for(x in idxX){
      idxP1 <- popInfo$popId[!(popInfo$popId %in% c(idxO, x))]
      f3[[x]] <- foreach(p1 = idxP1, .combine = "rbind") %dopar% {
        idxM <- cbind(idxO, x, p1)
        r <- doF3Test(counts, idxM, blockIdx)
        r
      }
    }
    f3 <- do.call("rbind", f3)
    
    #<PLOT>
    pdf("f3.outgroup.plot.pdf", width = 4, height = 3.5)
    for(x in idxX){
      d <- f3[f3$p2 == x,]
      print(plotFStat(d$f3, d$se, d$p3, d$p3, idxC, intercept = max(d$f3), 
                      orderValue = TRUE, showLegend = FALSE) + 
              ylab(expression(f[3])) + ggtitle(x)) 
      ## add ylab and title to plot object returned from function
    }
    dev.off()
    message("> OUTPUT: Module [5. F3 STATISTICS] > f3.outgroup.plot.pdf")
    #</PLOT>
  }
  
  # [5.2 Paired f3 Statistics] --------------------------------------------------
  if(5.2 %in% run.modules){  
    # Plot pairwise f3 for popd1 / popd2
    message("> PROCESSING: Module [5.2 PAIRED F3 STATISTICS]")
    d1 <- f3[f3$p2 == popd1 & f3$p3 != popd2,]
    d2 <- f3[f3$p2 == popd2 & f3$p3 != popd1,]
    
    pdf("f3.outgroup.pairs.plot.pdf", width = 4, height = 4)
    print(plotFPairs(d1$f3, d2$f3, d1$se, d2$se, d1$p3, idxC, z = 1, shape = 21, 
                     alpha = 0.9, showLegend = FALSE) + xlab(popd1) + ylab(popd2)) 
    ## add axis labels to plot object returned from function
    print(plotFPairs(d1$f3, d2$f3, d1$se, d2$se, d1$p3, label = d1$p3, idxC, z = 1, 
                     size = 0, alpha = 0.9, showLegend = FALSE) + xlab(popd1) + 
            ylab(popd2)) ## add axis labels to plot object returned from function
    dev.off()
    
    message("> OUTPUT: Module [5. F3 STATISTICS] > f3.outgroup.pairs.plot.pdf")
  }

  # [5.3 Grouped f3 Statistics] ------------------------------------------------
  if(5.3 %in% run.modules){ 
    # Plot popd results grouped by outgroup / derived population
    message("> PROCESSING: Module [5.3 GROUPED F3 STATISTICS]")
    d <- f3[f3$p2 == popd,]
    
    idxC1 <- c(outgroup = "blue", derived = "red")
    pdf("f3.outgroup.grouped.plot.pdf", width = 4, height = 3.5)
    print(plotFStat(d$f3, d$se, d$p3, groups, idxC1, intercept = max(d$f3), 
                    grouped = TRUE, showLegend = FALSE) + ylab(expression(f[3])) 
          + ggtitle(popd)) 
    ## add ylab and title to plot object returned from function
    dev.off()
    
    message("> OUTPUT: Module [5. F3 STATISTICS] > f3.outgroup.grouped.plot.pdf")
  }
  
  message("> COMPLETED: Module [5. F3 STATISTICS]")
} else {
  message("> SKIPPED: Module [5. F3 STATISTICS]")
}
#[/F3STATS]

#[/F3STATS]

#[6. D STATISTICS] ==============================================================

if (6 %in% run.modules){
  message("> STARTING: Module [6. D STATISTICS]")
  idxM <- cbind(ds.outgroup, ds.test, 
                popInfo$popId[!(popInfo$popId %in% 
                                  c(ds.outgroup,ds.test,ds.pop2))], ds.pop2)
  registerDoMC(4)
  d <- foreach(i = 1:nrow(idxM), .combine = "rbind") %dopar% {
      r <- doDTest(freq, idxM[i, , drop = FALSE], blockIdx)
      r
  }
  
  #<PLOT>
  pdf("D.plot.pdf", width = 4, height = 3)
  print(plotFStat(d$D, d$se, d$p3, d$p3, idxC, horizontal = TRUE, 
                  showLegend = FALSE) + xlab("Population") + ylab("D") +
          ggtitle(paste("D(",ds.outgroup,",",ds.test,")(X,",ds.pop2,")",
                        sep=""))) 
          ## Add xlab and title to plot object returned from function
  dev.off()
  message("> OUTPUT: Module [6. D STATISTICS] > D.plot.pdf")
  #</PLOT>
  message("> COMPLETED: Module [6. D STATISTICS]")
} else {
  message("> SKIPPED: Module [6. D STATISTICS]")
}

#[/DSTATS]

#[7. ADMIXTURE] ================================================================
if (any(c(7.1,7.2) %in% run.modules)){
  message("> STARTING: Module [7. ADMIXTURE]")
#[7.1 Admixture Plot] ----------------------------------------------------------
  if (7.1 %in% run.modules){
    message("> PROCESSING: Module [7.1 ADMIXTURE PLOT]")
    qcols <- count.fields(qfile, sep = " ")[1]
    qMatrix <- matrix(scan(qfile), ncol = qcols, byrow = TRUE)
    rownames(qMatrix) <- fam$V2
    
    #<PLOT>
    d <- reshape2::melt(qMatrix)
    colnames(d) <- c("sampleId", "k", "value")
    d$sampleId <- factor(d$sampleId, levels = unique(d$sampleId))
    d$popId <- sampleInfo$popId[match(d$sampleId, sampleInfo$sampleId)]
    d$popId <- factor(d$popId, levels = unique(d$popId))
    idxCLabel <- popInfo$color[match(levels(d$popId), popInfo$popId)]
    idxCCluster <- rainbow(qcols)
    pdf("admixture.full.plot.pdf", width = 8, height = 2.5)
    plotAdmixture(sampleId = d$sampleId, popId = d$popId, k = d$k, value = d$value, 
                  colors = idxCCluster, labColors = idxCLabel, alpha = 1, width = 1,
                  showLegend = TRUE, rot = 90)
    dev.off()
    message("> OUTPUT: Module [7.1 ADMIXTURE PLOT] > admixture.full.plot.pdf")
    #</PLOT>
  }

#[7.2 Admixture Projection ] ---------------------------------------------------
  if (7.2 %in% run.modules){
    message("> PROCESSING: Module [7.2 ADMIXTURE PROJECTING]")
    projqcols <- count.fields(proj.qfile, sep = " ")[1]
    qMatrix <- matrix(scan(projfile), ncol = projqcols, byrow = TRUE)
    fam <- read.table(proj.fam, stringsAsFactors = FALSE)
    rownames(qMatrix) <- fam$V2
    
    # Prepare data for individuals to be projected
    idxS <- sampleInfo$sampleId[sampleInfo$popId %in% proj.pops]
    bim <- read.table(proj.bim, stringsAsFactors = FALSE)
    idx <- match(bim$V2, snpInfo$posId) ## find indices of SNPs used in ref panel
    gt1 <- gt[idx, idxS] ## subset of genotype matrix for individuals to project
                         ## and ref panel SNPs
    
    idx <- bim$V5 == 1 ## find SNPs where alleles are switched in recoded subset 
                       ## plink file of ref panel; genotype matrix has copies of 
                       ## minor allele, which is coded as allele "2" in the file
    gt1[idx,] <- 2 - gt1[idx,]
    gt1 <- 2 - gt1 ## final step is to reverse all snps for estimation, to match 
                   ## coding of P matrix from admixture (has major allele instead
                   ## of minor allele coding)
    
    # Do projection
    projpcols <- count.fields(proj.pfile, sep = " ")[1]
    pMatrix <- matrix(scan(proj.pfile), ncol = projpcols, byrow = TRUE)
    registerDoMC(4) ## run 4 threads in parallel
    qMProj <- foreach(s = idxS, .combine = "rbind") %dopar% {
        cat(s, "\r")
        doAdmixtureProjection(gt1[, s], pMatrix)
    }
    rownames(qMProj) <- idxS
    qMatrix <- rbind(qMatrix, qMProj)
    qMatrix <- qMatrix[sampleInfo$sampleId,]
    
    #<PLOT>
    d <- reshape2::melt(qMatrix)
    colnames(d) <- c("sampleId", "k", "value")
    d$sampleId <- factor(d$sampleId, levels = unique(d$sampleId))
    d$popId <- sampleInfo$popId[match(d$sampleId, sampleInfo$sampleId)]
    d$popId <- factor(d$popId, levels = unique(d$popId))
    idxCLabel <- popInfo$color[match(levels(d$popId), popInfo$popId)]
    idxCCluster <- rainbow(projpcols)
    pdf("admixture.projected.plot.pdf", width = 8, height = 2.5)
    plotAdmixture(sampleId = d$sampleId, popId = d$popId, k = d$k, 
                  value = d$value, colors = idxCCluster, labColors = idxCLabel, 
                  alpha = 1, width = 1, showLegend = TRUE, rot = 90)
    dev.off()
    message("> OUTPUT: Module [7.2 ADMIXTURE PROJECTING] > admixture.projected.plot.pdf")
    #</PLOT>
  }
  message("> COMPLETED: Module [7. ADMIXTURE]")
} else {
  message("> SKIPPED: Module [7. ADMIXTURE]")
}

#[/ADMIXTURE]

message(">>> FINISHED: ADMIXR+")

#<END> #########################################################################