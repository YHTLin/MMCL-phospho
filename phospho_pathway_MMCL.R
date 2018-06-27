########################################################################
### KINASE ACTIVITY INFERENCE ON MASS SPECTROMETRY-BASED PHOSPHO DATA
########################################################################

############################################################
# Step 1a - Data Filtering Functions
############################################################
## Group filtering to remove -Inf values (sourced from pulsed_SILAC_ozlem.R)
filter_valids = function(df, conditions, min_count, at_least_one = FALSE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  
  log2.names = grep("^LOG2", names(df), value = TRUE)   # Extract LOG2 column names
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}



############################################################
# Step 1b - Data Filtering Functions
############################################################
## NORMALIZATION METHOD 2: Global median normalization for each sample (center median at 0)
median_centering = function(df, use_keep = TRUE) {
  # df = data frame containing LOG2 columns for normalization
  # use_keep = filter rows using KEEP column prior to computing the median for normalization
  LOG2.names = grep("^LOG2", names(df), value = TRUE)
  
  if (use_keep) {
    df[, LOG2.names] = lapply(LOG2.names, 
                              function(x) {
                                LOG2 = df[df$KEEP, x]
                                LOG2[!is.finite(LOG2)] = NA   # Exclude 0 intensity values from median calculation
                                gMedian = median(LOG2, na.rm = TRUE)
                                tmp = df[, x] - gMedian
                                tmp[!is.finite(tmp)] = NA
                                tmp
                              }
    )
  } else {
    df[, LOG2.names] = lapply(LOG2.names, 
                              function(x) {
                                LOG2 = df[, x]
                                LOG2[!is.finite(LOG2)] = NA   # Exclude 0 intensity values from median calculation
                                gMedian = median(LOG2, na.rm = TRUE)
                                tmp = df[, x] - gMedian
                                tmp[!is.finite(tmp)] = NA
                                tmp
                              }
    )
  }
  return(df)
}


## NORMALIZATION METHOD 1: Quantile normalization for each cell line (alternative to global median normalization)
quantile_norm = function(df, conditions, use_keep = TRUE) {
  # df = data frame containing LOG2 columns for normalization
  # conditions = a character vector dictating the grouping
  # use_keep = filter rows using KEEP column prior to normalization
  require(dplyr)
  require(limma)
  
  log2.names = grep("^LOG2", names(df), value = TRUE)   # Extract LOG2 columns
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  if (use_keep) df = df[df$KEEP, ]
  LOG2_df = as.matrix(df[log2.names])
  LOG2_df[!is.finite(LOG2_df)] = NA
  
  for(group in cond.names) LOG2_df[, group] = normalizeQuantiles(LOG2_df[, group])
  
  df[log2.names] = LOG2_df[, log2.names]
  
  return(df)
}



############################################################
# Step 1c - Data Imputation Functions
############################################################
## Impute missing data METHOD 1: Downshift and tighten the sampling distribution of missing values
# Does not remove KEEP = FALSE in output
impute_data = function(df, width = 0.3, downshift = 1.8, use_keep = TRUE) {
  # df = data frame containing filtered 
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  # use_keep = filter rows using KEEP column prior to imputation
  
  LOG2.names = grep("^LOG2", names(df), value = TRUE)
  impute.names = sub("^LOG2", "impute", LOG2.names)
  
  # Create new columns indicating whether the values are imputed
  df[impute.names] = lapply(LOG2.names, function(x) !is.finite(df[, x]))
  
  # Imputation
  set.seed(1)
  df[LOG2.names] = lapply(LOG2.names,
                          function(x) {
                            temp = df[[x]]
                            temp[!is.finite(temp)] = NA
                            
                            if (use_keep) {  # calculate sd and mean from filtered rows
                              temp.sd = width * sd(temp[df$KEEP], na.rm = TRUE)   # shrink sd width
                              temp.mean = mean(temp[df$KEEP], na.rm = TRUE) - 
                                downshift * sd(temp[df$KEEP], na.rm = TRUE)   # shift mean of imputed values
                            } else {  # calculate sd and mean from all rows
                              temp.sd = width * sd(temp, na.rm = TRUE)   # shrink sd width
                              temp.mean = mean(temp, na.rm = TRUE) - 
                                downshift * sd(temp, na.rm = TRUE)   # shift mean of imputed values
                            }
                              
                            n.missing = sum(is.na(temp))
                            temp[is.na(temp)] = rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                            return(temp)
                          })
  return(df)
}


## Impute missing data METHOD 2: Hybrid missing data imputation (sourced from Jae_surfaceome_analysis.R)
# Removes KEEP = FALSE in output
hybrid_impute = function(df, conditions, use_keep = TRUE) {
  # df = data frame containing filtered 
  # conditions = a character vector dictating the grouping
  # use_keep = filter rows using KEEP column prior to imputation (WILL AFFECT IMPUTATION OUTCOME IF DATA NOT FILTERED ALREADY)
  require(imputeLCMD)
  require(dplyr)
  
  # Apply KEEP filter
  if (use_keep) df = df[df$KEEP, ]
  
  # Group column names by condition
  log2DF = dplyr::select(df, starts_with("LOG2"))   # Extract LOG2 dataframe
  DF = dplyr::select(df, -starts_with("LOG2"))
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, names(log2DF), value = TRUE, perl = TRUE))
  
  # Create new columns indicating whether the values are imputed
  impute.names = sub("^LOG2", "IMP", unlist(cond.names))
  DF[impute.names] = lapply(unlist(cond.names), function(x) !is.finite(log2DF[[x]]))
  
  # Imputation
  set.seed(1)
  imputeDF = lapply(cond.names,   # Impute each group separately
                    function(x) {
                      tempDF = as.matrix(log2DF[x])
                      modelS = model.Selector(tempDF)
                      impute.MAR.MNAR(tempDF, modelS, method.MAR = "MLE", method.MNAR = "MinProb")
                    })
  imputeDF = do.call(cbind, imputeDF)   # Combine a list of data frames into a single data frame
  
  return(cbind(DF, imputeDF))
}



############################################################
# Step 1d - Data Imputation Function
############################################################
## Average biological replicates
average_reps = function(df, samples) {
  # df = data frame containing imputed data
  # samples = character vector containing sample names for averaging over LOG2 columns
  require(dplyr)
  
  log.names = grep("^LOG2", names(df), value = TRUE)
  group.names = lapply(samples, function(x) grep(x, log.names, value = TRUE, perl = TRUE))
  avg.names = sapply(group.names, function(x) sub("_.*", "", x[1]))

  df[avg.names] = sapply(group.names, function(x) rowMeans(df[x]))
  
  return(df)
}



############################################################
# Step 2 - Data clean-up and pre-processing
############################################################
## Set working directory
setwd("C:/Users/Tony Lin/Desktop/Wiita_Lab/Projects/Proteomics_project/Phosphoproteomics_MMCL/2 - MaxQuant Analysis")


## Read in phospho raw file
raw.data = read.delim("maxquant_output/Phospho (STY)Sites.txt", 
                      header = TRUE, stringsAsFactors = FALSE)


## Remove reverse proteins, contaminations, and localization probability < 0.75
require(dplyr)
data = dplyr::filter(raw.data, Reverse != "+") %>%
  dplyr::filter(Potential.contaminant != "+") %>%
  dplyr::filter(Localization.prob >= 0.75) %>%
  dplyr::filter(Delta.score >= 8)


## Calculate log2 of intensity
Intensity.names = grep("^Intensity", names(data), value = TRUE)
LOG.names = sub("^Intensity", "LOG2", Intensity.names)
data[, LOG.names] = lapply(data[, Intensity.names], function(x) log2(x))
rm(Intensity.names, LOG.names)


## Filter columns by names and extract protein IDs
filter_cols = function(df) {
  # df = data frame containing phosphoproteomics data
  require(stringr)
  
  # Extract UniProtID from Protein column
  regex = regexpr("(?<=\\|).*(?=\\|)", df$Protein, perl = TRUE)
  out = rep(NA, nrow(df))
  out[regex != -1] = regmatches(df$Protein, regex)
  df$UniProtID = out
  
  # Extract gene name (WORKS WITH SWISS-PROT DATABASE)
  regex2 = regexpr("((?<=\\|[[:alnum:]]{6}\\|).*(?=_HUMAN)|(?<=\\|[[:alnum:]]{10}\\|).*(?=_HUMAN))", 
                   df$Protein, perl = TRUE)
  out2 = rep(NA, nrow(df))
  out2[regex2 != -1] = regmatches(df$Protein, regex2)
  df$Gene.Name = out2
  
  LOG2.names = grep("^LOG2\\..*(?<!___[0-9])$", names(df), value = TRUE, perl = TRUE)
  
  keepCols = c("UniProtID", "Gene.Name", "Position", "Amino.acid", LOG2.names)
  return(df[keepCols])
}
dataC = filter_cols(data) %>%
  filter(!is.na(UniProtID)) %>%   # Eliminate rows without valid identifiers
  dplyr::select(-LOG2.INA6_bR2) %>%    # INA6_bR2 correlate poorly with other bioreplicates so is removed from analysis
  dplyr::select(-(LOG2.U266_bR1:LOG2.U266_bR3))   # Remove U266 from analysis


## Filter rows by valid values
dataCF = filter_valids(dataC, 
                       c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34", "L363", "MM1S", "RPMI8226"),
                       min_count = rep(2, 8),
                       at_least_one = TRUE)


## Normalization
#dataCFN = quantile_norm(dataCF,
#                        c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34", "L363", "MM1S", "RPMI8226"),
#                        use_keep = TRUE)
dataCFN = median_centering(dataCF, TRUE)


## Imputation (remove KEEP = FALSE after imputation)
#dataCFNI = impute_data(dataCFN, use_keep = TRUE)
#dataCFNI = dplyr::filter(dataCFNI, KEEP)
dataCFNI = hybrid_impute(dataCFN, 
                         c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34", "L363", "MM1S", "RPMI8226"),
                         use_keep = TRUE)


## Calculate ANOVA for significance of differences in phosphosite intensity across cell lines
#testDF = select(dataCFNI, starts_with("LOG2"))
#pvalues = sapply(1:nrow(testDF), function(i) {
#  y = melt(testDF[i, ])
#  y$variable = sub("_.*", "", y$variable)
#  pvalues[i] = oneway.test(formula = value ~ variable, data = y, var.equal = F)$p.value
#})


## Average replicates
dataCFNIA = average_reps(dataCFNI, c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34",
                                     "L363", "MM1S", "RPMI8226"))


## Organize phospho data
phospho = dataCFNIA[grep("^LOG2.*(?<!_bR[0-9])$", names(dataCFNIA), value = TRUE, perl = TRUE)]

# Normalize intensities such that the mean intensity is 0 for each site (relative changes in phosphosite intensity)
phospho = as.data.frame(t(apply(as.matrix(phospho), 1, function(x) x - mean(x))))

# Annotate rows
phospho$ID = paste0(dataCFNIA[["UniProtID"]], "_",
                    dataCFNIA[["Amino.acid"]],
                    dataCFNIA[["Position"]]) # Organize phospho matrix ([UniProtID]_[Amino Acid][Position])
phospho$ID = toupper(phospho$ID)   # Raise IDs to upper case
phospho$Gene.name = dataCFNIA[["Gene.Name"]]

# Organize into list with 2 elements, value (phospho intensities) and lab (phospho site/annotation)
phospho = list(values = dplyr::select(phospho, starts_with("LOG2")),
               lab = dplyr::select(phospho, ID:Gene.name))
colnames(phospho$values) = sub("^LOG2.", "", colnames(phospho$values))


## Save workspace
#save.image("checkpoint_phospho_v2.RData")



############################################################
# Step 3a - Kinase Set Enrichment Analysis + PhosphoSitePlus
############################################################
#load("checkpoint_phospho_v2.RData")
## Load PhosphoSitePlus Kinase-Substrate table
pSitePlus = read.table("PhosphoSitePlus_Kinase_Substrate_Dataset_180205", header = TRUE,
                       colClasses = "character", skip = 3, fill = TRUE, sep = "\t")
pSitePlus = pSitePlus[pSitePlus$KIN_ORGANISM == "human" &
                        pSitePlus$SUB_ORGANISM == "human", ]   #Filter for interactions in human
# Define "Site" column as UniProtID_Psite
pSitePlus$SUB_ACC_ID = sub("-.*", "", pSitePlus$SUB_ACC_ID)   # Remove extraneous identifier in substrate UniProtID
pSitePlus$KIN_ACC_ID = sub("-.*", "", pSitePlus$KIN_ACC_ID)   # Remove extraneous identifier in kinase UniProtID
pSitePlus$Site = paste0(pSitePlus$SUB_ACC_ID, "_", pSitePlus$SUB_MOD_RSD)   # Contains 3 p-His
pSitePlus$Site = toupper(pSitePlus$Site)   # Raise IDs to upper case


## Construct adjacency matrix from pSitePlus data
adjMat_pSitePlus = matrix(0, 
                          nrow = length(unique(pSitePlus$Site)),   # Substrate
                          ncol = length(unique(pSitePlus$GENE)),   # Kinases
                          dimnames = list(row = unique(pSitePlus$Site),
                                          col = unique(pSitePlus$GENE)))
for (i in 1:nrow(pSitePlus)) {
  adjMat_pSitePlus[pSitePlus[i, "Site"], pSitePlus[i, "GENE"]] = 1   #Assign 1 if interaction exists
}


## Calculate kinase activity matrix (m samples by n kinases)
KSEA = function(phospho, pMatrix) {
  # phospho = a list of two dataframes with same number of columns (values = phospho data; lab = site information)
  # pMatrix = an adjacency matrix detailing relationship between kinase and substrate
  # Method described in Chapter 6 of Springer Cancer Systems Biology Protocol 
  store = matrix(NA, nrow = ncol(pMatrix), ncol = ncol(phospho$values),
                 dimnames = list(colnames(pMatrix), colnames(phospho$values)))   #Initialize matrix
  zScore = store   #Initialize matrix
  pValue = store   #Initialize matrix
  subMatch = store   #initialize matrix (number of matched substrates for each kinase)
  
  require(reshape2)
  for (j in 1:ncol(pMatrix)) {
    # Iterate through each kinase
    sites = rownames(pMatrix)[pMatrix[, j] == 1]   # Extract substrate ID associated with kinase
    phospho_id = phospho$values[phospho$lab$ID %in% sites, ]    # Extract log2FCs for all associated substrates
    for (i in 1:ncol(phospho$values)) {
      # Iterate through each sample
      if (nrow(phospho_id) != 0) {
        store[j, i] = mean(phospho_id[, i])   #Calculate KSEA score
        zScore[j, i] = (store[j, i] - mean(phospho$values[, i])) * sqrt(nrow(phospho_id)) / sd(phospho$values[, i])   #Calculate zscore
        subMatch[j, i] = nrow(phospho_id)   #Save number of matched substrates
      }
    }
    
    # Calculate P-value from ANOVA
    #if (nrow(phospho_id) > 1) {
    #  tmpDF = melt(phospho_id)
    #  # Welch's vs normal ANOVA
    #  #pValue[j, ] = summary(aov(formula = value ~ variable, data = tmpDF))[[1]][[1,"Pr(>F)"]]
    #  pValue[j, ] = oneway.test(formula = value ~ variable, data = tmpDF, var.equal = F)$p.value
    #}
  }
  pValue = 2*pnorm(-abs(zScore))   # 2-sided test
  
  # Correct for multiple testing
  pValue = apply(pValue, 2, function(x) p.adjust(x, method = "BH"))
  
  # Keep valid rows
  KEEP = !is.na(store[, 1])
  
  return(list(score = store[KEEP, ], 
              zScore = zScore[KEEP, ], 
              pValue = pValue[KEEP, ], 
              subMatch = subMatch[KEEP, ]))
}
kinAct_pSitePlus  = KSEA(phospho, adjMat_pSitePlus)



############################################################
# Step 3b - Kinase Set Enrichment Analysis + NetworKIN
############################################################
## Produce a table of sites for NetworKIN to predict possible kinase associations
#write.table(dataCFNIA[dataCFNIA$KEEP, c("UniProtID", "Position", "Amino.acid")],
#            file = "all-filtered-sites-quantileNorm-hybridImpute.tsv", 
#            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


## Load NetworKIN phosphosite-kinase predictions
networKIN = read.table(file = 'NetworKIN_predictions_180404.tsv', sep = '\t', header = TRUE,
                       stringsAsFactors = F, quote = "", fill = TRUE)
# Set threshold for networKIN score
networKIN = networKIN[networKIN$NetworKIN.score > 2, ]   # CHANGE CUTOFF HERE
# Organize site information
networKIN$Site = paste0(networKIN$Name, "_", networKIN$Position)   # Contains 3 p-His
networKIN$Site = toupper(networKIN$Site)   # Raise IDs to upper case


## Construct adjacency matrix from NetworKIN predictions data (sites in rows and kinase in columns)
adjMat_networKIN = matrix(0, 
                          nrow = length(unique(networKIN$Site)),
                          ncol = length(unique(networKIN$Kinase.Phosphatase.Phospho.binding.domain)),
                          dimnames = list(row = unique(networKIN$Site),
                                          col = unique(networKIN$Kinase.Phosphatase.Phospho.binding.domain)))
for (i in 1:nrow(networKIN)) {
  adjMat_networKIN[networKIN[i, "Site"], networKIN[i, "Kinase.Phosphatase.Phospho.binding.domain"]] = 1   #Assign 1 if interaction exists
}


## Calculate kinase activity matrix (m samples by n kinases)
kinAct_networKIN = KSEA(phospho, adjMat_networKIN)    # Requires phospho data and KSEA function from "Step 6a"



############################################################
# Step 3c - Kinase Activity Prediction using Phosphosites on Kinases
############################################################
## Retrieve list of kinases
kinome = read.table("kinome_swissprot_Feb18_release.txt", header = FALSE, fill = TRUE)
# Clean up kinase data
kinome = kinome[1:3]
colnames(kinome) = c("name", "gene", "uniprotID")
kinome$uniprotID = sub("^\\(", "", kinome$uniprotID)


## Calculate mean phosphosite intensity (Kinase phosphorylation scoring or KPS) on kinases
KPS = function(phospho, kinome) {
  # phospho = a list of two dataframes with same number of columns (values = phospho data; lab = site information) 
  # kinome = a dataframe containing kinase information
  store = matrix(NA, nrow = nrow(kinome), ncol = ncol(phospho$values),
                 dimnames = list(row = kinome$name, col = colnames(phospho$values)))   #Initialize matrix
  subMatch = store   #initialize matrix (number of matched substrates for each kinase)
  
  # Character vector of UniProtID from phospho
  phosphoID = sub("_.*", "", phospho$lab$ID)
  
  for (i in 1:nrow(kinome)) {
    # Loop through each kinase in kinome data frame
    siteIndex = which(phosphoID %in% kinome$uniprotID[i])  # match phospho table to kinome table by UniProtID
    if (length(siteIndex) != 0) {
      store[i, ] = sapply(phospho$values[siteIndex, ], mean)
      subMatch[i, ] = sapply(phospho$values[siteIndex, ], length)
    }
  }
  
  return(list(score = store, subMatch = subMatch))
}
kinAct_phospho = KPS(phospho, kinome)
kinAct_phospho = lapply(kinAct_phospho, function(x) x[complete.cases(x), ])



############################################################
# Step 4 - Data Visualization
############################################################
## Visualization: pairs plot (sourced from pulsed_SILAC_ozlem.R)
plot_pairs = function(df, pattern, use_keep = FALSE) {
  # df = data frame carrying data of interest
  # pattern = regex pattern to select column of interest
  # use_keep = TRUE means to construct plot on filtered values; FALSE uses all available values
  require(gpairs)
  require(scales)
  if (use_keep) {
    plot.df = df[df$KEEP, grep(pattern, names(df), value = TRUE, perl = TRUE)]
  } else {
    plot.df = df[grep(pattern, names(df), value = TRUE, perl = TRUE)]
  }
  # AUTOMATICALLY EXCLUDE -Inf VALUES
  plot.df = as.matrix(plot.df)
  plot.df[is.infinite(plot.df)] = NA
  
  gpairs(plot.df,
         upper.pars = list(scatter = "lm"),
         scatter.pars = list(pch = 20,
                             col = alpha("black", 0.3)),
         lower.pars = list(scatter = "stats"),
         stat.pars = list(verbose = FALSE, fontsize = 15))
}
#plot_pairs(dataCF, "^LOG2.*AMO1", TRUE)
#plot_pairs(dataCF, "^LOG2.*INA6", TRUE)
#plot_pairs(dataCF, "^LOG2.*MM1S", TRUE)
#plot_pairs(dataCF, "^LOG2.*RPMI8226", TRUE)
#plot_pairs(dataCF, "^LOG2.*L363", TRUE)
#plot_pairs(dataCF, "^LOG2.*KMS11", TRUE)
#plot_pairs(dataCF, "^LOG2.*KMS34", TRUE)


## Visualization: Heatmap on sample distances in finalized phospho data (sourced form PROGENy_RNAseq)
sampleHeatmap = function(matr, col = TRUE) {
  # matr = numeric matrix for clustering
  # col = logical indicating whether samples are organized by columns
  require(dendsort)
  require(RColorBrewer)
  require(pheatmap)
  if (col) {
    sampleDists = dist(t(matr), method = "euclidean")
  } else {
    sampleDists = dist(matr, method = "euclidean")
  }
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))
  
  sampleDistMatrix = as.matrix(sampleDists)
  colors = colorRampPalette(rev(brewer.pal(9, "Blues")))(300)
  pheatmap(sampleDistMatrix,
           cluster_cols = sort_hclust(hclust(sampleDists)),
           cluster_rows = sort_hclust(hclust(sampleDists)),
           show_colnames = FALSE,
           col = colors)
}
##sampleHeatmap(as.matrix(dataCFMIA[dataCFMIA$KEEP, 56:64]))
#sampleHeatmap(as.matrix(dataCFNIA[50:57]))


## Visualization: histogram on imputed data
hist_impute = function(df, use_keep = TRUE) {
  # df = data frame containing imputed data
  # use_keep = logical indicating whether to use KEEP column to filter rows
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(stringr)
  
  if (use_keep) {
    df = filter(df, KEEP)
  }
  
  LOG2.df = dplyr::select(df, starts_with("LOG2"))
  impute.df = dplyr::select(df, starts_with("IMP"))
  
  # Reshape data into columns
  LOG2.df = gather(LOG2.df, "sample", "intensity")
  impute.df = gather(impute.df, "sample", "IMP")
  
  # Combine data
  combine.df = bind_cols(LOG2.df, impute.df["IMP"])
  
  # Create labels
  combine.df = mutate(combine.df, sample = sub("^LOG2\\.", "", sample)) %>%
    mutate(replicate = str_extract(sample, "_.*")) %>%
    mutate(replicate = sub("_", "", replicate)) %>%
    mutate(sample = sub("_.*", "", sample))
  
  ggplot(combine.df, aes(x = intensity, fill = IMP)) +
    geom_histogram(alpha = 0.3, binwidth = 0.2, position = "identity") +
    labs(x = expression("log"[2]*"-transformed Intensity"), y = "Frequency") +
    facet_grid(sample ~ replicate) +
    scale_fill_discrete(name = "Imputed",
                        breaks = c("FALSE", "TRUE"),
                        labels = c("-", "+"))
}
#hist_impute(dataCFNI, use_keep = TRUE)

## HEATMAP VISUALIZATION FUNCTION: superheat on kinase activities
superheatKIN = function(KSEA_lst) {
  # KSEA_lst = list output from KSEA function (list with three elements: score, zScore, pValue, subMatch)
  require(superheat)
  require(dendsort)
  require(RColorBrewer)
  
  superheat(X = t(KSEA_lst$score),
            
            # Set heatmap color
            #heat.pal = rev(brewer.pal(9, "RdBu")),
            heat.pal = c("lightskyblue", "black", "red"),
            heat.pal.values = c(0, 0.5, 1),
            heat.lim = c(-max(abs(KSEA_lst$score)), 
                         max(abs(KSEA_lst$score))),
            
            
            # Draw row (sample) dendrogram
            row.dendrogram = TRUE,
            
            # Hierarchial clustering
            pretty.order.cols = TRUE,
            pretty.order.rows = TRUE,
            dist.method = "euclidean",
            
            # Top bar plot data: number of sites matched to each kinase
            yt = KSEA_lst$subMatch[, 1],  # Take the first column in subMatch matrix as every column is identical
            yt.plot.type = "bar",
            yt.axis.name = "Number of\nassociated\nsubstrates",
            yt.axis.name.size = 15,
            
            # Aesthetics
            grid.hline = FALSE,
            grid.vline = FALSE,
            left.label = "variable",   # label samples
            bottom.label = "variable",   # label kinases
            yr.plot.size = 0.1,   # dendrogram size
            bottom.label.text.size = 1.7,   # default to 5
            left.label.text.size = 4,   # default to 5
            bottom.label.size = 0.15,   # size of bottom panel
            left.label.size = 0.04,   # size of left panel
            bottom.label.col = NULL,   # bottom panel color
            left.label.col = "white",   # bottom panel color
            force.bottom.label = TRUE,   # display bottom labels
            bottom.label.text.angle = 90,
            padding = 0.5,   # margins around plot (in cm)
            column.title.size = 0, 
            row.title.size = 0
  )
}


## HEATMAP VISUALIZATION FUNCTION: pheatmap on Kinase Activities
pheatmapKIN = function(KSEA_lst, fontsize_col = 5.5) {
  # KSEA_lst = list output from KSEA function (list with three elements: score, zScore, pValue, subMatch)
  # fontsize_col = numeric of length 1 indicating the size of column fonts
  require(pheatmap)
  require(dendsort)
  require(RColorBrewer)
  score = as.matrix(KSEA_lst$score)
  sampleDist = dist(t(score), method = "euclidean")   # Euclidean distances between samples
  kinaseDist = dist(score, method = "euclidean")   # Euclidean distances between kinases
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  
  #colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(300)
  colors = colorRampPalette(c("lightskyblue", "black", "red"))(300)
  breakPTS = seq(from = -max(abs(KSEA_lst$score)),
                 to = max(abs(KSEA_lst$score)),
                 length.out = 300)
  
  # Significance
  pVals = as.data.frame(KSEA_lst$pValue)
  pVals = pVals %>% 
    mutate_all(funs(case_when(. >= 0.05 ~ "",
                              . < 0.05 & . >= 0.01 ~ '*',
                              . < 0.01 & . >= 0.001 ~ '**',
                              . < 0.001 ~ '***')))
  
  pheatmap(t(score),
           cluster_cols = sort_hclust(hclust(kinaseDist)),
           cluster_rows = sort_hclust(hclust(sampleDist)),
           col = colors,
           breaks = breakPTS,
           fontsize_row = 10,
           fontsize_col = fontsize_col,
           display_numbers = t(pVals),
           number_color = "yellow", 
           fontsize_number = 15)
}


## Heatmap visualization: Kinase activities across cell lines
#superheatKIN(kinAct_networKIN)
#pheatmapKIN(kinAct_networKIN)
#superheatKIN(kinAct_pSitePlus)
#pheatmapKIN(kinAct_pSitePlus)
#superheatKIN(kinAct_phospho)
#pheatmapKIN(kinAct_phospho)


## Visualization of kinase-substrate relationships
if (FALSE) {
  # Plot historgram of # of substrates per kinase in PhosphoSitePlus
  cairo_ps("Substrates_per_kinase_pSitePlus.eps")
  hist(table(pSitePlus$GENE), breaks = 20, main = "Distribution of Number of Substrates per Kinase", 
       xlab = "Number of Substrates")
  dev.off()
  
  # Plot historgram of # of substrates per kinase from NetworKIN PREDICTIONS THAT PASSED THRESHOLD!
  cairo_ps("Substrates_per_kinase_NetworKIN.eps")
  hist(table(networKIN$Kinase.Phosphatase.Phospho.binding.domain), breaks = 30, 
       main = "Distribution of Number of Substrates per Kinase", 
       xlab = "Number of Substrates")
  dev.off()
  
  # Plot top 12 kinases with most substrates in PhosphoSitePlus
  cairo_ps("top12_Substrates_per_kinase_barplot_pSitePlus.eps")
  bar.data = table(pSitePlus$GENE)[table(pSitePlus$GENE) > 200]
  bar.data = bar.data[order(bar.data)]
  barplot(bar.data, main = "Top 12 Annotated Kinases", xlab = "", 
          ylab = "Number of Substrate", las=2)
  dev.off()
  
  # Plot top 12 kinases with most substrates in NetworKIN PREDICTIONS THAT PASSED THRESHOLD!
  cairo_ps("top12_Substrates_per_kinase_barplot_NetworKIN.eps")
  bar.data = table(networKIN$Kinase.Phosphatase.Phospho.binding.domain)[table(networKIN$Kinase.Phosphatase.Phospho.binding.domain) > 130]
  bar.data = bar.data[order(bar.data)]
  barplot(bar.data, main = "Top 12 Annotated Kinases", xlab = "", 
          ylab = "Number of Substrate", las=2)
  dev.off()
}



############################################################
# Step 5 - Map Kinases to Pathways
############################################################
## Convert kinase IDs to KEGG ID
require(httr)
# Map uniprotID to KEGG identifiers
GET(url = "http://rest.kegg.jp", path = "conv/uniprot/hsa",
    write_disk("KEGG_uniprot_ID_map.txt", overwrite = TRUE))
IDmap = read.table("KEGG_uniprot_ID_map.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(IDmap) = c("kid", "uniprotID")   # name columns
IDmap$kid = sub("^hsa:", "", IDmap$kid)   # clean characters
IDmap$uniprotID = sub("^up:", "", IDmap$uniprotID)   # clean characters


## Associate KEGG IDs (kid) to KEGG pathways
GET(url = "http://rest.kegg.jp", path = "link/hsa/pathway",
    write_disk("KEGGid_pathway_map.txt", overwrite = TRUE))   
pathMAP = read.table("KEGGid_pathway_map.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(pathMAP) = c("path", "kid")   # name columns
pathMAP$path = sub("^path:hsa", "", pathMAP$path)   # clean characters
pathMAP$kid = sub("^hsa:", "", pathMAP$kid)   # clean characters


## Extract names of KEGG human pathways
GET(url = "http://rest.kegg.jp", path = "list/pathway/hsa",
    write_disk("pathway_names.txt", overwrite = TRUE))
pathNames = read.table("pathway_names.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(pathNames) = c("path", "path_name")   # name columns
pathNames$path = sub("^path:hsa", "", pathNames$path)   # clean characters


## Filter for human signaling pathways
require(dplyr)
signalID = c("04150", "04014", "04015", "04010", "04012", "04310", "04330", "04340", 
             "04350", "04390", "04370", "04371", "04630", "04064", "04668", "04066", 
             "04068", "04020", "04070", "04072", "04071", "04024", "04022", "04151", 
             "04152", "04115", "05200")
# Merge human signaling pathways with their components (KEGG id) and respective uniprotIDs
sigPath = filter(pathNames, path %in% signalID) %>%
  left_join(pathMAP, by = "path") %>%
  left_join(IDmap, by = "kid")
rm(signalID)


## Visualize pathways using PathViewer
# Gene ID converter using gProfiler
gene_convert = function(genes, organism = "hsapiens", target = "UNIPROTSWISSPROT") {
  # gene = character vector of gene list
  # organism = character vector of length 1 indicating organism of interest
  # target = character vector of length 1 indicating target database
  require(httr)
  
  webPATH = paste0("gprofiler/gconvert.cgi?organism=",
                   organism, 
                   "&target=", 
                   target, 
                   "&output=mini&query=",
                   paste(genes, collapse = "+"))   # Construct path name
  GET(url = "https://biit.cs.ut.ee", path = webPATH,
      write_disk("gconvert.txt", overwrite = TRUE))   # Save to gconvert.txt
  read.table("gconvert.txt", sep = "\t", stringsAsFactors = FALSE)   # return gconvert.txt
}

# pSitePlus kinase name to KEGG ID conversion table (required for pathway mapping)
MAPpSitePlus = dplyr::select(pSitePlus, c(GENE, KIN_ACC_ID)) %>%
  dplyr::rename(name = GENE, uniprotID = KIN_ACC_ID)  %>% # Extract Kinase names and associated uniprotIDs
  dplyr::left_join(IDmap, by = "uniprotID")
MAPpSitePlus = unique(MAPpSitePlus)   # retain unique rows

# networKIN kinase name to KEGG ID conversion table (required for pathway mapping)
MAPnetworKIN = dplyr::select(networKIN, c(Kinase.Phosphatase.Phospho.binding.domain, 
                                          Kinase.Phosphatase.Phospho.binding.domain.STRING.ID)) %>%
  dplyr::rename(name = Kinase.Phosphatase.Phospho.binding.domain,   # Extract Kinase names and associated ENSEMBL IDs
                ENSP = Kinase.Phosphatase.Phospho.binding.domain.STRING.ID)
MAPnetworKIN = unique(MAPnetworKIN)   # retain unique rows
ENSPconvert = gene_convert(MAPnetworKIN$ENSP)    # convert ENSP to uniprotID
MAPnetworKIN = dplyr::select(ENSPconvert, V2, V4) %>%
  dplyr::rename(ENSP = V2, uniprotID = V4) %>%
  dplyr::right_join(MAPnetworKIN, by = "ENSP")
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000306043", "uniprotID"] = "P06493"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000363277", "uniprotID"] = "Q15349"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000339151", "uniprotID"] = "O14920"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000284384", "uniprotID"] = "P17252"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000374323", "uniprotID"] = "Q9UF33"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000321410", "uniprotID"] = "J3KNK1"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000352157", "uniprotID"] = "D6RJF9"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000263551", "uniprotID"] = "Q9H2X6"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000391518", "uniprotID"] = "P49761"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000383210", "uniprotID"] = "P51956"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000262211", "uniprotID"] = "Q96BR1"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000376688", "uniprotID"] = "P07948"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000351997", "uniprotID"] = "P52564"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000365012", "uniprotID"] = "P08631"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000351486", "uniprotID"] = "J3KP20"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000368124", "uniprotID"] = "Q8IU85"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000377204", "uniprotID"] = "P43250"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000313644", "uniprotID"] = "G5E948"
MAPnetworKIN[MAPnetworKIN$ENSP == "ENSP00000376684", "uniprotID"] = "O15197"
MAPnetworKIN = dplyr::left_join(MAPnetworKIN, IDmap, by = "uniprotID")


# phospho kinase name to KEGG ID conversion table (required for pathway mapping)
MAPphospho = mutate(phospho$lab, ID = sub("_.*", "", ID)) %>%
  dplyr::rename(name = Gene.name, uniprotID = ID) %>%
  dplyr::left_join(IDmap, by = "uniprotID")
MAPphospho = unique(MAPphospho)


## Map Kinases to KEGG pathway
#source("http://bioconductor.org/biocLite.R")
#biocLite("pathview")
KEGGmapper = function(mat, convertDF, pathwayID) {
  # mat = matrix containing kinAct scores from KSEA analysis
  # convertDF = data frame containing two columns "name" and "uniprotID" linking row names to uniprotID
  # pathwayID = a character vector of length 1 indicating KEGG pathway map for visualization
  require(pathview)
  
  rownames(mat) = sapply(rownames(mat),
                         function(x) convertDF[which(convertDF$name == x)[1], "kid"])
  pathview(gene.data = mat, pathway.id = pathwayID, 
           species = "hsa", 
           limit = list(gene = max(abs(mat)), cpd = 1),
           bins = list(gene = 10, cpd = 10),
           low = list(gene = "blue", cpd = NA),
           mid = list(gene = "gray", cpd = NA),
           high = list(gene = "red", cpd = NA))
}
#KEGGmapper(kinAct_networKIN$score, MAPnetworKIN, "05212")
#KEGGmapper(kinAct_pSitePlus$score, MAPpSitePlus, "05200")
#KEGGmapper(kinAct_phospho$score, MAPphospho, "05200")


# Kinase activity heatmap for KEGG pathway-specific members
KEGGheatmap = function(kinACT, convertDF, sigPath, pathwayID) {
  # kinACT = output from KSEA analysis
  # convertDF = data frame containing two columns "name" and "uniprotID" linking row names to uniprotID
  # sigPath = data frame containing KEGG pathways and associated KEGG members
  # pathwayID = a character vector of length 1 indicating KEGG pathway map for visualization
  # Requires "pheatmapKIN" function!
  
  # Convert row names of scoring matrix to KEGG IDs
  score = kinACT$score
  rowID = sapply(rownames(score),
                 function(x) convertDF[which(convertDF$name == x)[1], "kid"])
  
  # Filter for relevant pathways
  sigPath = dplyr::filter(sigPath, path == pathwayID)
  
  # Filter for rows in scoring matrix that associated with the pathway
  filter_kinACT = lapply(kinACT, function(x) x[which(rowID %in% sigPath$kid), ])   # select rows in the desired pathway
  #return(t(filter_kinACT$pValue))
  
  # Plot heatmap
  pheatmapKIN(filter_kinACT, fontsize_col = 10)
}
#KEGGheatmap(kinAct_networKIN, MAPnetworKIN, sigPath, "04630")


# Kinase activity heatmap for KEGG pathway-specific members + P-VALUE
KEGGheatmapPVAL = function(kinACT, convertDF, sigPath, pathwayID) {
  # kinACT = output from KSEA function (list containing scores and P-values)
  # convertDF = data frame containing two columns "name" and "uniprotID" linking row names to uniprotID
  # sigPath = data frame containing KEGG pathways and associated KEGG members
  # pathwayID = a character vector of length 1 indicating KEGG pathway map for visualization
  
  # Convert row names of scoring matrix to KEGG IDs
  rowID = sapply(rownames(kinACT$score),
                 function(x) convertDF[which(convertDF$name == x)[1], "kid"])
  
  # Filter for relevant pathways
  sigPath = dplyr::filter(sigPath, path == pathwayID)
  
  # Filter for rows in scoring matrix that associated with the pathway
  mat = kinACT$score[which(rowID %in% sigPath$kid), ]   # select rows in the desired pathway
  pval = kinACT$pValue[which(rowID %in% sigPath$kid), ]   # select rows in the desired pathway
  
  # Aggregate P-values using Fisher method
  # PICK P-VALUE AGGREGATION METHOD!
  require(metap)
  pval = data.frame(`P-value` = apply(pval, 1, function(x) sumlog(x)$p))
  #pval = data.frame(`P-value` = apply(pval, 1, mean))   

  # Plot heatmap
  require(pheatmap)
  require(dendsort)
  require(RColorBrewer)
  sampleDist = dist(t(mat), method = "euclidean")   # Euclidean distances between samples
  kinaseDist = dist(mat, method = "euclidean")   # Euclidean distances between kinases
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  
  #colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(300)
  colors = colorRampPalette(c("lightskyblue", "black", "red"))(300)
  breakPTS = seq(from = -max(abs(mat)),
                 to = max(abs(mat)),
                 length.out = 300)
  
  pVals = as.data.frame(kinACT$pValue[which(rowID %in% sigPath$kid), ])
  pVals = pVals %>% 
    mutate_all(funs(case_when(. >= 0.05 ~ "",
                              . < 0.05 & . >= 0.01 ~ '*',
                              . < 0.01 & . >= 0.001 ~ '**',
                              . < 0.001 ~ '***')))
  
  pheatmap(t(mat),
           cluster_cols = sort_hclust(hclust(kinaseDist)),
           cluster_rows = sort_hclust(hclust(sampleDist)),
           col = colors,
           breaks = breakPTS,
           annotation_col = pval,
           #annotation_colors = annot_colors,
           fontsize_row = 10,
           fontsize_col = 10,
           display_numbers = t(pVals),
           number_color = "yellow", 
           fontsize_number = 10)
}
KEGGheatmapPVAL(kinAct_networKIN, MAPnetworKIN, sigPath, "04014")   # Ras
KEGGheatmapPVAL(kinAct_networKIN, MAPnetworKIN, sigPath, "04010")   # MAPK
KEGGheatmapPVAL(kinAct_networKIN, MAPnetworKIN, sigPath, "04630")   # JAK-STAT
KEGGheatmapPVAL(kinAct_networKIN, MAPnetworKIN, sigPath, "04151")   # PI3K-AKT
KEGGheatmapPVAL(kinAct_networKIN, MAPnetworKIN, sigPath, "04150")   # mTOR
KEGGheatmapPVAL(kinAct_networKIN, MAPnetworKIN, sigPath, "04668")   # TNF
KEGGheatmapPVAL(kinAct_networKIN, MAPnetworKIN, sigPath, "04064")   # NF-kB


# Kinase activity heatmap for KEGG pathway-specific members + ANOVA P-VALUE
#KEGGheatmapAOV = function(kinACT, convertDF, sigPath, pathwayID) {
  # kinACT = output from KSEA function (list containing scores and P-values)
  # convertDF = data frame containing two columns "name" and "uniprotID" linking row names to uniprotID
  # sigPath = data frame containing KEGG pathways and associated KEGG members
  # pathwayID = a character vector of length 1 indicating KEGG pathway map for visualization
  # Requires "pheatmapKIN" function!
  
  # Convert row names of scoring matrix to KEGG IDs
  rowID = sapply(rownames(kinACT$score),
                 function(x) convertDF[which(convertDF$name == x)[1], "kid"])
  
  # Filter for relevant pathways
  sigPath = dplyr::filter(sigPath, path == pathwayID)
  
  # Filter for rows in scoring matrix that associated with the pathway
  mat = kinACT$score[which(rowID %in% sigPath$kid), ]   # select rows in the desired pathway
  pval = kinACT$pValue[which(rowID %in% sigPath$kid), ]   # select rows in the desired pathway
  pval = data.frame(`P-value` = pval[, 1])
  
  # Plot heatmap
  require(pheatmap)
  require(dendsort)
  require(RColorBrewer)
  sampleDist = dist(t(mat), method = "euclidean")   # Euclidean distances between samples
  kinaseDist = dist(mat, method = "euclidean")   # Euclidean distances between kinases
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  
  #colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(300)
  colors = colorRampPalette(c("lightskyblue", "black", "red"))(300)
  breakPTS = seq(from = -max(abs(mat)),
                 to = max(abs(mat)),
                 length.out = 300)
  
  pheatmap(t(mat),
           cluster_cols = sort_hclust(hclust(kinaseDist)),
           cluster_rows = sort_hclust(hclust(sampleDist)),
           col = colors,
           breaks = breakPTS,
           annotation_col = pval,
           #annotation_colors = "blue",
           fontsize_row = 10,
           fontsize_col = 10)
}
#KEGGheatmapAOV(kinAct_networKIN, MAPnetworKIN, sigPath, "04668")


# Kinase activity heatmap for all kinases + P-value
KEGGheatmapALL = function(kinACT) {
  # kinACT = output from KSEA function (list containing scores and P-values)
  # convertDF = data frame containing two columns "name" and "uniprotID" linking row names to uniprotID
  # sigPath = data frame containing KEGG pathways and associated KEGG members
  # pathwayID = a character vector of length 1 indicating KEGG pathway map for visualization
  mat = kinACT$score
  pval = kinACT$pValue
  
  # Aggregate P-values using Fisher method
  # PICK P-VALUE AGGREGATION METHOD!
  require(metap)
  pval = data.frame(`P-value` = apply(pval, 1, function(x) sumlog(x)$p))
  #pval = data.frame(`P-value` = apply(pval, 1, mean))   
  
  # Plot heatmap
  require(pheatmap)
  require(dendsort)
  require(RColorBrewer)
  sampleDist = dist(t(mat), method = "euclidean")   # Euclidean distances between samples
  kinaseDist = dist(mat, method = "euclidean")   # Euclidean distances between kinases
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  
  #colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(300)
  colors = colorRampPalette(c("lightskyblue", "black", "red"))(300)
  breakPTS = seq(from = -max(abs(mat)),
                 to = max(abs(mat)),
                 length.out = 300)
  
  pVals = as.data.frame(kinACT$pValue)
  pVals = pVals %>% 
    mutate_all(funs(case_when(. >= 0.05 ~ "",
                              . < 0.05 & . >= 0.01 ~ '*',
                              . < 0.01 & . >= 0.001 ~ '**',
                              . < 0.001 ~ '***')))
  
  pheatmap(t(mat),
           cluster_cols = sort_hclust(hclust(kinaseDist)),
           cluster_rows = sort_hclust(hclust(sampleDist)),
           col = colors,
           breaks = breakPTS,
           annotation_col = pval,
           #annotation_colors = annot_colors,
           fontsize_row = 10,
           fontsize_col = 7,
           display_numbers = t(pVals),
           number_color = "yellow", 
           fontsize_number = 5)
}
KEGGheatmapALL(kinAct_networKIN)

#save.image("workInProgress.RData")




############################################################
# Step 6 - Kinase P-value Analysis
############################################################
#load("workInProgress.RData")
## Extract expression of kinases in cell lines
RNA =  read.delim("HMCL66_CUFFLINKS_GENE_FPKM.txt", stringsAsFactors = FALSE)
# Modify sample names for conciseness
require(dplyr)
colnames(RNA)[1] = "ENSG"
colnames(RNA)[23] = "Karpas929-10"
colnames(RNA)[24] = "Karpas929-15"
colnames(RNA)[29] = "KMS11-adh"
colnames(RNA)[30] = "KMS11-sus"
colnames(RNA) = sub("_.*", "", names(RNA))
RNA = RNA[ c("ENSG", "GENE", "MM1S", "KMS11-sus", "KMS11-adh", #"U266",
             "AMO1", "INA6", "L363", "KMS34", "RPMI8226")]
# Convert FPKM to TPM
RNA[-(1:2)] = sapply(RNA[-(1:2)], function(x) x / sum(x) * 10^6)
RNA$MeanTPM = apply(RNA[-(1:2)], 1, mean)


## Convert ENSP to ENSG in "MAPnetworKIN"
ENSEMBLdf = read.csv("ENSEMBL-id-convert-180531.csv", stringsAsFactors = F)
colnames(ENSEMBLdf) = c("ENSG", "ENSP")
MAPnetworKIN = left_join(MAPnetworKIN, ENSEMBLdf, by = "ENSP")
MAPnetworKIN[MAPnetworKIN$name == "CDK1", "ENSG"] = "ENSG00000170312"
MAPnetworKIN[MAPnetworKIN$name == "RSK3", "ENSG"] = "ENSG00000071242"
MAPnetworKIN[MAPnetworKIN$name == "HIPK2", "ENSG"] = "ENSG00000064393"
MAPnetworKIN[MAPnetworKIN$name == "CLK3", "ENSG"] = "ENSG00000179335"
MAPnetworKIN[MAPnetworKIN$name == "NEK3", "ENSG"] = "ENSG00000136098"
MAPnetworKIN[MAPnetworKIN$name == "SGK3", "ENSG"] = "ENSG00000104205"
MAPnetworKIN[MAPnetworKIN$name == "Lyn", "ENSG"] = "ENSG00000254087"
MAPnetworKIN[MAPnetworKIN$name == "CaMKIdelta", "ENSG"] = "ENSG00000183049"
MAPnetworKIN[MAPnetworKIN$name == "EphB6", "ENSG"] = "ENSG00000106123"


## Compute kinase P-value using Fisher's method
kinase_pVal = apply(kinAct_networKIN$pValue, 1, function(x) sumlog(x)$p)
pValDF = data.frame(name = names(kinase_pVal),
                    FDR = kinase_pVal,
                    Phosphosites = kinAct_networKIN$subMatch[, 1])
# Annotate P-value table with ENSG ID
MAPtemp = dplyr::select(MAPnetworKIN, c("name", "ENSG"))
MAPtemp = unique(MAPtemp)
pValDF = left_join(pValDF, MAPtemp, by = "name")
# Join P-value table with RNA expression table by "ENSG" column
pValDF = left_join(pValDF, RNA, by = "ENSG")


# Clean P-value table and output
pValDF = dplyr::select(pValDF, c("ENSG", "name", "Phosphosites", "FDR", "MeanTPM", 
                                 "MM1S", "KMS11-sus", "KMS11-adh", "AMO1", "INA6", "L363", 
                                 "KMS34", "RPMI8226"))
write.csv(pValDF, "Kinase_P-values_180531.csv", quote = F, row.names = F)












### DRUG ANALYSIS
drugs = read.csv("CCLE_NP24.2009_Drug_data_2015.02.24.csv", stringsAsFactors = F)
drugs = filter(drugs, grepl("(AMO1|L363|KMS11|KMS34)", CCLE.Cell.Line.Name))


GDSC = read.csv("v17.3_fitted_dose_response.csv", stringsAsFactors = F)
GDSC = filter(GDSC, grepl("(MM1S|AMO-1|KMS-11|L-363|RPMI-8226)", CELL_LINE_NAME))
