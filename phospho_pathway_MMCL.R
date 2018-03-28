########################################################################
### PATHWAY ACTIVITY INFERENCE ON MASS-SPECTROMETRY BASED PHOSPHO DATA
########################################################################

############################################################
# Step 1 - Set working directory
############################################################
setwd("C:/Users/Tony Lin/Desktop/Wiita_Lab/Projects/Proteomics_project/Phosphoproteomics_MMCL/2 - MaxQuant Analysis")



############################################################
# Step 2 - Read in protein groups CSV data file
############################################################
raw.data = read.delim("maxquant_output/Phospho (STY)Sites.txt", header = TRUE, stringsAsFactors = FALSE)



############################################################
# Step 3 - Data clean-up
############################################################
## Remove reverse proteins, contaminations, and localization probability < 0.75
data = raw.data
data = data[data$Reverse != "+", ]
data = data[data$Potential.contaminant != "+", ]
data = data[as.numeric(data$Localization.prob) >= 0.75, ]
data = data[as.numeric(data$Delta.score) >= 8, ]

## Calculate log2 of intensity
Intensity.names = grep("^Intensity", names(data), value = TRUE)
LOG.names = sub("^Intensity", "LOG2", Intensity.names)
data[, LOG.names] = lapply(data[, Intensity.names], function(x) log2(x))

## Filter columns by names and extract protein IDs
filter_cols = function(df, group) {
  # df = data frame containing phosphoproteomics data
  # group = a vector of characters indicating the sample names
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
require(dplyr)
dataC = filter_cols(data)
dataC = filter(dataC, !is.na(UniProtID))   # Eliminate rows without valid identifiers
dataC = dplyr::select(dataC, -LOG2.INA6_bR2)    # INA6_bR2 correlate poorly with other bioreplicates so is removed from analysis
dataC = dplyr::select(dataC, -(LOG2.U266_bR1:LOG2.U266_bR3))   # Remove U266 from analysis


############################################################
# Step 4 - Data Filtering / Normalization
############################################################
## Group filtering to remove -Inf values (sourced from pulsed_SILAC_ozlem.R)
filter_valids = function(df, conditions, min_count, 
                         is_infinite = TRUE, at_least_one = FALSE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # is_infinite = Boolean indicating the nature of missing data
  #     if TRUE, counts Inf/-Inf values as missing values
  #     if FALSE, counts NaN values as missing values
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  require(dplyr)
  
  log2.names = names(dplyr::select(df, starts_with("LOG2")))   # Extract LOG2 columns
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  if (is_infinite) {
    cond.filter = sapply(1:length(cond.names), function(i) {
      df2 = df[cond.names[[i]]]   # Extract columns of interest
      df2 = as.matrix(df2)   # Cast as matrix for the following command
      sums = rowSums(!is.infinite(df2)) # count the number of valid values for each condition
      sums >= min_count[i]   # Calculates whether min_count requirement is met
    })
  } else {
    cond.filter = sapply(1:length(cond.names), function(i) {
      df2 = df[cond.names[[i]]]   # Extract columns of interest
      df2 = as.matrix(df2)   # Cast as matrix for the following command
      sums = rowSums(!is.nan(df2)) # count the number of valid values for each condition
      sums >= min_count[i]   # Calculates whether min_count requirement is met
    })
  }
  
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}
dataCF = filter_valids(dataC, 
                       c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34", "L363", "MM1S", "RPMI8226"),
                       min_count = rep(2, 8),
                       is_infinite = TRUE,
                       at_least_one = TRUE)


# Visualization: pairs plot (sourced from pulsed_SILAC_ozlem.R)
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


## Quantile normalization for each cell line (alternative to global median normalization)
quantile_norm = function(df, conditions, use_keep = TRUE) {
  # df = data frame containing LOG2 columns for normalization
  # conditions = a character vector dictating the grouping
  # use_keep = filter rows using KEEP column prior to normalization
  require(dplyr)
  require(limma)
  
  log2.names = names(dplyr::select(df, starts_with("LOG2")))   # Extract LOG2 columns
  cond.names = sapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  if (use_keep) df = df[df$KEEP, ]
  LOG2_df = as.matrix(df[log2.names])
  LOG2_df[is.infinite(LOG2_df)] = NA
  
  for(group in cond.names) LOG2_df[, group] = normalizeQuantiles(LOG2_df[, group])
  
  df[log2.names] = LOG2_df[, log2.names]
  
  return(df)
}
dataCFQ = quantile_norm(dataCF,
                        c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34", "L363", "MM1S", "RPMI8226"))


## Global median normalization for each sample (center median at 0)
median_centering = function(df, use_keep = TRUE) {
  # df = data frame containing LOG2 columns for normalization
  # use_keep = filter rows using KEEP column prior to computing the median for normalization
  LOG2.names = grep("^LOG2", names(df), value = TRUE)
  if (use_keep) {
    df[, LOG2.names] = lapply(LOG2.names, 
                              function(x) {
                                LOG2 = df[df$KEEP, x]
                                LOG2[is.infinite(LOG2)] = NA   # Exclude 0 intensity values from median calculation
                                gMedian = median(LOG2, na.rm = TRUE)
                                df[, x] - gMedian
                                }
                              )
  } else {
    df[, LOG2.names] = lapply(LOG2.names, 
                              function(x) {
                                LOG2 = df[, x]
                                LOG2[is.infinite(LOG2)] = NA   # Exclude 0 intensity values from median calculation
                                gMedian = median(LOG2, na.rm = TRUE)
                                df[, x] - gMedian
                              }
    )
  }
  return(df)
}
#dataCFM = median_centering(dataCF, TRUE)



############################################################
# Step 5 - Data Imputation
############################################################
## Impute missing data METHOD 1: Downshift and tighten the sampling distribution of missing values
impute_data = function(df, width = 0.3, downshift = 1.8, use_keep = TRUE) {
  # df = data frame containing filtered 
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  # use_keep = filter rows using KEEP column prior to imputation
  
  LOG2.names = grep("^LOG2", names(df), value = TRUE)
  impute.names = sub("^LOG2", "impute", LOG2.names)
  
  # Create new columns indicating whether the values are imputed
  df[impute.names] = lapply(LOG2.names, function(x) is.infinite(df[, x]))
  
  # Imputation
  set.seed(1)
  df[LOG2.names] = lapply(LOG2.names, 
                          function(x) {
                            temp = df[, x]
                            temp[is.infinite(temp)] = NA
                            
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
#dataCFMI = impute_data(dataCFM, use_keep = TRUE)


## Impute missing data METHOD 2: Hybrid missing data imputation
hybrid_impute = function(df, conditions, width = 0.3, downshift = 1.8, use_keep = TRUE) {
  # df = data frame containing filtered 
  # conditions = a character vector dictating the grouping
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  # use_keep = filter rows using KEEP column prior to imputation
  require(dplyr)
  
  log2.names = names(dplyr::select(df, starts_with("LOG2")))   # Extract LOG2 columns
  impute.names = sub("^LOG2", "IMP", log2.names)
  cond.names = sapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE, perl = TRUE))
  
  if (use_keep) df = df[df$KEEP, ]
  LOG2_df = as.matrix(df[log2.names])
  LOG2_df[is.infinite(LOG2_df)] = NA
  
  # Create new columns indicating whether the values are imputed
  df[impute.names] = lapply(log2.names, function(x) is.na(df[, x]))p
  
}
dataCFQH = hybrid_impute(dataCFQ, use_keep = TRUE)
  


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
  impute.df = dplyr::select(df, starts_with("impute"))
  
  # Reshape data into columns
  LOG2.df = gather(LOG2.df, "sample", "intensity")
  impute.df = gather(impute.df, "sample", "impute")
  
  # Combine data
  combine.df = bind_cols(LOG2.df, impute.df["impute"])
  
  # Create labels
  combine.df = mutate(combine.df, sample = sub("^LOG2\\.", "", sample)) %>%
    mutate(replicate = str_extract(sample, "_.*")) %>%
    mutate(replicate = sub("_", "", replicate)) %>%
    mutate(sample = sub("_.*", "", sample))

  ggplot(combine.df, aes(x = intensity, fill = impute)) +
    geom_histogram(alpha = 0.3, binwidth = 0.2, position = "identity") +
    labs(x = expression("log"[2]*"-transformed Intensity"), y = "Frequency") +
    facet_grid(sample ~ replicate) +
    scale_fill_discrete(name = "Imputed",
                        breaks = c("FALSE", "TRUE"),
                        labels = c("-", "+"))
  
}
hist_impute(dataCFMI, use_keep = TRUE)


# Average biological replicates
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
dataCFMIA = average_reps(dataCFMI, c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34",
                                     "L363", "MM1S", "RPMI8226", "U266"))


# Visualization: Heatmap on finalized phospho data (sourced form PROGENy_RNAseq)
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
sampleHeatmap(as.matrix(dataCFMIA[dataCFMIA$KEEP, 56:64]))


# Visualization: Pairs plot on final sample data set (SLOW!)
#plot_pairs(dataCFMIA, "^LOG2.*(?<!_bR[0-9])$", use_keep = TRUE)

# Save workspace
#save.image("checkpoint_phospho.RData")


############################################################
# Step 6a - Kinase Set Enrichment Analysis + PhosphoSitePlus
############################################################
#load("checkpoint_phospho.RData")
## Organize phospho data
phospho = dataCFMIA[dataCFMIA$KEEP, grep("^LOG2.*(?<!_bR[0-9])$", names(dataCFMIA),
                                         value = TRUE, perl = TRUE)]  # Transfer filtered data
# Normalize intensities such that the mean intensity is 0 for each site (relative changes in phosphosite intensity)
phospho = as.data.frame(t(apply(as.matrix(phospho), 1, function(x) x - mean(x))))
phospho$ID = paste0(dataCFMIA[dataCFMIA$KEEP, "UniProtID"], "_",
                    dataCFMIA[dataCFMIA$KEEP, "Amino.acid"],
                    dataCFMIA[dataCFMIA$KEEP, "Position"]) # Organize phospho matrix ([UniProtID]_[Amino Acid][Position])
phospho$ID = toupper(phospho$ID)   # Raise IDs to upper case
phospho$Gene.name = dataCFMIA[dataCFMIA$KEEP, "Gene.Name"]
# Organize into list with 2 elements, value (phospho intensities) and lab (phospho site/annotation)
require(dplyr)
phospho = list(values = dplyr::select(phospho, starts_with("LOG2")),
               lab = dplyr::select(phospho, ID:Gene.name))
colnames(phospho$values) = sub("^LOG2.", "", colnames(phospho$values))


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
  
  mP = mean(as.matrix(phospho$values))
  sdP = sd(as.matrix(phospho$values))
  
  for (j in 1:ncol(pMatrix)) {
    # Iterate through each kinase
    sites = rownames(pMatrix)[pMatrix[, j] == 1]   # Extract substrate ID associated with kinase
    phospho_id = phospho$values[phospho$lab$ID %in% sites, ]    # Extract log2FCs for all associated substrates
    for (i in 1:ncol(phospho$values)) {
      # Iterate through each sample
      if (nrow(phospho_id) != 0) {
        store[j, i] = mean(phospho_id[, i])   #Calculate KSEA score
        zScore[j, i] = (store[j, i] - mP * sqrt(nrow(phospho_id))) / sdP   #Calculate zscore
        subMatch[j, i] = nrow(phospho_id)   #Save number of matched substrates
      }
    }
  }
  pValue = 2*pnorm(-abs(zScore))   # 2-sided test
  
  return(list(score = store, zScore = zScore, pValue = pValue, subMatch = subMatch))
}
kinAct_pSitePlus  = KSEA(phospho, adjMat_pSitePlus)
kinAct_pSitePlus = lapply(kinAct_pSitePlus, function(x) x[complete.cases(x), ])   #Remove rows with NA
###ONLY 700 QUANTIFIED SITES ARE ANNOTATED IN PhosphoSitePlus, ASSOCIATED WITH 173 KINASES (meaning we're tossing a lot of phospho data)
###IN TOTAL, 1219 KINASE-SUBSTRATE RELATIONSHIPS ARE ESTABLISHED USING PhosphoSitePlus




############################################################
# Step 6b - Kinase Set Enrichment Analysis + NetworKIN
############################################################
## Produce a table of sites for NetworKIN to predict possible kinase associations
#write.table(dataCFMIA[dataCFMIA$KEEP, c("UniProtID", "Position", "Amino.acid")],
#            file = "all-filtered-sites.tsv", 
#            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


## Load NetworKIN phosphosite-kinase predictions
networKIN = read.table(file = 'NetworKIN_predictions_180216.tsv', sep = '\t', header = TRUE,
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
kinAct_networKIN = lapply(kinAct_networKIN, function(x) x[complete.cases(x), ])   #Remove rows with NA



############################################################
# Step 6c - Kinase Activity Prediction using Phosphosites on Kinases
############################################################
## Retrieve list of kinases
kinome = read.table("kinome_swissprot_Feb18_release.txt", header = FALSE, fill = TRUE)
# Clean up kinase data
kinome = kinome[1:3]
colnames(kinome) = c("name", "gene", "uniprotID")
kinome$uniprotID = sub("^\\(", "", kinome$uniprotID)


## Calculate mean phosphosite intensity (Kinase phosphorylation scoring or KAS) on kinases
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
# Step 7 - Kinase Activity Heatmap
############################################################
## HEATMAP VISUALIZATION FUNCTION: superheat on kinase activities
superheatKIN = function(KSEA_lst) {
  # KSEA_lst = list output from KSEA function (list with three elements: score, zScore, pValue, subMatch)
  require(superheat)
  require(dendsort)
  require(RColorBrewer)
  
  superheat(X = t(KSEA_lst$score),
            
            # Set heatmap color
            heat.pal = rev(brewer.pal(9, "RdBu")),
            
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
pheatmapKIN = function(KSEA_lst) {
  # KSEA_lst = list output from KSEA function (list with three elements: score, zScore, pValue, subMatch)
  require(pheatmap)
  require(dendsort)
  require(RColorBrewer)
  score = as.matrix(KSEA_lst$score)
  sampleDist = dist(t(score), method = "euclidean")   # Euclidean distances between samples
  kinaseDist = dist(score, method = "euclidean")   # Euclidean distances between kinases
  
  sort_hclust <- function(...) as.hclust(dendsort(as.dendrogram(...)))  # For ordering dendrogram
  
  colors = colorRampPalette(rev(brewer.pal(11, "RdBu")))(300)
  pheatmap(t(score),
           cluster_cols = sort_hclust(hclust(kinaseDist)),
           cluster_rows = sort_hclust(hclust(sampleDist)),
           col = colors,
           fontsize_row = 10,
           fontsize_col = 5.5)
}


## Heatmap visualization: Kinase activities across cell lines
superheatKIN(kinAct_networKIN)
pheatmapKIN(kinAct_networKIN)
superheatKIN(kinAct_pSitePlus)
pheatmapKIN(kinAct_pSitePlus)
superheatKIN(kinAct_phospho)
pheatmapKIN(kinAct_phospho)


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

# Save workspace
#save.image("checkpoint2_phosphO.RData")



############################################################
# Step 8 - Map Kinases to Pathways
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
             "04152")
sigPath = filter(pathNames, path %in% signalID)
# Merge human signaling pathways with their components (KEGG id) and respective uniprotIDs
sigPath = left_join(sigPath, pathMAP, by = "path")
sigPath = left_join(sigPath, IDmap, by = "kid")
rm(signalID)


## Filter for pathways in human cancer (path:hsa05200)
cancerPath = filter(pathNames, path %in% "04115")
# Merge human cancer pathways with their components (KEGG id) and respective uniprotIDs
cancerPath = left_join(cancerPath, pathMAP, by = "path")
cancerPath = left_join(cancerPath, IDmap, by = "kid")
rm(pathNames, pathMAP, IDmap)


## Visualize pathways using PathViewer
#source("http://bioconductor.org/biocLite.R")
#biocLite("pathview")
require(pathview)
x = kinAct3$score
rownames(x) = sapply(rownames(x), function(y) kinome[which(kinome$name == y)[1], "uniprotID"])
rownames(x) = sapply(rownames(x), function(y) IDmap[which(IDmap$uniprotID == y)[1], "kid"])
pathview(x2, pathway.id = "04115", species = "hsa", limit = list(gene = max(abs(x2)), cpd = 1))

y = kinAct$score
rownames(y) = sapply(rownames(y), function(z) pSitePlus[which(pSitePlus$GENE == z)[1], "KIN_ACC_ID"])
rownames(y) = sapply(rownames(y), function(z) IDmap[which(IDmap$uniprotID == z)[1], "kid"])
pathview(y2, pathway.id = "04115", species = "hsa", limit = list(gene = max(abs(y2)), cpd = 1))

z = kinAct2$score
ENSEMBLid = data.frame(kinase = networKIN$Kinase.Phosphatase.Phospho.binding.domain,
                       ID = networKIN$Kinase.Phosphatase.Phospho.binding.domain.STRING.ID,
                       stringsAsFactors = FALSE)
ENSEMBLid$uniprotID2 = test$V4
ENSEMBLid$uniprotID1[ENSEMBLid$kinase == "RSK3"] = "Q15349"
ENSEMBLid$uniprotID1[ENSEMBLid$kinase == "IKKbeta"] = "O14920"
ENSEMBLid$uniprotID1[ENSEMBLid$kinase == "PKCalpha"] = "P17252"
rownames(z) = sapply(rownames(z), function(i) ENSEMBLid[which(ENSEMBLid$kinase == i)[1], "uniprotID1"])
rownames(z) = sapply(rownames(z), function(i) IDmap[which(IDmap$uniprotID == i)[1], "kid"])
pathview(z2, pathway.id = "04115", species = "hsa", limit = list(gene = max(abs(z2)), cpd = 1))



test = read.table("test.txt", sep = "\t", stringsAsFactors = F)

GET(url = "https://biit.cs.ut.ee", path = "gprofiler/gconvert.cgi?organism=hsapiens&target=AFFY_HG_U133_PLUS_2&output=mini&query=POU5F1+SOX2+NANOG&prefix=AFFY_HUGENE_1_0_ST_V1",
    write_disk("gconvert.txt", overwrite = TRUE))



save.image("workInProgress.RData")
