########################################################################
### KINASE ACTIVITY RANKING BY MASS SPECTROMETRY-BASED PHOSPHO DATA
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
  cond.names = sapply(conditions, # Group column names by conditions
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
raw.data = read.delim("maxquant_output/Phospho (STY)Sites.txt", header = TRUE, stringsAsFactors = FALSE)


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
dataCFN = quantile_norm(dataCF,
                        c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34", "L363", "MM1S", "RPMI8226"),
                        use_keep = TRUE)
#dataCFN = median_centering(dataCF, TRUE)


## Imputation (remove KEEP = FALSE after imputation)
#dataCFNI = impute_data(dataCFN, use_keep = TRUE)
#dataCFNI = dplyr::filter(dataCFNI, KEEP)
dataCFNI = hybrid_impute(dataCFN, 
                         c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34", "L363", "MM1S", "RPMI8226"),
                         use_keep = TRUE)


## Average replicates
dataCFNIA = average_reps(dataCFNI, c("AMO1(?=_)", "AMO1br", "INA6", "KMS11", "KMS34",
                                     "L363", "MM1S", "RPMI8226"))


## Organize phospho data
phospho = dataCFNIA[grep("^LOG2.*(?<!_bR[0-9])$", names(dataCFNIA), value = TRUE, perl = TRUE)]

# Normalize intensities (phosphosite intensity relative to total intensity in sample)
phospho = as.data.frame(sapply(phospho, function(x) x / sum(x)))

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



############################################################
# Step 3 - Calculate K-score using PhosphoSitePlus
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


## Calculate K-score (m samples by n kinases)
KARP = function(phospho, pMatrix) {
  # phospho = a list of two dataframes with same number of columns (values = phospho data; lab = site information)
  # pMatrix = an adjacency matrix detailing relationship between kinase and substrate
  # Method described in Wilkes et al. MCP 2017
  store = matrix(NA, nrow = ncol(pMatrix), ncol = ncol(phospho$values),
                 dimnames = list(colnames(pMatrix), colnames(phospho$values)))   #Initialize matrix
  m = store   #Initialize matrix (m = number of phosphorylation sites in dataset matched to kinase)
  t = store   #Initialize matrix (t = total number of known target phosphorylation sites in PsitePlus for kinase)
  
  require(reshape2)
  for (j in 1:ncol(pMatrix)) {
    # Iterate through each kinase
    sites = rownames(pMatrix)[pMatrix[, j] == 1]   # Extract substrate ID associated with kinase
    phospho_id = phospho$values[phospho$lab$ID %in% sites, ]    # Extract normalized intensity for all associated substrates
    for (i in 1:ncol(phospho$values)) {
      # Iterate through each sample
      if (nrow(phospho_id) != 0) {
        m[j, i] = nrow(phospho_id)
        t[j, i] = length(sites)
        store[j, i] = sum(phospho_id[, i]) * (m[j, i] / t[j, i])^(1/2) * 10^6   #Calculate K-score
      }
    }
  }
  # Keep valid rows
  KEEP = !is.na(store[, 1])
  
  return(list(score = store[KEEP, ], 
              subMatch = m[KEEP, ], 
              subTotal = t[KEEP, ]))
}
score_pSitePlus  = KARP(phospho, adjMat_pSitePlus)



############################################################
# Step 4 - Visualization
############################################################
## Pairwise scatter plot for correlation of K-scores
require(scales)
pairs(score_pSitePlus$score,
      log = "xy", pch = 20, col = alpha("blue", 0.2))

## Output table
require(dplyr)
require(gProfileR)
# Organize K-scores
scores = as.data.frame(score_pSitePlus$score)
scores$meanScore = apply(score_pSitePlus$score, 1, mean)   # Aggregate K-scores across cell lines
scores$sdScore = apply(score_pSitePlus$score, 1, sd)
scores$GENE = rownames(scores)
scores$Phosphosites = score_pSitePlus$subMatch[, 1]
scores = left_join(scores, pSitePlus[c("KIN_ACC_ID", "GENE")], by = "GENE") %>%
  filter(KIN_ACC_ID != "Q15300")
scores = unique(scores)

# Convert UniProtID into ENSG
ENSGtable = gconvert(scores$KIN_ACC_ID, target = "ENSG", mthreshold = 1)
ENSGtable = ENSGtable[c("alias", "target")]
ENSGtable[] = lapply(ENSGtable, as.character)
colnames(ENSGtable) = c("KIN_ACC_ID", "ENSG")
scores = left_join(scores, ENSGtable, by = "KIN_ACC_ID")
scores[scores$GENE == "MAP4K5", "ENSG"] = "ENSG00000012983"   # manual curation

# Read RNA expression table from Keats
RNA =  read.delim("HMCL66_CUFFLINKS_GENE_FPKM.txt", stringsAsFactors = FALSE)
# Modify sample names for conciseness
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

# Join scores and RNA dataframes
scores = left_join(scores, RNA[c("ENSG", "MeanTPM")], by = "ENSG")

# Clean P-value table and output
scores = dplyr::select(scores, c("ENSG", "GENE", "Phosphosites", "meanScore", "sdScore", "MeanTPM",
                                 "AMO1", "AMO1br", "INA6", "KMS11", "KMS34", "L363", "MM1S", "RPMI8226"))
write.csv(scores, "KARP-scores-190601.csv", quote = F, row.names = F)

