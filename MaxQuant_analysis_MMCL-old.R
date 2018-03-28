### THIS SUITE OF CODE CLEANS AND REDUCES MAXQUANT OUTPUT ON PHOSPHOPROTEOMICS DATA

# Step 1 - Set working directory and increase memory size for xlsx file writing
options(java.parameters = "- Xmx8000m")
setwd("C:/Users/Tony Lin/Desktop/Wiita_Lab/Projects/Proteomics_project/Phosphoproteomics_MMCL/2 - MaxQuant Analysis")



# Step 2 - Read in protein groups CSV data file
raw.data = read.csv("Phospho (STY)Sites.csv", header = TRUE, colClasses = "character")



# Step 3 - Data clean-up
## Remove reverse proteins, contaminations, and localization probability < 0.75
data = raw.data
data = data[data$Reverse != "+", ]
data = data[data$Potential.contaminant != "+", ]
data = data[as.numeric(data$Localization.prob) >= 0.75, ]



# Step 4 - Data-specific manipulations
## Cast peptide intensity columns as numeric type
intensity.names = grep("^Intensity", names(data), value = TRUE)
data[intensity.names] = lapply(data[intensity.names], function(x) as.numeric(x))



## Rearrange the phosphopeptide intensities into a single column for each sample
collapse_phospho = function(df) {
  # Collapses intensity columns for singly, doubly and triply phosphorylated peptides 
  # into a single column for each sample
  
  # Extract relevant column names
  sample.names = grep("^Intensity\\.", names(df), value = TRUE)
  mult_index = grep("___[[:digit:]]$", sample.names)
  
  # Applies operation on Intensity columns
  require(tidyr)
  df2 = lapply(sample.names[-mult_index], 
               function(x) gather(df, Multiplicity, ALLintensity,
                                  grep(x, sample.names[mult_index], value = TRUE)))
  
  # Rearrage data for output
  df3 = df2[[1]]
  df3 = df3[-grep("^Intensity.*___[[:digit:]]$", names(df3))]   # Removes *___1,2,3* columns
  df3 = df3[-c(ncol(df3))]   # df3 is the base dataframe without ALLintensity columns
  df3$Multiplicity = as.numeric(gsub(".*___", "", df3$Multiplicity))
  ALLintensity = lapply(df2, function(x) as.numeric(x[, ncol(x)]))
  ALLintensity = do.call(cbind.data.frame, ALLintensity)
  colnames(ALLintensity) = paste0("ALL", sample.names[-mult_index])
  
  return(cbind(df3, ALLintensity))
}
data = collapse_phospho(data)


## Calculate log2(ALLintensity)
ALLintensity.names = grep("^ALLIntensity\\.", names(data), value = TRUE)
LOG.names = sub("^ALLIntensity", "LOG2", ALLintensity.names)
data[, LOG.names] = lapply(data[, ALLintensity.names], function(x) log2(x))


## Filter columns by names
filter_cols = function(df, group) {
  
  
  # Extract protein name (WORKS WITH SWISS-PROT DATABASE)
  ProteinName = strsplit(df$Fasta.headers, ";")
  FastaHeader = unlist(lapply(ProteinName, function(x) x[1]))
  regex = regexpr("(?<=_HUMAN ).*(?= OS=.*)", FastaHeader, perl = TRUE)
  out = rep(NA, length(FastaHeader))   # This fills in NA for rows without matching patterns
  out[regex != -1] = regmatches(FastaHeader, regex)
  df$Protein.Name = out
  
  # Extract gene name (WORKS WITH SWISS-PROT DATABASE)
  regex2 = regexpr("(?<=\\|[[:alnum:]]{6}\\|).*(?=_HUMAN)", FastaHeader, perl = TRUE)
  out2 = rep(NA, length(FastaHeader))
  out2[regex2 != -1] = regmatches(FastaHeader, regex2)
  df$Gene.Name = out2
  
  # Extract UniProt ID
  regex3 = regexpr("(?<=\\|).*(?=\\|)", FastaHeader, perl = TRUE)
  out3 = rep(NA, length(FastaHeader))
  out3[regex3 != -1] = regmatches(FastaHeader, regex3)
  df$UniProtID = out3
  
  group.names = grep(group, names(df), value = TRUE)
  ALLintensity.names = grep("^ALLIntensity", group.names, value = TRUE)
  IDtype.names = grep("^Identification.type", group.names, value = TRUE)
  LOG2.names = grep("^LOG2", group.names, value = TRUE)
  df$PhosphoSite = paste0(df$Amino.acid, df$Positions.within.proteins)
  
  keepCols = c("UniProtID", "Gene.Name", "Protein.Name", "Phospho..STY..Probabilities", 
               "Position.in.peptide", "PhosphoSite", "Multiplicity", 
               LOG2.names, ALLintensity.names, IDtype.names, "Sequence.window", 
               "Peptide.window.coverage", "Localization.prob")
  return(df[, keepCols])
}
HS5.data = filter_cols(data, "HS5")
U266.data = filter_cols(data, "U266")


## Group filtering to remove -Inf values
NA_filter = function(df, groupList = list("DMSO", "Tof"), 
                     min_number = 1, impute = TRUE) {
  # df = data frame containing log2 intensity data
  # groupList = a list of vectors containing replicate groups
  # min_number = filter out rows that do not meet the minimum number of valid values 
  #       in each group
  
  #Count up all valid values in each group
  validSum = lapply(groupList,
                    function(group) {
                      group.names = grep(paste0("LOG2.*", group), 
                                         colnames(df), value = TRUE)
                      return(apply(data.frame(df[, group.names]), 1, 
                                          function(x) sum(x != -Inf)))
                    })
  #Mark rows that passes the min_number filter
  filter = apply(data.frame(validSum), 1, function(x) all(x >= min_number))
  
  return(df[filter, ])
}
HS5.data = NA_filter(HS5.data)
U266.data = NA_filter(U266.data)


## Global median normalization for each sample
gMedian_norm = function(df) {
  # Performs median normalization on LOG2 columns (subtract column median-- computed
  # based on valid values-- from each value)
  LOG2.names = grep("^LOG2", names(df), value = TRUE)
  df[, LOG2.names] = lapply(LOG2.names, 
                            function(x) {
                              LOG2 = df[, x]
                              LOG2[is.infinite(LOG2)] = NA
                              gMedian = median(LOG2, na.rm = TRUE)
                              return(df[, x] - gMedian)})
  return(df)
}
HS5.data = gMedian_norm(HS5.data)
U266.data = gMedian_norm(U266.data)


## Impute missing data
impute_data = function(df, width = 0.3, downshift = 1.8) {
  # Assumes missing data (in df) follows a narrowed and downshifted normal distribution
  LOG2.names = grep("^LOG2", names(df), value = TRUE)
  impute.names = sub("^LOG2", "impute", LOG2.names)
  
  # Create new columns indicating whether the values are imputed
  df[, impute.names] = lapply(LOG2.names, function(x) is.infinite(df[, x]))
  
  # Imputation
  set.seed(1)
  df[, LOG2.names] = lapply(LOG2.names, 
                            function(x) {
                              temp = df[, x]
                              temp[is.infinite(temp)] = NA
                              temp.sd = width * sd(temp, na.rm = TRUE)
                              temp.mean = mean(temp, na.rm = TRUE) - downshift * sd(temp, na.rm = TRUE)
                              n.missing = sum(is.na(temp))
                              temp[is.na(temp)] = rnorm(n.missing, mean = temp.mean, sd = temp.sd)                          
                              return(temp)
                            })
  return(df)
}
HS5.data = impute_data(HS5.data)
U266.data = impute_data(U266.data)


## Calculate LOG2 difference (ONLY WORKS WHEN LENGTH OF groupList EQUALS 2)
LOG2_diff = function(df, groupList = list("Tof", "DMSO"), identifier = "") {
  LOG2.names = grep("^LOG2", colnames(df), value = TRUE)
  mean.names = paste0("mean.LOG2.", identifier, ".", groupList)
  df[, mean.names] = lapply(groupList, 
                            function(x) {
                              group.names = grep(x, LOG2.names, value = TRUE)
                              return(rowMeans(data.frame(df[, group.names])))
                            })
  df[, paste0("LOG2.FC.", identifier)] = df[, mean.names[1]] - df[, mean.names[2]]
  return(df)
}
HS5.data = LOG2_diff(HS5.data, identifier = "HS5")
U266.data = LOG2_diff(U266.data, identifier = "U266")


# Step 5 - Output Excel file with processed data
if(FALSE) {
require(xlsx)
wb <- createWorkbook()

addDataFrame(raw.data, row.names = FALSE, 
             createSheet(wb, sheetName="raw_data"))
addDataFrame(data, row.names = FALSE, 
             createSheet(wb, sheetName="cleaned_data"))
addDataFrame(HS5.data, row.names = FALSE, 
             createSheet(wb, sheetName="HS5_data"))
addDataFrame(U266.data, row.names = FALSE, 
             createSheet(wb, sheetName="U266_data"))

saveWorkbook(wb, "HS5_U266_monoculture_phospho_processed.xlsx")
}