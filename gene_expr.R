# For this file we start with gse85331_diff_exp.R but before we do the 
#  analysis looking for differentially expressed genes we will put in 
#  simulated data rather than the real data - so we can see if the 
#  analysis properly identifies which genes/rows are significant
#
# For the parts of the file that are changed, search for "simulated"


fileName <- "GSE85331_all.gene.FPKM.output.replicates.txt"


# Read in the data
data <- read.csv(fileName, sep='\t')


dim(data)     # how many rows and columns
summary(data) # summary of each column
#View(data)    # GUI viewing the data - don't do this if data is too huge


# pick out just the data (not the metadata - gene symbols, etc.)
numeric_cols <- unlist(lapply(data, is.numeric))
data_nums <- data[,numeric_cols]


# boxplot of each sample - plotted with log2 because the data had not been
#  log-scaled (which we normally do with gene expression data)
#boxplot(log2(data_nums+1), las=2)


# a factor giving the t groups for each of the columns (from data_nums)
# data frame with information about the samples
samples <- data.frame(time = factor(rep(c("d0","d0","d2","d2","d4","d4","CM","CM"),4)),
                      cell_line = factor(c(rep("H1",8), rep("H9", 8), rep("C15", 8), rep("C20", 8))))
samples




# Identify differentially expressed genes
# From the supplementary information from - https://pubmed.ncbi.nlm.nih.gov/28663367/
# "Statistical analysis was performed for each cell line individually by pairwise 
#  comparisons across time-points and day 0 (control)."
data_H1 <- data_nums[,samples$cell_line == "H1"]
#View(data_H1)


# create data_H1_simn that is simulated data 
data_H1_sim <- data_H1  # start with data_H1
data_H1_sim[,1:8] <- t(sapply(1:nrow(data_H1), function(i) {
  m <- mean(as.numeric(data_H1[i,]))
  
  # same distribution for all samples, using the median of the real data as mean
  #  - any positives in this part will be false positives
  if (i <= 10000) {     
    rnorm(8, mean=m, sd=m/10)
  }
  
  # the rest will be generated with different distributions for different samples,
  #  so any negatives in these will be false negatives
  else if (i <= 15000) {   # d0 different than d2/d4/CM 
    c(rnorm(2, mean=2*m, sd=m/10), rnorm(6, mean=m, sd=m/10))
  }
  else if (i <= 20000) { # d2 different than d0/d4/CM
    c(rnorm(2, mean=m, sd=m/10), rnorm(2, mean=2*m, sd=m/10), rnorm(4, mean=m, sd=m/10))
  }
  else if (i <= 25000) { # d0/d2 different than d4/CM
    c(rnorm(4, mean=m, sd=m/10), rnorm(4, mean=2*m, sd=m/10))
  }
  else { # increasing d0 / d2 / d4 / CM
    c(rnorm(2, mean=m, sd=m/10), rnorm(2, mean=2*m, sd=m/10), rnorm(2, mean=4*m, sd=m/10), rnorm(2, mean=8*m, sd=m/10))
  }
}))


# replace data_H1 with data_H1_sim in the rest of the file to see the analysis
#  on the simulated data


# take a look at the simulated data
boxplot(log2(data_H1_sim+1))


# We are going to do an ANOVA** for 
#  each gene (row) independently.  Let's do just the first row to get the hang of it.
row <- as.vector(t(data_H1_sim[1,]))
row
d <- data.frame(expr = row, group = samples$time[samples$cell_line == "H1"])
colnames(d) <- c("expr", "group")
results_aov <- aov(expr ~ group, d)
summary(results_aov)[[1]]$'Pr(>F)' # we just want the p-value

# perform on every row
aov_by_row <- sapply(1:(nrow(data_H1_sim)), function(i) {
  row <- as.vector(t(data_H1_sim[i,]))
  d <- data.frame(expr=row, group=samples$time[samples$cell_line == "H1"])
  colnames(d) <- c("expr", "group")
  results_aov <- aov(expr ~ group, d)
  summary(results_aov)[[1]]$'Pr(>F)' # we just want the p-value
})


# pull out p value from list
aov_by_row_p <- aov_by_row[1,]  
#View(aov_by_row_p)
hist(aov_by_row_p, breaks=20)


# false positives and false negatives
cutoff <- 0.05 # cutoff for p_value


# false positives
sum(aov_by_row_p[1:10000] <= cutoff, na.rm = TRUE)


# false negatives
sum(aov_by_row_p[10001:15000] > cutoff, na.rm = TRUE)
hist(aov_by_row_p[10001:15000], breaks=100)
sum(aov_by_row_p[15001:20000] > cutoff, na.rm = TRUE)
sum(aov_by_row_p[20001:25000] > cutoff, na.rm = TRUE)
sum(aov_by_row_p[25001:length(aov_by_row_p)] > cutoff, na.rm = TRUE)
nrow(p_values)


hist(aov_by_row_p[10001:length(aov_by_row_p)], breaks=100)



# total # declared positives
pos_count <- sum(aov_by_row_p <= cutoff, na.rm = TRUE)
457/pos_count  # false positive rate




# adjust for running multiple tests
aov_by_row_p_adj <- p.adjust(aov_by_row_p, method="fdr")


# total # declared positives
sum(aov_by_row_p_adj <= cutoff, na.rm = TRUE)


# false positives
sum(aov_by_row_p_adj[1:10000] <= cutoff, na.rm = TRUE)


# false negatives
sum(aov_by_row_p_adj[10001:15000] > cutoff, na.rm = TRUE)
sum(aov_by_row_p_adj[15001:20000] > cutoff, na.rm = TRUE)
sum(aov_by_row_p_adj[20001:25000] > cutoff, na.rm = TRUE)
sum(aov_by_row_p_adj[25001:length(aov_by_row_p)] > cutoff, na.rm = TRUE)



# put it together with gene symbol
p_values <- data.frame(p_value_adj = aov_by_row_p_adj)
rownames(p_values) <- data$X
#View(p_values)




# Can put the p values in with the data
data_H1_sim_p <- cbind(data_H1_sim, p_values, data$gene_id)
View(data_H1_sim_p)



# <= .001 (an arbitrary cutoff)
significant <- data_H1_sim_p[p_values$p_value_adj <= 0.001 & !is.na(p_values$p_value_adj),1,drop=FALSE]
View(significant)