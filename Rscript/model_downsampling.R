# AIM: To identify downsampling total reads

# import libraries

library(ggplot2)
library(reshape2)
options("scipen"=100, "digits"=4)
# read initQC
init.qc <- read.csv("~/Documents/GSI/Downsampling/initQC.txt", sep = "\t", as.is = T)
head(init.qc)
total.number.of.samples <- dim(init.qc)[1]
# plot total reads Vs PF reads

p <- ggplot(init.qc) + 
  geom_point(aes(x = TOTAL_READS, y = PF_UNIQUE_READS + MEAN_TARGET_COVERAGE))
p

p <- ggplot(init.qc) + 
  geom_point(aes(x = TOTAL_READS, y = PF_UQ_READS_ALIGNED))
p

p <- ggplot(init.qc) + 
  geom_point(aes(x = TOTAL_READS, y = MEAN_TARGET_COVERAGE))
p


heatmap(z)
# z <- data.frame(z)

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

total.reads.mean <- rowSums(z)/length(downsample_ratio)
total.reads.var <- RowVar(z)
summ <- data.frame(summary(t(z)))


zplot <- melt(z)
zplot$log10Var1 <- round(log10(zplot$Var1),3)
zpl <- ggplot(zplot) + 
  geom_point(aes(x = reorder(log10Var1, -Var1), y = value / (10 ^ 6), color = Var2)) + 
  xlab("log10 downsampling fraction") + ylab("(Total reads x 10^6) after downsampling") + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 90)) 
  # geom_line(data = melt(z), aes(x = Var1, y = value, color = Var2, group = Var2))

zpl



zpl <- ggplot(zplot, aes(x = reorder(Var1, -Var1), y = value / (10 ^ 6), group = Var2, color = Var2)) + 
  geom_point() + 
  xlab("Downsampling fraction") + ylab("(Total reads x 10^6) after downsampling") + 
  theme(legend.position = "none") 
  
# geom_line(data = melt(z), aes(x = Var1, y = value, color = Var2, group = Var2))

zpl


# fitting function here 

mod <- nls(I~exp(a+b*s+c*s^2), start=list(a=0, b=0, c=0))
summary(mod)


zpl <- ggplot(zplot, aes(x = reorder(Var1, -Var1), y = value / (10 ^ 6), color = Var1, group = Var1)) + 
  # geom_point() +
  geom_boxplot() +
  # geom_point() +
  xlab("Downsampling fraction") + ylab("(Total reads x 10^6) after downsampling") + 
  theme(legend.position = "none") 

# geom_line(data = melt(z), aes(x = Var1, y = value, color = Var2, group = Var2))

zpl



# read mean target coverage 

QC_dir <- "/Users/prath/Documents/GSI/Downsampling"

QC.files <- list.files(QC_dir, full.names = TRUE, pattern = ".txt")

QC <- rep(NA, 60)
for (fl in QC.files){
  qc <- read.csv(fl, sep = "\t", as.is = T)
  QC <- rbind(QC, qc)
}
QC <- data.frame(QC[-1,])

QC$downsamplingRatio <- as.numeric(paste(stringr::str_split_fixed(QC$sample_name, "_", 11)[,10],
                              stringr::str_split_fixed(QC$sample_name, "_", 11)[,11], 
                              sep = "."))

QC$sample_name <- apply(stringr::str_split_fixed(QC$sample_name, "_", 11)[,1:7], 1, paste, collapse="_")
head(QC)

# 

#  calculate the expected number of total reads, PF unique reads and PF uniq reads aligned from downsampling ratios
expectedDSTotalReads <- c()
total.reads <- QC[QC$downsamplingRatio == 1.0,]

downsample_ratio <- unique(QC$downsamplingRatio) # read csv from file
z <- matrix(NA, nrow=length(downsample_ratio), ncol = length(total.reads$sample_name))
row.names(z) <- as.character(downsample_ratio)
colnames(z) <- total.reads$sample_name
for (ds in downsample_ratio){
  sample <- total.reads$sample_name
  z[as.character(ds), sample] <- total.reads$TOTAL_READS * ds
}

zpl <- melt(z)
colnames(zpl) <- c("downsamplingRatio",
                   "sample_name", "expected_total_reads")
zpl <- data.frame(zpl)
head(zpl)


QC <- merge(QC, zpl, by = c("downsamplingRatio", "sample_name"), all.x= TRUE)

# plot Total reads vs doansampling ratio

plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = expected_total_reads, size = downsamplingRatio), color = "red", alpha = 0.5) +
  # geom_line(data = QC, aes(x = 1-downsamplingRatio, y = expected_total_reads, group = sample_name), color = "red", alpha = 0.5) +
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = TOTAL_READS), color = "blue", alpha = 0.25) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = TOTAL_READS, group = sample_name), color = "blue", alpha = 0.25) 
plt

# plot total reads vs mean target coverage
plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(log10(TOTAL_READS), -TOTAL_READS), y = MEAN_TARGET_COVERAGE, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_line(data = QC, aes(x = reorder(log10(TOTAL_READS), -TOTAL_READS), y = MEAN_TARGET_COVERAGE, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("log(10) Total reads") + scale_y_log10() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) 
  
plt


plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(TOTAL_READS, -TOTAL_READS), y = MEAN_TARGET_COVERAGE, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_line(data = QC, aes(x = reorder(TOTAL_READS, -TOTAL_READS), y = MEAN_TARGET_COVERAGE, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Total reads") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) 

plt


plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(TOTAL_READS, -TOTAL_READS), y = MEDIAN_TARGET_COVERAGE, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_line(data = QC, aes(x = reorder(TOTAL_READS, -TOTAL_READS), y = MEDIAN_TARGET_COVERAGE, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Total reads") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) 

plt


plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(TOTAL_READS, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_line(data = QC, aes(x = reorder(TOTAL_READS, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Total reads") + ylab("% mapping") +  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) 

plt

plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(TOTAL_READS, -TOTAL_READS), y = PCT_PF_UQ_READS_ALIGNED, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_line(data = QC, aes(x = reorder(TOTAL_READS, -TOTAL_READS), y = PCT_PF_UQ_READS_ALIGNED, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Total reads") +  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2)) 

plt

# downsampling ratio vs rest

plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_PF_UQ_READS_ALIGNED, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_PF_UQ_READS_ALIGNED, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Downsamping ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

plt



plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = MEAN_TARGET_COVERAGE, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = MEAN_TARGET_COVERAGE, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Downsamping ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

plt

# report the following plot (plot1)
# plt <- ggplot() +
#   geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, size = downsamplingRatio), color = "red", alpha = 0.5) +
#   geom_line(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, group = sample_name), color = "blue", alpha = 0.5) + 
#   xlab("Downsampling ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("% mapping")
# 
# plt


plt <- ggplot() +
  # geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_boxplot(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, group = downsamplingRatio), fill = "grey", alpha = 0.5) + 
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, size = downsamplingRatio), color = "red", alpha = 0.25) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Downsampling ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("% bases with MAPQ > 20")

plt


plt <- ggplot() +
  # geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_boxplot(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_USABLE_BASES_ON_TARGET, group = downsamplingRatio), fill = "grey", alpha = 0.5) + 
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_USABLE_BASES_ON_TARGET, size = downsamplingRatio), color = "red", alpha = 0.25) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_USABLE_BASES_ON_TARGET, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Downsampling ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
plt


# for each downsamping ratio compute sd
qc.sd <- c()
for (ds in downsample_ratio){
  qc.sd[[as.character(ds)]] <- sd(QC[QC$downsamplingRatio == ds,]$PCT_USABLE_BASES_ON_TARGET)
}
QC$PCT_USABLE_BASES_ON_TARGET_sd <- qc.sd[as.character(QC$downsamplingRatio)]

estimated_downsampling_cuttoff <- 1

v.init <- qc.sd["1"]
v.cutoff <- qc.sd - v.init
v.cutoff <- data.frame(v.cutoff)
v.cutoff$downsamplingRatio <- as.numeric(row.names(v.cutoff))

# finding the most optimal cut off 
v.cutoff$opt <- v.cutoff$v.cutoff * (10 ^ (-1 * log10(downsample_ratio)))
v.cutoff$opt_round <- round(v.cutoff$opt)
v.cutoff$decision <- ifelse(v.cutoff$opt_round > 0, "cut", "opt")

# OUTPUT
ds.cutoff.intercept <- max(v.cutoff[v.cutoff$decision == "cut",]$downsamplingRatio)
sd.cutoff.intercept <- v.cutoff[v.cutoff$decision == "cut" &
                                  v.cutoff$downsamplingRatio == ds.cutoff.intercept ,]$v.cutoff


ggplot(v.cutoff, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = v.cutoff))+ 
  geom_point(aes(color = decision), size=3) +
  geom_line(linetype="dashed", group = 1) + xlab("Downsamping ratio") + ylab("Standard deviation") +
  geom_hline(yintercept = sd.cutoff.intercept, linetype = "dotted", color = "red")



plt <- ggplot(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_USABLE_BASES_ON_TARGET)) +
  geom_point(aes(group = downsamplingRatio), color = "black", size = 1) +
  geom_line(aes(group = sample_name), color = "blue", alpha = 0.5)+
  geom_boxplot(aes(group = downsamplingRatio), fill = "pink", alpha = 0.5)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),
               geom="errorbar", color="red", width=0.5) +
  stat_summary(fun.y=mean, geom="point", color="red", ) +
  theme_classic() + 
  xlab("Downsampling ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "none") 
 
plt


ds.cutoff.intercept

## now go back to the reads vs downsampling plot 

plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = expected_total_reads, size = downsamplingRatio), color = "red", alpha = 0.5) +
  # geom_line(data = QC, aes(x = 1-downsamplingRatio, y = expected_total_reads, group = sample_name), color = "red", alpha = 0.5) +
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = TOTAL_READS), color = "blue", alpha = 0.25) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = TOTAL_READS, group = sample_name), color = "blue", alpha = 0.25) 
plt

# obtain / replot for downsampling cut off

cut_off_subset <- QC[QC$downsamplingRatio == ds.cutoff.intercept,]

plt <- ggplot(data = cut_off_subset, aes(x = TOTAL_READS, y = 1 - PCT_EXC_MAPQ)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_point(data = cut_off_subset, 
             aes(x = TOTAL_READS, y = PCT_USABLE_BASES_ON_TARGET), color = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(paste0("Total downsampled reads at downsamping ratio ", ds.cutoff.intercept)) + 
  ylab("% bases aligned with MQ > 20 (blue) and % usable bases on target (red)")
plt


plt <- ggplot(data = cut_off_subset) +
  geom_point(aes(x = sample_name, y = 1 - PCT_EXC_MAPQ), color = "blue", alpha = 0.5) +
  geom_point(aes(x = sample_name, y = PCT_USABLE_BASES_ON_TARGET), color = "red") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(paste0("samples at ", ds.cutoff.intercept)) + 
  ylab("% bases aligned with MQ > 20 (blue) and % usable bases on target (red)") + 
  ylim(c(0,1))
plt







# get a summary of total reads in this subset
summ.cutoff.reads <- summary(cut_off_subset$TOTAL_READS)

# obtain all metrics for reads within IQR of this range
threshold.check <- QC[QC$TOTAL_READS > summ.cutoff.reads["1st Qu."] &
                             QC$TOTAL_READS < summ.cutoff.reads["3rd Qu."],]

# re-check for all 




plt <- ggplot() +
  # geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_boxplot(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_USABLE_BASES_ON_TARGET, group = downsamplingRatio), fill = "grey", alpha = 0.5) + 
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_USABLE_BASES_ON_TARGET, size = downsamplingRatio), color = "red", alpha = 0.25) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = PCT_USABLE_BASES_ON_TARGET, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Downsampling ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
plt






# for each 


# from the plot below infer expected mean target coverage for the given target bed file
plt <- ggplot() +
  # geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_boxplot(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = MEAN_TARGET_COVERAGE, group = downsamplingRatio), fill = "grey", alpha = 0.5) + 
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = MEAN_TARGET_COVERAGE, size = downsamplingRatio), color = "red", alpha = 0.25) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = MEAN_TARGET_COVERAGE, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Downsampling ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))  + scale_y_log10()

plt

plt <- ggplot() +
  # geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = 1-PCT_EXC_MAPQ, size = downsamplingRatio), color = "red", alpha = 0.5) +
  geom_boxplot(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = HET_SNP_SENSITIVITY, group = downsamplingRatio), fill = "grey", alpha = 0.5) + 
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = HET_SNP_SENSITIVITY, size = downsamplingRatio), color = "red", alpha = 0.25) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -TOTAL_READS), y = HET_SNP_SENSITIVITY, group = sample_name), color = "blue", alpha = 0.5) + 
  xlab("Downsampling ratio") +  theme(axis.text.x = element_text(angle = 90, hjust = 1))

plt

# calculate variance for the each downsampling



# identify the downsampling ratio where variance increases

# report the following plot (plot2)
plt <- ggplot() +
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = expected_total_reads, size = downsamplingRatio), color = "red", alpha = 0.5) +
  # geom_line(data = QC, aes(x = 1-downsamplingRatio, y = expected_total_reads, group = sample_name), color = "red", alpha = 0.5) +
  geom_point(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = TOTAL_READS), color = "blue", alpha = 0.25) +
  geom_line(data = QC, aes(x = reorder(downsamplingRatio, -downsamplingRatio), y = TOTAL_READS, group = sample_name), color = "blue", alpha = 0.25) +
  ylab("Total reads") + xlab("Downsampling ratio")
plt


# identify the total reads from downsampling ratio obtained from plot1




# since no duplicates were marked so ignoring the duplicate pct check


# read median target coverage 

# read sd target coverage 

# plot downsampling ratio Vs above 

# plot total reads vs above

# bootstrapping method to tune the optimal solution



