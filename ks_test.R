# ks_test.r

get_significance_label <- function(p_value) {
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")  # ns = not significant
  }
}


library(stringr)




###CG
data_CG <- read.table("sorted_Mp_WT_embryo_merged_trimmed_bismark_bt2.deduplicated.genome1.genome2.overlapped.CpG.w100.thallus05.gff.windowed_by_tak1_20190222_repeatmodeler_withdata.JW3.sorted.filtered.short.gff", header=FALSE)
data_CG <- data_CG %>%rename(chrp=V1,startp=V2,endp=V3,paternal=V4,chrm=V5,startm=V6,endm=V7,maternal=V8)
paternal_CG <- data_CG$paternal
maternal_CG <- data_CG$maternal

result_CG <- wilcox.test(paternal_CG, maternal_CG, paired = TRUE, alternative = "less")
print(result_CG)
significance_label_CG <- get_significance_label(result_CG$p.value)


###CHH
data_CHH <- read.table("sorted_Mp_WT_embryo_merged_trimmed_bismark_bt2.deduplicated.genome1.genome2.overlapped.CHH.w100.thallus05.gff.windowed_by_tak1_20190222_repeatmodeler_withdata.JW3.sorted.filtered.short.gff", header=FALSE)
data_CHH <- data_CHH %>%rename(chrp=V1,startp=V2,endp=V3,paternal=V4,chrm=V5,startm=V6,endm=V7,maternal=V8)

paternal_CHH <- data_CHH$paternal
maternal_CHH <- data_CHH$maternal

result_CHH <- wilcox.test(paternal_CHH, maternal_CHH, paired = TRUE, alternative = "less")
print(result_CHH)
significance_label_CHH <- get_significance_label(result_CHH$p.value)


###CHG
data_CHG <- read.table("sorted_Mp_WT_embryo_merged_trimmed_bismark_bt2.deduplicated.genome1.genome2.overlapped.CHG.w100.thallus05.gff.windowed_by_tak1_20190222_repeatmodeler_withdata.JW3.sorted.filtered.short.gff", header=FALSE)
data_CHG <- data_CHG %>%rename(chrp=V1,startp=V2,endp=V3,paternal=V4,chrm=V5,startm=V6,endm=V7,maternal=V8)
paternal_CHG <- data_CHG$paternal
maternal_CHG <- data_CHG$maternal


result_CHG <- wilcox.test(paternal_CHG, maternal_CHG, paired = TRUE, alternative = "less")
print(result_CHG)
significance_label_CHG <- get_significance_label(result_CHG$p.value)


## plots


# Reshape data for plotting
data_long_CG <- data_CG %>% gather(key="Genome", value="Methylation_Score", paternal, maternal) %>% mutate(Context = "CG")
data_long_CHH <- data_CHH %>% gather(key="Genome", value="Methylation_Score", paternal, maternal) %>% mutate(Context = "CHH")
data_long_CHG <- data_CHG %>% gather(key="Genome", value="Methylation_Score", paternal, maternal) %>% mutate(Context = "CHG")

combined_data_long <- bind_rows(data_long_CG, data_long_CHH, data_long_CHG)

# Boxplot

# Define a position slightly above the maximum Methylation_Score for the bar and asterisk
y_pos <- max(combined_data_long$Methylation_Score) * 1.05

significance_data <- data.frame(
  Context = c("CG", "CHG", "CHH"),
  label = c(significance_label_CG, significance_label_CHG, significance_label_CHH)
)

title_text <- "Methylation scores of TEs between parental genomes in the embryo"
wrapped_title <- str_wrap(title_text, width = 40)  # Adjust the width as needed


boxplot <- ggplot(combined_data_long, aes(x=Genome, y=Methylation_Score, fill=Genome)) +
  geom_boxplot(outlier.color="black", outlier.shape=16, outlier.size=1) +
  labs(title=title_text, subtitle="Filtered by thallus CG>0.5 windowed by TEs", x="Genome", y="Methylation Score") +
  annotate("segment", x = 1, xend = 2, y = y_pos, yend = y_pos, size=1) +
  geom_text(data=significance_data, aes(x=1.5, y=y_pos * 1.02, label=label), inherit.aes=FALSE) +
  theme_classic() +
  facet_wrap(~ Context, nrow=1)

print(boxplot)

# Violin plot
violin_plot <- ggplot(data_long, aes(x=Genome, y=Methylation_Score)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, outlier.shape=NA) + # Adding a boxplot inside the violin plot for clarity
  labs(title="Violin Plot of Methylation Scores", x="Genome", y="Methylation Score") +
  theme_minimal()

print(violin_plot)



