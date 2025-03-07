
# Load Required Libraries
library(dplyr)
library(furrr)
library(data.table)
library(future)
library(parallel)
library(ggplot2)
library(tidyr)

file_type<-"png" #option to toggle between pdf or png outputs, pdf being faster
# Read Command-Line Arguments
#args <- commandArgs(trailingOnly = TRUE)
slurm_job_id <- as.character(args[1])

output_dir <- file.path("output")

# load combined_lrs
#combined_lrs = fread(paste0("output/",slurm_job_id,"/summary_genotypes.csv"))
combined_lrs = fread(paste0("sim_combined_genotypes_focal.csv"))
combined_lrs$LR<-as.numeric(combined_lrs$LR)
combined_lrs$relationship_tested <- factor(combined_lrs$relationship_tested, levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))
combined_lrs$population_known<-as.factor(combined_lrs$population)
combined_lrs$population_tested<-as.factor(combined_lrs$population_tested)

# Function to calculate cut-off values for 1%, 0.1%, and 0.01% FPR
calculate_cutoffs <- function(input_df, fp_rates, hypothesis) {
  input_df<-input_df[input_df$relationship_known == "unrelated" & input_df$relationship_tested == hypothesis,]

  #LR assumes known and tested populations are the same, and that known and tested relationships are accurate
  #LR columns specifying a different predicted population assume that the predicted relationship is accurate
  #LR columns specifying a different predicted relationship assume that the predicted population is accurate

  #LR column should be assigned cutoffs based on the original code- calculate proportions where population is the tested population
  #LR_x_population columns should be assigned separate cutoffs grouped by the known population

  cutoffs1 <- input_df[input_df$population_known==input_df$population_tested,
                       list(
                         fixed_cutoff = 1.00,
                         cutoff_1 = quantile(LR, probs = 1 - fp_rates[1] / 100, na.rm = TRUE),
                         cutoff_0_1 = quantile(LR, probs = 1 - fp_rates[2] / 100, na.rm = TRUE),
                         cutoff_0_01 = quantile(LR, probs = 1 - fp_rates[3] / 100, na.rm = TRUE),
                         n_unrelated = .N),
                       by=c("population_tested", "loci_set")] #Define cutoffs for population_tested, loci set
  return(cutoffs)
}

calculate_proportions_exceeding_cutoffs <- function(input_df, cutoffs, hypothesis) {
  input_df<-input_df[input_df$relationship_tested == hypothesis & input_df$population_known==input_df$population_tested,]
  df_with_cutoffs <- left_join(input_df, cutoffs, by = c("population_tested", "loci_set"))
  df_with_cutoffs <- df_with_cutoffs [,list(
    population_tested=population_tested,
    relationship_tested = relationship_tested,
    relationship_known = relationship_known,
    loci_set = loci_set,
    exceeds_fixed_cutoff = LR > fixed_cutoff,
    exceeds_cutoff_1 = LR > cutoff_1,
    exceeds_cutoff_0_1 = LR > cutoff_0_1,
    exceeds_cutoff_0_01 = LR > cutoff_0_01
  )]

  df_no_cutoffs<-df_with_cutoffs[!df_with_cutoffs$relationship_tested == "unrelated" & input_df$population_known==input_df$population_tested,]
  proportions_exceeding <- df_no_cutoffs[,    list(
    proportion_exceeding_fixed = sum(exceeds_fixed_cutoff) / .N,
    proportion_exceeding_1 = sum(exceeds_cutoff_1, na.rm = TRUE) / .N,
    proportion_exceeding_0_1 = sum(exceeds_cutoff_0_1, na.rm = TRUE) / .N,
    proportion_exceeding_0_01 = sum(exceeds_cutoff_0_01, na.rm = TRUE) / .N,
    n_related = .N),
    by=c("population","population_tested","relationship_tested","loci_set", "relationship_known")]
  return(proportions_exceeding)
}


# make sure loci set is in the order we want
combined_lrs$loci_set_factor <- factor(combined_lrs$loci_set, levels=c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"))

# make sure relationships are in the order we want
combined_lrs$relationship_known_factor <- factor(combined_lrs$relationship_known, levels=c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))
combined_lrs$relationship_tested_factor <- factor(combined_lrs$relationship_tested, levels=c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))

# race labels (can be changed)
combined_lrs <- combined_lrs %>%
  mutate(population_known_label = case_when(population_known == "AfAm" ~ "African-American",
                                            population_known == "Asian" ~ "Asian",
                                            population_known == "Cauc" ~ "Caucasian",
                                            population_known == "Hispanic" ~ "Hispanic"),
         population_tested_label = case_when(population_tested == "AfAm" ~ "African-American",
                                             population_tested == "Asian" ~ "Asian",
                                             population_tested == "Cauc" ~ "Caucasian",
                                             population_tested == "Hispanic" ~ "Hispanic"))

summary_stats <- combined_lrs[,    list(
  mean_LR = mean(LR),   #LR assumes relationship tested and population tested are correct
  lower_95 = quantile(LR, 0.025),
  upper_95 = quantile(LR, 0.975)#,
),
by=c("relationship_tested","population_tested_label","loci_set_factor")]
summary_stats<-summary_stats[summary_stats$population_known_label==summary_stats$population_tested_label,]
color_palette_race = c("#AA4499", "#DDCC77", "#88CCEE", "#117733")

plot0<-ggplot(summary_stats, aes(x = loci_set_factor, y = mean_LR, group = population_known_label, color = population_known_label)) +
  geom_line(size = 2) +
  facet_wrap(~ relationship_tested,  ncol = 2) +
  scale_y_log10() +
  labs(
    title = "Mean LR Across Populations and Relationship Types",
    x = "Loci Set",
    y = "Combined LR",
    color = "Population"
  ) +
  scale_color_manual(values = color_palette_race) +
  theme_minimal(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )
if(!file_type=="png"){
  pdf(paste(output_dir, "/Mean_LR_Across_Populations_and_Relationship_Types.pdf", sep = ""))
  invisible(print(plot0))
  dev.off()
} else {
  ggsave(plot=plot0,path=paste0(output_dir, "/"),filename="Mean_LR_Distributions_Across_Populations_and_Relationship_Types.png", width = 16, height = 12)
}



#combined_lrs$LR = ifelse(combined_lrs$LR < 1e-32, 1e-32, combined_lrs$LR)
plota<-ggplot(combined_lrs, aes(x = relationship_tested, y = LR, fill = population_tested_label, color = population_tested_label)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  facet_wrap(~ loci_set_factor, scales = "fixed") +
  labs(
    title = "LR Distributions Across Populations and Relationship Types",
    x = "Relationship Tested",
    y = "LR",
    fill = "Population",
    color = "Population"
  ) +
  scale_color_manual(values = color_palette_race) +
  scale_fill_manual(values = color_palette_race) +
  theme_minimal(base_size = 24) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_log10()

if(!file_type=="png"){
  pdf(paste(output_dir, "/LR_Distributions_Across_Populations_and_Relationship_Types.pdf", sep = ""))
  invisible(print(plota))
  dev.off()
} else {
  #  setwd(output_dir)
  ggsave(plot=plota,path=paste0(output_dir, "/"),filename="LR_Distributions_Across_Populations_and_Relationship_Types.png", width = 16, height = 12)
}

gb_df = combined_lrs %>%
  group_by(relationship_known, relationship_tested, loci_set) %>%
  summarize(mean = mean(LR))

gb_df$indicator = ifelse(gb_df$mean == 0, "*", "")

ann_text1 <- tribble(
  ~relationship_tested, ~LR, ~relationship_known, ~loci_set, ~population_label,
  "parent_child", 0, "cousins", "expanded_20", "African-American",
  "parent_child", 0, "cousins", "supplementary", "African-American",
  "parent_child", 0, "cousins", "autosomal_29", "African-American",
  "parent_child", 0, "second_cousins", "expanded_20", "African-American",
  "parent_child", 0, "second_cousins", "supplementary", "African-American",
  "parent_child", 0, "second_cousins", "autosomal_29", "African-American",
  "parent_child", 0, "unrelated", "expanded_20", "African-American",
  "parent_child", 0, "unrelated", "supplementary", "African-American",
  "parent_child", 0, "unrelated", "autosomal_29", "African-American",
  "parent_child", 0, "unrelated", "identifiler_15", "African-American"
)
ann_text1$loci_set_factor <- factor(ann_text1$loci_set, levels=c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"))
ann_text1$relationship_known_factor <- factor(ann_text1$relationship_known, levels=c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))

combined_lrs$trunc <- ifelse(combined_lrs$LR == 0, 1e-16, combined_lrs$LR)
n = nrow(combined_lrs[which(relationship_known == "unrelated" & relationship_tested == "parent_child"),])

plotb<-ggplot(combined_lrs, aes(x = relationship_tested, y = trunc, fill = population_label)) +
  geom_boxplot(position = position_dodge(width = 0.9)) +
  facet_grid(relationship_known_factor ~ loci_set_factor, scales = "free_y") +
  labs(
    title = "LR Distributions for Relationship Comparisons",
    subtitle = paste0("(n=", n, " related/unrelated pairs)"),
    x = "Tested Relationship Type",
    y = "LR",
    fill = "Population"
  ) +
  #geom_text(data = ann_text1, aes(y = LR), label = "*", size = 10) +
  theme_minimal(base_size = 24) +
  scale_fill_manual(values = color_palette_race) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_y_log10()+
  coord_cartesian(ylim=c(1e-15,NA))

if(!file_type=="png"){
  pdf(paste(output_dir, "/LR_Distributions_for_Relationship_Comparisons.pdf", sep = ""))
  invisible(print(plotb))
  dev.off()
} else {
  #  setwd(output_dir)
  ggsave(plot=plotb,path=paste0(output_dir, "/"),filename="LR_Distributions_for_Relationship_Comparisons.png", width = 16, height = 12)
}

#Plots showing the distribution of likelihood ratios relative to the cutoffs for each hypothesis type (cousins, full siblings, half siblings, parent_child, and #second cousins).
#Here, we use the following cutoffs: false positive rate (FPR) of 0.01%, 0.1%, 1%, and a fixed cutoff of 1.
#Values that are not present in the plot correspond to a likelihood ratio of 0.

# color palettte
my_colors = colorRampPalette(c("#FFB000","#F77A2E","#DE3A8A","#7253FF","#5E8BFF"))(6)

hypotheses = c("parent_child","full_siblings","cousins","half_siblings","second_cousins")
plot_list = c()
i=1
for(hypothesis in hypotheses){
  cutoffs <- calculate_cutoffs(combined_lrs, c(1, 0.1, 0.01), hypothesis = hypothesis)
  combined_lrs_hyp_cutoffs = merge(combined_lrs, cutoffs, by = c("population","loci_set"), all.x = TRUE)
  combined_lrs_hyp_cutoffs$relationship_known_factor <- factor(combined_lrs_hyp_cutoffs$relationship_known, levels=c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))

  combined_lrs_hyp_cutoffs$loci_set_factor <- factor(combined_lrs_hyp_cutoffs$loci_set, levels=c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"))

  combined_lrs_hyp_cutoffs = combined_lrs_hyp_cutoffs[which(relationship_tested == hypothesis),]

  plot1 = ggplot(data = combined_lrs_hyp_cutoffs) +
    geom_density(aes(x = log(LR), fill =relationship_known_factor),alpha = 0.5) +
    geom_vline(aes(xintercept = log(fixed_cutoff), col = "Fixed cutoff (1)"),linetype = "dashed") +
    geom_vline(aes(xintercept = log(cutoff_0_1), col = "0.1% FPR"),linetype = "dashed") +
    geom_vline(aes(xintercept = log(cutoff_0_01), col = "0.01% FPR"),linetype = "dashed") +
    geom_vline(aes(xintercept = log(cutoff_1), col = "1% FPR"),linetype = "dashed") +
    ggtitle(paste0("Log Likelihood Ratios for the ", hypothesis,  " hypothesis")) +
    scale_fill_manual(values = my_colors, labels = c("parent_child","full_siblings","half_siblings","cousins","second_cousins","unrelated")) +
    scale_color_manual(name = paste0("Hypothesis: ",hypothesis), labels = c("0.1% FPR", "0.01% FPR","1% FPR","Fixed cutoff (1)"), values = c("red","blue","forestgreen","black")) +
    facet_grid(loci_set_factor~population_label, scales = "free") +
    theme(axis.text = element_text(size = 28), axis.title = element_text(size = 28)) +
    guides(fill=guide_legend(title="Known Relationship")) +
    xlab("Log Likelihood Ratio (LLR)") + ylab("Frequency") +
    theme_minimal(base_size = 24)

  if(!file_type=="png"){
    pdf(paste(output_dir, "/Log_Likelihood_Ratios_for_hypothesis.pdf", sep = ""))
    invisible(print(plot1))
    dev.off()
  } else {
    ggsave(plot=plot1,path=paste0(output_dir, "/"),filename="Log_Likelihood_Ratios_for_hypothesis.png", width = 16, height = 12)
  }
  plot2 = ggplot(data = combined_lrs_hyp_cutoffs) +
    geom_histogram(aes(x = log(LR), fill =relationship_known_factor),position="identity", alpha = 0.5) +
    geom_vline(aes(xintercept = log(fixed_cutoff), col = "Fixed cutoff (1)"),linetype = "dashed") +
    geom_vline(aes(xintercept = log(cutoff_0_1), col = "0.1% FPR"),linetype = "dashed") +
    geom_vline(aes(xintercept = log(cutoff_0_01), col = "0.01% FPR"),linetype = "dashed") +
    geom_vline(aes(xintercept = log(cutoff_1), col = "1% FPR"),linetype = "dashed") +
    ggtitle(paste0("Log Likelihood Ratios for the ", hypothesis,  " hypothesis")) +
    scale_fill_manual(values = my_colors, labels = c("parent_child","full_siblings","half_siblings","cousins","second_cousins","unrelated")) +
    scale_color_manual(name = paste0("Hypothesis: ",hypothesis), labels = c("0.1% FPR", "0.01% FPR","1% FPR","Fixed cutoff (1)"), values = c("red","blue","forestgreen","black")) +
    facet_grid(loci_set_factor~population_label, scales = "free") +
    theme(axis.text = element_text(size = 28), axis.title = element_text(size = 28)) +
    guides(fill=guide_legend(title="Known Relationship")) +
    xlab("Log Likelihood Ratio (LLR)") + ylab("Frequency") +
    theme_minimal(base_size = 24)


  if(!file_type=="png"){
    pdf(paste(output_dir, "/Log_Likelihood_Ratios_for_hypothesis2.pdf", sep = ""))
    invisible(print(plot2))
    dev.off()
  } else {
    ggsave(plot=plot2,path=paste0(output_dir, "/"),filename="Log_Likelihood_Ratios_for_hypothesis2.png", width = 16, height = 12)
  }

  plot_list[[i]] = plot1
  plot_list[[i+1]] = plot2

  i = i+2

}

plot_list[[1]]
plot_list[[2]]

plot_list[[3]]
plot_list[[4]]

plot_list[[5]]
plot_list[[6]]

plot_list[[7]]
plot_list[[8]]

plot_list[[9]]
plot_list[[10]]


#Below, for each hypothesis (parent_child, full siblings, half siblings, cousins, and second cousins),
#we plot the proportion of pairs with a likelihood ratio exceeding the cutoffs for a false positive rate (FPR) of
#0.01%, 0.1%, 1% and fixed cutoff of 1. Results are faceted by loci set and known relationship type on the x-axis.


hypotheses = c("parent_child","full_siblings","cousins","half_siblings","second_cousins")
hypotheses_label = c("parent_child","Full Siblings","Cousins","Half Siblings","Second Cousins")

plot_list = c()
i=1

for(hypothesis in hypotheses){
  cutoffs <- calculate_cutoffs(combined_lrs, c(1, 0.1, 0.01), hypothesis = hypothesis)
  proportion_exceeding = calculate_proportions_exceeding_cutoffs(input_df = combined_lrs, cutoffs = cutoffs, hypothesis = hypothesis)

  # race labels (can be changed)
  proportion_exceeding <- proportion_exceeding %>%
    mutate(population_label = case_when(population == "AfAm" ~ "African-American",
                                        population == "Asian" ~ "Asian",
                                        population == "Cauc" ~ "Caucasian",
                                        population == "Hispanic" ~ "Hispanic"))

  proportions_exceeding_cutoffs_long = proportion_exceeding %>%
    pivot_longer(
      cols = c(proportion_exceeding_1,proportion_exceeding_0_1,proportion_exceeding_0_01,proportion_exceeding_fixed),
      names_to = c("Cutoff"),
      values_to = "Value")

  proportions_exceeding_cutoffs_long$loci_set_factor <- factor(proportions_exceeding_cutoffs_long$loci_set, levels=c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"))

  proportions_exceeding_cutoffs_long = proportions_exceeding_cutoffs_long %>%
    mutate(relationship_tested_label = case_when(relationship_tested == "cousins" ~ "Cousins",
                                                 relationship_tested == "full_siblings" ~ "Full Siblings",
                                                 relationship_tested == "half_siblings" ~ "Half Siblings",
                                                 relationship_tested == "parent_child" ~ "Parent-Child",
                                                 relationship_tested == "second_cousins" ~ "Second Cousins"))

  proportions_exceeding_cutoffs_long = proportions_exceeding_cutoffs_long %>%
    mutate(cutoff_label = case_when(Cutoff == "proportion_exceeding_0_01" ~ "0.01% FPR",
                                    Cutoff == "proportion_exceeding_0_1" ~ "0.1% FPR",
                                    Cutoff == "proportion_exceeding_1" ~ "1% FPR",
                                    Cutoff == "proportion_exceeding_fixed" ~ "Fixed (1)"))

  proportions_exceeding_cutoffs_long$loci_set_factor <- factor(proportions_exceeding_cutoffs_long$loci_set, levels=c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"))

  # make sure relationships are in the order we want
  proportions_exceeding_cutoffs_long$relationship_known_factor <- factor(proportions_exceeding_cutoffs_long$relationship_known, levels=c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))

  #proportions_exceeding_cutoffs_long_parent_child = proportions_exceeding_cutoffs_long[proportions_exceeding_cutoffs_long$relationship_tested == "parent_child",]

  plot1 = ggplot(proportions_exceeding_cutoffs_long, aes(x = relationship_known_factor, y = Value, fill = population_label, color = population_label)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_grid(cutoff_label ~ loci_set_factor, scales = "fixed") +
    labs(
      title = "Proportions Exceeding Likelihood Cut-offs",
      subtitle = paste0(hypotheses_label[i], " hypothesis"),
      x = "Relationship Known",
      y = "Proportion Exceeding Cut-off",
      fill = "Population",
      color = "Population"
    ) +
    scale_color_manual(values = color_palette_race) +
    scale_fill_manual(values = color_palette_race) +
    theme_minimal(base_size = 24) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  plot_list[[i]] = plot1
  i = i+1

  if(!file_type=="png"){
    pdf(paste(output_dir, "/Proportions_Exceeding_Likelihood_Cut-offs.pdf", sep = ""))
    invisible(print(plot1))
    dev.off()
  } else {
    ggsave(plot=plot1,path=paste0(output_dir, "/"),filename="Proportions_Exceeding_Likelihood_Cut-offs.png", width = 16, height = 12)
  }
}
