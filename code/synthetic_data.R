#Create a synthetic dataset with 1100 alleles (604,400 pairwise combinations) from our test data to test function run times

set.seed(444)
require(dplyr)

df_allelefreq <- fread("data/df_allelefreq_combined.csv")
df_allelefreq$population<-as.factor(df_allelefreq$population)
syn_data<-df_allelefreq %>%
  group_by(population) %>%
  do(sample_n(., 1100, replace = TRUE))
