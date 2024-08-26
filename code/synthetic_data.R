#Create a synthetic dataset with 1100 alleles (604,400 pairwise combinations) from our test data to test function run times

set.seed(444)
require(dplyr)

df_allelefreq <- fread("data/df_allelefreq_combined.csv")
df_allelefreq$population<-as.factor(df_allelefreq$population)
df_allelefreq$marker<-as.factor(df_allelefreq$marker)

syn_data<-df_allelefreq %>%
  group_by(population) %>%
  do(sample_n(., 1100, replace = TRUE))

syn_data<-syn_data %>%
  group_by(population,marker) %>%
  mutate(COUNTER = 1:n()) %>%  
  ungroup()
syn_data$marker<-paste(syn_data$COUNTER,syn_data$marker,sep="_")
syn_data<-syn_data %>%
  select(-COUNTER)

  
