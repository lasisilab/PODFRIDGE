#Create a synthetic dataset with 1100 alleles (604,400 pairwise combinations) from our pilot test data
#This is not a randomly resampled dataset, it just repeats markers with new names in order to test function run times over more rows

#Create a synthetic dataset with 1100 alleles (604,400 pairwise combinations) from our test data

set.seed(444)
require(dplyr)

df_allelefreq <- fread("data/df_allelefreq_combined.csv")
df_allelefreq$population<-as.factor(df_allelefreq$population)
df_allelefreq$marker<-as.factor(df_allelefreq$marker)
t2<-round(1100/length(unique(df_allelefreq$marker)))

#syn_data<-df_allelefreq %>%
 # slice(rep(1:n(), each = t2))
syn_data<-do.call("rbind", replicate(t2, df_allelefreq, simplify = FALSE))

syn_data<-syn_data %>%
  group_by(allele,population,marker) %>%
  mutate(COUNTER = 1:n()) %>%
  ungroup()

syn_data$marker<-paste(syn_data$COUNTER,syn_data$marker,sep="_")
syn_data<-syn_data %>%
  select(-COUNTER)

fwrite(syn_data,"data/syn_data.csv")
