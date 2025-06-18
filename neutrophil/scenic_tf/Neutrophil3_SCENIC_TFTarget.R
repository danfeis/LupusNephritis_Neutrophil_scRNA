## Extract TF-target result from SCENIC, and organize into a more readable table
## The table contains TF, Target, Importance of TF regulation of each target gene, and context of what database used

library(tidyverse)

## read reg.csv from SCENIC, which contains TF and their regulated modules
## process headers of the file
df <- read.csv("pySCENIC_FinalVersion/Output/reg.csv", skip = 1, header = TRUE, check.names = FALSE) # skipping the first row for duplication removal
colnames(df)[1:2] <- c("TF", "MotifID") # Rename the first two unnamed columns
df <- df %>% select(where(~ !all(is.na(.)))) # Remove unnecessary NA columns
df <- df[-1,] # remove the 1st row of NA

# Convert TF, Context, and TargetGenes into a dataframe
TF_Target <- df %>% 
    select(TF, Context, TargetGenes) %>% mutate(TargetGenes = str_extract_all(df$TargetGenes, "\\(([^)]+)\\)")) %>%
    unnest(TargetGenes) %>%
    mutate(TargetGenes = str_remove_all(TargetGenes, "[()'\"]")) %>%  # Remove brackets & quotes
    separate(TargetGenes, into = c("TargetGene", "Coef"), sep = ", ", convert = TRUE) %>%  # Separate into columns
    mutate(Coef = as.numeric(Coef)) %>% # Ensure Coef is numeric
    drop_na(Coef) # Remove rows where Coef is NA
    

write_tsv(TF_Target, file = 'Genelist/Genelist3_TFTarget_SCENIC_All.tsv')
