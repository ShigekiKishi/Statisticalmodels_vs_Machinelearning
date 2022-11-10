## file processing

library(tidyverse)
library(magrittr)

# in this case, data files were placed in "dat_files" folder
file_list <- list.files(path="dat_files/")
file_list <- file_list %>% str_replace(pattern=".csv", "")
Nfile <- length(file_list) # 203

for(i in 1:Nfile){ #
  file_name <- file_list[i]
  dat <- read.csv(paste0("dat_files/", file_name, ".csv"), header=T, row.names = 1) %>% as.data.frame() %>% na.omit()
  pref_num <- dat$Pref %>% unique() %>% length()
  if(pref_num > 1){
    source("bayes_calc.R")
  }
}

