library(tidyverse)

# parameter identifiability checks
sims.list <- readRDS("analysis-output/ms-simulation-simslist2021-02-08.rds")
sims.list <- as.data.frame(sims.list)

theme_set(theme_bw())

for (i in colnames(sims.list)){
  
  png(filename = paste0("analysis-output/parameter-identifiability/2021-02-08/",
                        i,"-","check.png"),
      width=4, height=3, units="in", res=600)
  
  print(ggplot(sims.list, aes(sims.list[,i])) +
          geom_density() +
          geom_hline(yintercept = 1, linetype = "dashed") +
          xlab(i))
  
  dev.off()
  
}
