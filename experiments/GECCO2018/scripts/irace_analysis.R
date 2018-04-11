library(irace)

names <- c("PSO", "ABC","CS","FFA")
files <- c(
  "raw/irace_tuning_PSO/irace.Rdata",
  "raw/irace_tuning_ABC/irace.Rdata",
  "raw/irace_tuning_CS/irace.Rdata",
  "raw/irace_tuning_FFA/irace.Rdata"
)

for (i in seq_along(names)) {
  load(files[i])
  assign(paste(names[i],".irace",sep=""), iraceResults)
}

best.configuration <- function(result.data, n = 1)
{
  best.config <- getFinalElites(iraceResults = result.data, n = 1)
  print("Best Configuration Details: ")
  print(best.config)
  id <- best.config$.ID.
  all.exp <- result.data$experiments[, as.character(id)]
  print("Best Configuration Results: ")
  all.exp[!is.na(all.exp)]
}

plot.convergence <- function(result.data)
{
  d <- apply(result.data$experiments, MARGIN = 2, FUN = mean, na.rm = T)
  plot(d, ylab = "Total Violations", xlab = "Configurations")
}

plot.elites <- function(result.data)
{
  best.exp <- result.data$experiments[,result.data$iterationElites]
  #conf <- gl(ncol(best.exp), nrow(best.exp), labels = colnames(best.exp))
  #pairwise.wilcox.test(as.vector(best.exp), conf, paired = TRUE, p.adj = "bonf")
  boxplot(best.exp)
}

plot.parameters <- function(result.data)
{
  parameterFrequency(result.data$allConfigurations, result.data$parameters)
}

