library(irace)
library(EvoGCP)
library(parallel)

# Generating problems
gen.problems <- function()
{
  for (i in 1:10) {
    P <- EvoGCP::problemset(50,90,2.5,"random")
    save("P", file = paste("problemset",i,".Rda",sep=""))
  }
}

# Testing irace
test.irace <- function()
{
  irace::checkIraceScenario(scenario, parameters)
}

irace.param <- 'name "" c (hdpso)
pop "" i (5,40)
w "" r (0,0.5)
c1 "" r (0,7)
c2 "" r (0,1)'

parameters <- irace::readParameters(text = irace.param)

# Build scenario
scenario                <- irace::defaultScenario()
scenario$seed           <- 42 # Seed for the experiment
scenario$targetRunner   <- "target.runner" # Runner function (def. below)
scenario$instances      <- sapply(1:10,FUN = function(x) paste("problemset",x,".Rda",sep=""))
scenario$instances.extra.params <- rep(100000,10)

scenario$maxExperiments <- 20000 # Tuning budget

# Number of cores to be used by irace (set with caution!)
nc                      <- floor(parallel::detectCores()/2) + 1
scenario$parallel <- nc

target.runner <- function(experiment, scenario)
{
  load(experiment$instance)

  res <- EvoGCP::tester(P, nfe = experiment$extra.params, solverpar = experiment$configuration)

  print(res$violation)
  return(list(cost = sum(res$violation)))
}

# Uncomment below to run the full training session
irace.output <- irace::irace(scenario, parameters)
saveRDS(irace.output, "RESULTS.rds")
file.copy(from = "./irace.Rdata", to = "./irace-tuning.Rdata")


