
library(rstan)
library(bayesplot)

setwd("~/MCSim_under_R")
source("MCSim/function.R")

model <- "crosstalk_compinhib.model.R" # the model file put in the model folder
input <- "crosstalk_compinhib.in.R"
#input <- "compinhib_montecarlo.in.R" # the input file put in the infile folder

Modeling_dir = "Crosstalk"


# Create the executable file
makemcsim(model = model, deSolve = F, dir = Modeling_dir)
out = mcsim(model = model, input = input, dir = Modeling_dir, parallel = T)
