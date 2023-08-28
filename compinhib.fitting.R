rm(list = ls())
library(dplyr)
library(tidyr)
library(bayesplot)
library(rstan)
library(StanHeaders)
library(ggplot2)
parallel::detectCores()
library(rstudioapi)

setwd("~/MCSim_under_R")
source("MCSim/function.R")

model <- "crosstalk_compinhib.model.R" # the model file put in the model folder
input <- "crosstalk_compinhib.in.R" # the input file put in the infile folder

Modeling_dir = "Crosstalk"

makemcsim(model, dir = Modeling_dir) # 
# it is good to wait for 5 seconds inbetween running different chains
set.seed(7034)
jobRunScript(workingDir = getwd(), path = "Crosstalk/crosstalkfit.R", importEnv = T, exportEnv = "job_1")
set.seed(70004)
jobRunScript(workingDir = getwd(), path = "Crosstalk/crosstalkfit.R", importEnv = T, exportEnv = "job_2")
set.seed(728400)
jobRunScript(workingDir = getwd(), path = "Crosstalk/crosstalkfit.R", importEnv = T, exportEnv = "job_3")
set.seed(345400)
jobRunScript(workingDir = getwd(), path = "Crosstalk/crosstalkfit.R", importEnv = T, exportEnv = "job_4")
set.seed(43532)
jobRunScript(workingDir = getwd(), path = "Crosstalk/crosstalkfit.R", importEnv = T, exportEnv = "job_5")
set.seed(1451)
jobRunScript(workingDir = getwd(), path = "Crosstalk/crosstalkfit.R", importEnv = T, exportEnv = "job_6")

sims <- mcmc_array(list(job_1$out,job_2$out,job_3$out,job_4$out,job_5$out, job_6$out))  #,

endparms = length(names(job_1$out))-3
parms_name_1 <- paste((names(job_1$out)[2:endparms]),sep = ",")

# Files = list.files(path = paste0(getwd()),pattern = "^OSRfittingNrf2.*\\.out$",full.names=TRUE)
# myfiles = lapply(Files, read.delim)
# sims <- mcmc_array(list(myfiles[[1]],myfiles[[4]],myfiles[[5]],myfiles[[6]]),start_sampling = 0)
# endparms = length(names(myfiles[[1]]))-3

str <- ceiling(nrow(sims)/2) + 1
end <- nrow(sims)
j <- c(str:end)
length(j)
color_scheme_set("mix-blue-red")

#Trace and density plots of MCMC chains for specific parameters of interest

mcmc_trace(sims[j,,], pars = parms_name_1) + ggplot2::scale_color_discrete()

monitor_file = (monitor(sims[,,parms_name_1], digit=4))

mcmc_dens_overlay(x = sims[j,,], pars = parms_name_1)

mcmc_dens_overlay(x = sims[j,,], pars = "LnData")

dataf = data.frame((monitor_file))


write.csv(dataf[ ,c("mean","sd","X2.5.","X97.5.","Rhat")], file = paste0(getwd(),"/",Modeling_dir, "/compinhib_crosstalk_Srxn1.csv"))

X <- sims[j,,] %>% matrix(nrow = length(j)*6)

parms_name_2 <- paste((names(job_1$out)),sep = ",")

#####################

colnames(X) = parms_name_2

head(X)

#save the posterior parameter distribution file in your current directory

write.table(X, file = paste0(getwd(),"/",Modeling_dir, "/", Modeling_dir,"crosstalk_post_compinhibSrxn1.out"), row.names = F, sep = "\t")

#read the posterior parameter distribution file in your current directory

posteriordisfile <- read.delim(paste0(getwd(),"/",Modeling_dir, "/", Modeling_dir,"crosstalk_post_compinhibSrxn1.out"), header = TRUE, sep = "")

head(posteriordisfile)


# Below is the code to create simulation file using the posteriordisfile

params = data.frame(colnames(posteriordisfile))
colnames(params) = "parametersname"
file = params %>% separate(col = parametersname, c("Parms", "Subj", "levels"),sep = "([.])",fill = "right")
drug = c("Sulforaphane")
drugs = drug[1]
level = ifelse(drugs == "Sulforaphane",1, "NA") #here no role of level as this is non heirarchical methods 
a = file %>% filter(Subj ==1)
b = file %>% filter(Subj ==1 & levels == "")
c = file %>% filter(Subj ==1&levels ==level)
#c$levels = paste0(c$levels, ".")      #to add a dot as it is required  #required when you run heirarchical simulation
d = c(Parmshier = "iter",Parms = "iter", Subj = "", levels = " ")
e = c(Parmshier = "LnPrior",Parms = "LnPrior", Subj = "", levels = " ")
f = c(Parmshier = "LnData",Parms = "LnData", Subj = "", levels = " ")
g = c(Parmshier = "LnPosterior",Parms = "LnPosterior", Subj = "", levels = " ")

file1 = rbind(c,b)
file2 =distinct(file1,Parms,Subj, .keep_all= TRUE) %>% unite(.,"parmshier", Parms:levels, sep= ".", remove = FALSE) %>% rbind(d,.,e,f,g)
filename = file2$parmshier  #
randomdata = sample_n(posteriordisfile, 500)
hierSul_setpoint = randomdata %>% select(filename) %>% write.table(., file = paste0(getwd(),"/", Modeling_dir,"/",drugs, Modeling_dir,"setSim", ".out"), row.names = F, sep = "\t")
filesetout <- read.delim(paste0(getwd(),"/", Modeling_dir,"/",drugs, Modeling_dir,"setSim", ".out"), header = TRUE, sep = "")

#check if names are matching or not
colnames(filesetout)

Setsim = paste0(drugs, Modeling_dir,"setSim", ".out")
#chemical specific inputs and outputs
input = paste0(drugs, Modeling_dir,"set.in", ".R")  # this file will contain the chemical inputs. 
output = paste0(drugs, Modeling_dir,"set", ".out")  #check this output file name if it is the same or not

paramsSet = paste(file2$Parms[-c(1,length(file2$Parms)-2,(length(file2$Parms)-1), length(file2$Parms))], collapse = ",") 


# Read the experimental data file

File= read.delim(paste0(getwd(),"/",Modeling_dir, "/", "datafileA20.out"),header = TRUE, sep = "")
#File_A20 = read.delim(paste0(getwd(),"/",Modeling_dir, "/", "datafileA20.out"),header = TRUE, sep = "")



#Extract the information from experimental data

myfiles1 = File
myfiles1$replID = myfiles1$variable    
myfiles1$variable = 'A20'

myfiles1 = myfiles1%>%
  filter(dose_uM > 10)


# myfiles2 = File_A20
# myfiles2$replID = myfiles2$variable
# myfiles2$variable = 'A20'


Timeid = paste0(unique(myfiles1$timeID, incomparables = FALSE),sep = ",",collapse = " ") %>% gsub(",$", "", .)  #gsub to remove the last comma
#Timeid = paste0(seq(0,1320, by = 10),sep = ",",collapse = " ") %>% gsub(",$", "", .)  #gsub to remove the last comma
proteins = unique(myfiles1$variable, incomparables = FALSE)  # variables data that you want to fit
Final_Rdose = unique(myfiles1$dose_uM, incomparables = FALSE) # dose or stress amount

#Timeid1 = paste0(unique(myfiles2$timeID, incomparables = FALSE),sep = ",",collapse = " ") %>% gsub(",$", "", .)  #gsub to remove the last comma
#Timeid = paste0(seq(0,1320, by = 10),sep = ",",collapse = " ") %>% gsub(",$", "", .)  #gsub to remove the last comma
#proteins1 = unique(myfiles2$variable, incomparables = FALSE)  # variables data that you want to fit
#Final_Rdose1 = unique(myfiles2$dose_uM, incomparables = FALSE) # dose or stress amount
#myfiles1$replID = myfiles1$variable


### parameters in the output file need to be in simulation file with inbetween space
#gsub(",([A-Za-z])", ", \\1", .)  #to add a space after comma 

setsimparm = paste0( 'SetPoints ','(', '"', output , '"',',','"', Setsim ,'"', ',' , 0 , ',', paramsSet,')',';') %>% gsub(',([A-Za-z])', ', \\1', .)


simulationprint = function(dose, biomarkers, Timeid) {
  variables = ""
  for (i in 1:length(biomarkers)) {
    if (biomarkers[i] != "") {
      variable = paste0("Print(", biomarkers[i], ",", Timeid, ");")
      variables = paste(variables, variable, sep = "\n")
    }
  }
  
  return(paste0("Simulation", "{", "\n", "Dosing = ", dose, ";", "\n", variables, "\n", "}"))
}

#remove the input file if alread exist and create new one
file.remove(paste0(getwd(),"/",Modeling_dir, "/",input))

f = file(paste0(getwd(), "/", Modeling_dir, "/", input), open = 'a')
write(setsimparm, f)

for (i in 1:length(Final_Rdose)) {
  write(simulationprint(dose = Final_Rdose[i], biomarkers = proteins, Timeid), f)
}


close(f)

#Monte Carlo
MC_inhib = read.delim("MC.default.out")
head(MC_inhib)

proteins = c("A20")

vars <- names(MC_inhib)
Simulationfile = data.frame()


for (i in 1:length(proteins)){

  finder =  grep(pattern = paste0("^",proteins[i],collapse = ""),vars)
  file <- apply(MC_inhib[finder], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
  x = data.frame(file)
  X_mean <- apply(MC_inhib[finder],2, mean)   # can apply for all the doses
  x$mean = X_mean
  x$variable = proteins[i]
  #x$dose = Final_Rdose[j]
  Simulationfile <- rbind(Simulationfile,x)
}


X_setpts <- mcsim(model, input, Setsim = Setsim, dir = Modeling_dir)
invisible(file.remove(Setsim))    # to remove the output file
X_setpts1 = X_setpts
head(X_setpts1)
vars <- names(X_setpts1)

# Keap1modified = grep(pattern = paste0("^",proteins[4],collapse = ""),vars)
# names(X_setpts1)[Keap1modified] = "modifiedKeap1"

#updated the vars again
vars <- names(X_setpts1)
Simulationfile = data.frame()

for (i in 1:length(proteins)){

  finder =  grep(pattern = paste0("^",proteins[i],collapse = ""),vars)
  file <- apply(X_setpts1[finder], 2, quantile,  c(0.5, 0.025, 0.975)) %>% t()
  x = data.frame(file)
  X_mean <- apply(X_setpts1[finder],2, mean)   # can apply for all the doses
  x$mean = X_mean
  x$variable = proteins[i]
  #x$dose = Final_Rdose[j]
  Simulationfile <- rbind(Simulationfile,x)
  }



nrow(Simulationfile)

#############################
#prepration of dosing matrix

simul = length(Final_Rdose) #if we include dose zero to check the steady state
biomarker = length(proteins)
Times = unique(myfiles1$time) #experimental data time point
Simulationfile$Simulation = (rep(paste0(rep("Simu",simul),1:simul),each = length(Times)*biomarker))  # as we have 3 variables
Simulationfile$dose_uM = as.numeric(rep(paste0(rep(Final_Rdose)),each = length(Times)))
Simulationfile$Chemical = rep(paste0("NFT"),each = nrow(Simulationfile))
Simulationfile$Times <- Times
#Simulationfile$variables = (rep(paste0(rep(proteins)),each = length(Times)*simul))
colnames(Simulationfile) <- c("median", "LCL", "UCL","mean","variables","Simulations","dose_uM","Chemicals","Times")

Expdata = myfiles1
#Expdata1 = myfiles2


nrow(Expdata)
head(Simulationfile)

# plot_function = function (variable1) {
#   Simulation = Simulationfile %>% filter(.,variables == variable1)
#   #Simulation = Simulationfile %>% filter(.,variables == variable & Simulations != "Simu1")
#   exp_data = Expdata #%>% filter(variable == variable2)
#   
#   ggplot() +
#     geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_mean"),size=1) +
#     geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1) +
#     geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p95"),size=1) +
#     geom_point(data =exp_data,aes(x = timeID, y = value,color= "Exp_mean"),size=1) +
#     
#     #geom_point(data =exp_data,aes(x = time, y = value,color= "Exp_mean"),size=3) +
#     
#     facet_wrap(~dose_uM) +
#     labs(x = "Time (hr)", y = "GFP expression", sec.x="First exposure (uM)",
#          sec.y="Second exposure (uM)",title = paste0(variable1)) +
#     scale_y_continuous(limits=c(0,3)) +
#     scale_color_manual(name = "type",
#                        breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p95", "Exp_mean"),
#                        values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p95" = "red","Exp_mean" = "black"))+
#     theme(legend.title = element_blank()) + theme(legend.position="bottom") +
#     theme(legend.text = element_text(size= 8,face="bold")) 
# }




plot_function = function (variable1,variable2) {
  Simulation = Simulationfile %>% filter(.,variables == variable1)
  #Simulation = Simulationfile %>% filter(.,variables == variable & Simulations != "Simu1")
  exp_data = Expdata %>% filter(variable == variable2)
  
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_mean"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1.5) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p95"),size=1) +
    geom_point(data =exp_data,aes(x = timeID, y = value,color= "Exp_mean"),size=1) +
    
    facet_wrap(~dose_uM)+ #, scales = "free") + 
    labs(x = "Time (min)", y = "GFP expression", sec.x="First exposure (uM)",title = "Competitive Inhibition Model - A20") +
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p95", "Exp_mean"),
                       values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p95" = "red","Exp_mean" = "black"))
    #theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    #theme(legend.text = element_text(size= 8,face="bold")) 
}

plot_function_1 = function (variable1) {
  Simulation = Simulationfile %>% filter(.,variables == variable1)
  #Simulation = Simulationfile %>% filter(.,variables == variable & Simulations != "Simu1")
  #exp_data = Expdata %>% filter(variable == variable2)
  
  ggplot() +
    geom_line(data =Simulation,aes(x = Times, y = median, color="Simulated_mean"),size=1.2) +
    geom_line(data =Simulation,aes(x = Times, y = LCL, color="Simulated_p2.5"),size=1.5) +
    geom_line(data =Simulation,aes(x = Times, y = UCL, color="Simulated_p95"),size=1) +
    #geom_point(data =exp_data,aes(x = timeID, y = value,color= "Exp_mean"),size=1) +
    
    facet_wrap(~dose_uM)+ #, scales = "free") + 
    labs(x = "Time (min)", y = "GFP expression", sec.x="First exposure (uM)",title = "Competitive Inhibition Model - A20") +
    # scale_y_continuous(limits=c(0,3)) +
    scale_color_manual(name = "type",
                       breaks = c("Simulated_mean","Simulated_p2.5","Simulated_p95"),
                       values = c("Simulated_mean" = "blue","Simulated_p2.5" = "green","Simulated_p95" = "red"))+
    theme(legend.title = element_blank()) + theme(legend.position="bottom") +
    theme(legend.text = element_text(size= 8,face="bold"))
}
proteins

# plot_function(variable1 = "Srxn1", variable2 = "Srxn1")
# plot_function_1(variable1 = "Srxn1")
# plot_function_1(variable1 = "Nrf2")
plot_function_1(variable1 = "A20")
#plot_function(variable1 = "A20", variable2 = "A20")

#plot_function(variable1 = "Srxn1")#, variable2 = "repmean")



max(Simulation)



