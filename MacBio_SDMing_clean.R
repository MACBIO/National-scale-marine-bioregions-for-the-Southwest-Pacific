# SDM IUCN bioregionalisation
# by Maria Beger, m.beger@leeds.ac.uk
# final version, 29 September 2019

rm(list=ls()) 
setwd("//<directory>")

require(mgcv)
require(splines)
require(ReporteRs)
library(sqldf)


# Load data: sites with sitecodes, and fish data (PA or abundances adjusted across methods)
  ## all fish records
All.Fish = read.csv("fish_data_4561sites.csv", header = T, sep=",", dec=".", stringsAsFactors=FALSE)
   ## all sites where fish were recorded with environmental variables
Allenv = read.csv("env_data_4561sites_final_model8.csv", header = T, 
                  sep=",", dec=".", stringsAsFactors=FALSE)

### subset fishes by family by observability (cryptics, families are covered by all data providers to a high level)

# 1. grab cryptic species from all and stick into P.only dataset
cryptics = c("Holocentridae", "Ophichthidae", "Muraenidae", "Carcharhinidae", "Gobiidae", "Blenniidae", "Callionymidae", "Pseudochromidae", "Apogonidae", "Cirrhitidae", "Stegostomadae", "Scombridae", "Congridae", "Monacanthidae", "Hemiramphidae", "Mugilidae", "Tripterygiidae", "Belonidae", "Dasyatidae", "Syngnathidae", "Synodontidae", "Synanceiidae", "Engraulidae", "Atherinidae", "Clupeidae", "Plesiopidae", "Sphyrnidae", "Scorpaenidae", "Echeneidae", "Priacanthidae", "Plotosidae", "Ephippidae", "Pholidichthyidae", "Pempheridae", "Opistognathidae", "Ginglymostomatidae", "Mobulidae", "Malacanthidae", "Kuhliidae", "Synanceiidae", "Gerreidae", "PLATYCEPHALIDAE", "Echeneidae", "Gobiesocidae", "Chanidae", "Centriscidae", "Caracanthidae", "Eleotridae", "Anarhichadidae", "Myliobatidae", "Bothidae", "Myliobatidae", "Kyphosidae", "Monacanthidae", "Tetraodontidae", "Belonidae", "Pinguipedidae", "Scorpaenidae", "Dasyatidae" )
Fish.cryptics = All.Fish[All.Fish$Family %in% cryptics,] 

# 3. grab Presense - Absence/ Abundance fishes
good.fish = c("Carangidae",  "Labridae","Lethrinidae", "Lutjanidae", "Mullidae", "Serranidae","Chaetodontidae", "Fistulariidae", "Acanthuridae", "Kyphosidae", "Scaridae", "Siganidae", "Balistidae", "Ephippidae", "Caesionidae", "Zanclidae", "Cheilidactylidae", "Oplegnathidae", "Pomacanthidae", "Pomacentridae", "Cirrhitidae", "Haemulidae", "Nemipteridae", "Microdesmidae", "Diodontidae", "Sphyraenidae", "Aulostomidae", "Ostraciidae")
Fish.PA = All.Fish[All.Fish$Family %in% good.fish,] 


# @@@@@@@@  Presence/ Absense analysis  @@@@@@@
## convert all data to PA
All.Fish$Ab_adjust = 1
Fish.PA$Ab_adjust = 1

## subset sites to all providers but not Living Oceans Foundation (who have a specific set of target species). 
# Sites are either all sites, or all sites MINUS LOF sites (for non-LOF groups).  
# Assumes that the PA species groups occurr everywhere. 
All_nonLOF_PA_s = All.Fish[!All.Fish$Method %in% c("LOF.transect"),]
All_nonLOF_PA_sYES = as.data.frame(unique(All_nonLOF_PA_s$SurveyID))
colnames(All_nonLOF_PA_sYES) = c("SurveyID")

# First, run the predictions for the LOF species subset (with all sites), 
# then run the rest of species with non.LOF sites to predict. 
# Because PA species/ families should be everywhere, we treat absenses 
# as real across the project region, ie use all env data.  

PAfisheslist = as.data.frame(unique(All.Fish$Species))
LOFPAfishes_list = read.csv("LOF_spp.csv", header = T, sep=",", dec=".", stringsAsFactors=FALSE)
LOFPAfishes_list = as.data.frame(LOFPAfishes_list$Scientific.name); colnames(LOFPAfishes_list) = c("LOFSpecies")
PAFishnoLOF_list =  sqldf("select * from PAfisheslist EXCEPT select * from LOFPAfishes_list")

# load all environmental data for nearshore cells (9x9km)
newdata = read.csv("AllInshoreCells_envPar_model8.csv", header = T, 
                   sep=",", dec=".", stringsAsFactors=FALSE) #predictor data

# build dataframes to put results
fits = data.frame(cbind(newdata$cellID))
colnames(fits)[1] = c("cellID")
fit_ses = data.frame(cbind(newdata$cellID))
colnames(fit_ses) = c("cellID")
rejects = data.frame(Species=as.character(), n = numeric(), stringsAsFactors=FALSE)
models = data.frame(Species=as.character(), Intercept = numeric(), calcite = numeric(), CHL_MEAN = numeric(), 
                    nitrate = numeric(), parmax = numeric(), parmean = numeric(), ph = numeric(), 
                    SST_MEAN = numeric(), SST_RANGE = numeric(), R2 = numeric(), df = numeric(), logLik = numeric(), AICc = numeric(), 
                    delta = numeric(), weight = numeric(),stringsAsFactors=FALSE)

List_m = list()  # Initial model list
List_avm = list()  # Initial model averaging list

# start building the Word Document for model outputs
doc = docx(title = "Pacific Nearshore Species Models - Final Run LOF")
doc = addParagraph( doc , "2nd run of LOF Species", Level = 1)

#### model LOF species first, with ALL sites 
for (iPAspp in 1:nrow(LOFPAfishes_list)) {  
    iSp = Fish.PA[Fish.PA$Species == LOFPAfishes_list[iPAspp,],]  # select the spp we want
    
    if (nrow(iSp) > 30) {
      idata = sqldf("select * from Allenv LEFT JOIN iSp USING(SurveyID)")  
      idata$Ab_adjust[is.na(idata$Ab_adjust)] = 0 
      imod = gam(Ab_adjust ~ s(calcite, bs="cr") + s(CHL_MEAN, bs="cr") + s(nitrate, bs="cr") + s(parmax,bs="cr") + s(parmean, bs="cr") + s(ph, bs="cr") + s(SST_MEAN, bs="cr") + s(SST_RANGE, bs="cr"), data = idata, family = binomial, method="REML", na.action = na.fail, select=T)

      
      List_m[[length(List_m)+1]] = list(imod)
      names(List_m) = sprintf("m%i", iPAspp)
      
      # predict into all cells
      p_model = predict.gam(imod, newdata, type = 'response', se = T)
      fits = data.frame(fits, p_model$fit); colnames(fits)[ncol(fits)] = paste0(LOFPAfishes_list[iPAspp,], "_m", iPAspp)
      fit_ses = data.frame(fit_ses, p_model$se.fit);  colnames(fit_ses)[ncol(fit_ses)] = paste0(LOFPAfishes_list[iPAspp,], "_m", iPAspp)
      
      # generate plot of the model
      name = paste0("./jpegs/", LOFPAfishes_list[iPAspp,], ".jpeg")
      jpeg(file = name)
      plot(imod, residuals = TRUE, rug=FALSE, col = "blue", lwd = 2, shade=TRUE, shade.col="lightblue", seWithMean = TRUE, scale=0, cex.lab=1.2, cex.axis=1.2, pages = 1)
      dev.off()
      
      #add outputs to word doc
      fishname = paste0(LOFPAfishes_list[iPAspp,], ", n = ", nrow(iSp), " observations") 
      sum.imod = summary(imod)
      doc = addTitle(doc, fishname, level = 1)
      table = FlexTable(data = sum.imod$s.table, add.rownames = TRUE)
      doc = addFlexTable (doc, table)
      doc = addImage(doc, name, width = 6, height = 6)
      } else {
      # add to a reject table as the spp is too rare to model
      rejects [iPAspp,] =  c(LOFPAfishes_list[iPAspp,],nrow(iSp))
    }
  }
  
writeDoc(doc, "PA_FishModelsLOF_YAY_endRUN.docx" )
write.csv(fits, "PA_predictionsLOF_YAY_endRUN.csv")
write.csv(fit_ses, "PA_half_errorsLOF_YAY_endRUN.csv")
write.csv(rejects, "PA_rejectsLOF_YAY_endRUN.csv2")



### all other PA species, with no-LOF sites
Allenv_noLOF = sqldf("select * from Allenv LEFT JOIN All_nonLOF_PA_sYES USING(SurveyID)")
doc = addParagraph( doc , "Non-LOF PA Species", Level = 1)

for (iPAspp in 1:nrow(PAFishnoLOF_list)) {
    iSp = Fish.PA[Fish.PA$Species == PAFishnoLOF_list[iPAspp,],]  # select the spp we want
  
  if (nrow(iSp) > 30) {
    idata = sqldf("select * from Allenv_noLOF LEFT JOIN iSp USING(SurveyID)")  
    idata$Ab_adjust[is.na(idata$Ab_adjust)] = 0 
    imod = gam(Ab_adjust ~ s(calcite, bs="cr") + s(CHL_MEAN, bs="cr") + s(nitrate, bs="cr") + s(parmax,bs="cr") + s(parmean, bs="cr") + s(ph, bs="cr") + s(SST_MEAN, bs="cr") + s(SST_RANGE, bs="cr"), data = idata, family = binomial, method="REML", na.action = na.fail, select=T)
    
    List_m[[length(List_m)+1]] = list(imod)
    names(List_m) = sprintf("m%i", iPAspp)
  
    p_model = predict.gam(imod, newdata, type = 'response', se = T)
    fits = data.frame(fits, p_model$fit); colnames(fits)[ncol(fits)] = paste0(PAFishnoLOF_list[iPAspp,], "_m", iPAspp)
    fit_ses = data.frame(fit_ses, p_model$se.fit);  colnames(fit_ses)[ncol(fit_ses)] = paste0(PAFishnoLOF_list[iPAspp,], "_m", iPAspp)
    
    name = paste0("./jpegs/", PAFishnoLOF_list[iPAspp,], ".jpeg")
    jpeg(file = name)
    plot(imod, residuals = TRUE, rug=FALSE, col = "blue", lwd = 2, shade=TRUE, shade.col="lightblue", seWithMean = TRUE, scale=0, cex.lab=1.2, cex.axis=1.2, pages = 1)
    dev.off()
    
    fishname = paste0(PAFishnoLOF_list[iPAspp,], ", n = ", nrow(iSp), " observations") 
    sum.imod = summary(imod)
    doc = addTitle(doc, fishname, level = 1)
    table = FlexTable(data = sum.imod$s.table, add.rownames = TRUE)
    doc = addFlexTable (doc, table)
    doc = addImage(doc, name, width = 6, height = 6)
    } else {
    # add something to a table saying the spp is out
    rejects [iPAspp,] =  c(PAFishnoLOF_list[iPAspp,], nrow(iSp) )
  }
}


## Write PAFishnoLOF things to file:
writeDoc(doc, "PA_nonLOF_FishModels_YAY_endRUN.docx" )
write.csv(fits, "PA_nonLOF_Fishpredictions_YAY_endRUN.csv")
write.csv(fit_ses, "PAnon_nonLOF_Fish_half_errors_YAY_endRUN.csv")
write.csv(rejects, "PAnon_2nd_lot_nonLOF_rejects_YAY_endRUN.csv")

  