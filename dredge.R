library(MASS)
library(lme4)
library(MuMIn)
library(ggeffects)

## Set up ####

veg = read.csv("data/veg01.csv")
names(veg)
veg$year = as.character(veg$year)

env1 = read.csv("data/env01.csv")
colnames(env1) = c("node", "year", "airtemp_avg", "MEAN_ST_5", "soiltemp_30cm_avg", 
                   "MEAN_SM_5", "soilmoisture_a_30cm_avg", "airtemp_sd", 
                   "VAR_ST_5_tc", "soiltemp_30cm_sd", "VAR_SM_5_tc", 
                   "soilmoisture_a_30cm_sd", "soilmoisture_a_30cm_dr", "soiltemp_30cm_dr", 
                   "soilmoisture_a_5cm_dr", "soiltemp_5cm_dr", "PRED_ST_5", 
                   "st_acf_30cm_14", "PRED_SM_5", "sm_acf_30cm_14")
env1 = env1[env1$node != 15, ]
env1 = env1[env1$node != 12 | env1$year != 2021, ]


env2 = read.csv("data/env02.csv")
env3 = read.csv("data/env03.csv")

veg = veg[veg$year != 2017, ]

veg$multi.hit_or_top.hit = ifelse(veg$multi.hit_or_top.hit == "top-hit-of-multihit-protocol", 
                                  "top-hit", veg$multi.hit_or_top.hit)
veg$multi.hit_or_top.hit = ifelse(veg$multi.hit_or_top.hit == "non-top-hit-of-multihit-protocol", 
                                  "multi-hit", veg$multi.hit_or_top.hit)

veg = veg[veg$multi.hit_or_top.hit != "present", ]
veg = veg[veg$abundance != 0.5, ]

str(veg_imp)
veg_imp = veg[, c("species", "abundance", "sensor_node", "year")]
veg_agg = aggregate(veg_imp$abundance, by = list(veg_imp$year, veg_imp$sensor_node, veg_imp$species), FUN = sum)
veg_tot = aggregate(veg_imp$abundance, by = list(veg_imp$year, veg_imp$sensor_node), FUN = sum)

colnames(veg_tot) = c("Year", "node", "total.hits")
colnames(veg_agg) = c("Year", "node", "species", "abundance")

all_combinations = as.data.frame(expand.grid(c(2018,2019,2020,2021), c(10, 11, 14, 18, 20, 21, 15, 16, 8, 12, 13, 19, 
                                                                       6, 7, 17L, 9, 1)))
colnames(all_combinations) = c("Year", "node")

species = c("ARTSCO", "ACOROS", "DESCAE", "KOBMYO", "TRIPAR", "BISBIS")

holder = as.data.frame(matrix(nrow = 0, ncol = 8))
colnames(holder) = c("Species","(Intercept)", "MEAN_SM_5", "MEAN_ST_5", "PRED_SM_5", "PRED_ST_5", 
                     "VAR_SM_5_tc", "VAR_ST_5_tc")

env1 = env1[c(1:23, 25:51), ]
env3 = env3[c(1:23, 25:51), ]

## functions ####

get.coefs = function(mod.object, species){
  x0 = coef(mod.object)
  x1 = as.data.frame(x0$node)
  x1$species = species
  
  x2 = aggregate(x1[1:ncol(x1) - 1], by = list(x1$species), FUN = mean)
  return(x2)
  
}

ex.dredge = function(species, df){
  global.mod1 = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = df, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
  x0 = dredge(global.mod1, beta = "n" )
  holder[nrow(holder) +1, 2:8] = x0[1, 1:7]
  holder[nrow(holder), "Species"] = species
  return(holder)
  
}

make.df = function(species){
  veg_tripar = veg_agg[veg_agg$species == species, ]
  veg_tripar = merge(veg_tripar, all_combinations, all = T, by = c("Year", "node"))
  veg_tripar$species = species
  veg_tripar$abundance = ifelse(is.na(veg_tripar$abundance), 0, veg_tripar$abundance)
  
  tripar_all = merge(veg_tripar, env3, by.x = c("Year", "node"), by.y = c("year", "node"), all.y = T)
  tripar_all = merge(tripar_all, veg_tot, by = c("Year", "node"))
  tripar_all$Year = as.numeric(tripar_all$Year)
  return(tripar_all)
  
}


## Do all species dredge ####

for(i in species){
  tmp = make.df(i)
  holder = ex.dredge(i, tmp)
}

par(mfrow = c(3,2))

plot.slopes.points = function(term){
  as.df = make.df("ARTSCO")
  ar.df = make.df("ACOROS")
  dc.df = make.df("DESCAE")
  km.df = make.df("KOBMYO")
  tp.df = make.df("TRIPAR")
  bb.df = make.df("BISBIS")
  spp.list = list(as.df = as.df, ar.df = ar.df, dc.df = dc.df, km.df = km.df, tp.df = tp.df, bb.df = bb.df)
  print(length(spp.list))
  
  plot.new()
  plot.window(xlim = c(range(as.df[, term])), ylim = c(0,100))
  box(bty = "l")
  axis(1)
  axis(2)
  
  for(i in 1:length(spp.list)){
    if(!is.na(holder[i, term])){
      points(spp.list[[i]]$abundance ~ spp.list[[i]][, term], pch = 19, col = alpha(palette("Dark2")[i], 0.5))
      abline(a = holder[i, "(Intercept)"], b = holder[i, term], col = palette("Dark2")[i], lwd = 2, untf = TRUE)
      
    }

    
  }
  
  
}

plot.slopes.points("MEAN_SM_5")
plot.slopes.points("MEAN_ST_5")
plot.slopes.points("VAR_SM_5_tc")
plot.slopes.points("VAR_ST_5_tc")
plot.slopes.points("PRED_SM_5")
plot.slopes.points("PRED_ST_5")


dredge.em = function(env){
  as.df = make.df("ARTSCO", env)
  ar.df = make.df("ACOROS", env)
  dc.df = make.df("DESCAE", env)
  km.df = make.df("KOBMYO", env)
  tp.df = make.df("TRIPAR")
  bb.df = make.df("BISBIS")
  spp.list = list(as.df = as.df, ar.df = ar.df, dc.df = dc.df, km.df = km.df, tp.df = tp.df, bb.df = bb.df)
  mod.list = list()
  
  
  
  for(i in 1:length(spp.list)){
    df = spp.list[[i]]
    global.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = df, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
    x0 = dredge(global.mod)
    x1 = get.models(x0, subset = 1)
    mod.list[[i]] = x1
    
 
  }
  return(mod.list)
  
}

daily_mods = dredge.em()
seaso_mods = dredge.em()


plot.prediction = function(term, mod.list) {

  as.df = make.df("ARTSCO")
  ar.df = make.df("ACOROS")
  dc.df = make.df("DESCAE")
  km.df = make.df("KOBMYO")
  tp.df = make.df("TRIPAR")
  bb.df = make.df("BISBIS")
  spp.list = list(as.df = as.df, ar.df = ar.df, dc.df = dc.df, km.df = km.df, tp.df = tp.df, bb.df = bb.df)
  
  
  plot.new()
  plot.window(xlim = c(range(as.df[, term])), ylim = c(-5, 5))
  box(bty = "l")
  axis(1, cex.lab = 2)
  axis(2)
  
  for(i in 1:length(mod.list)){
    tmp = mod.list[[i]]
    tmp = tmp[[1]]
    if(any(names(coef(tmp)$node) %in% term)){
      x2 = as.data.frame(ggpredict(tmp, terms = paste(term, "[sample=10]", sep = "")))
      print(x2)
      points(log(spp.list[[i]]$abundance + 0.01) ~ spp.list[[i]][, term], pch = 19, col = alpha(palette("Dark2")[i], 0.5))
      lines(x2[, 1], log(x2[, 2]), lwd = 2, col = palette("Dark2")[i])
    }   
  }
}

#paste(term, "[sample=6]", sep = "")

set.seed(777)
library(ggplot2)

png("output/mod_plots_ln_daily.png", height = 10, width = 8, res = 200, units = "in")
par(mfrow = c(3,2), mar = c(2,3,2,3), cex.axis = 1.3)
plot.prediction("MEAN_SM_5", daily_mods)
plot.prediction("MEAN_ST_5", daily_mods)  
plot.prediction("VAR_SM_5_tc", daily_mods) 
plot.prediction("VAR_ST_5_tc", daily_mods)
plot.prediction("PRED_SM_5", daily_mods)
plot.prediction("PRED_ST_5", daily_mods)
dev.off()


### change make.df env dataframe to env3!!!! 
png("output/mod_plots_ln_seaso.png", height = 10, width = 8, res = 200, units = "in")
par(mfrow = c(3,2), mar = c(2,3,2,3), cex.axis = 1.3)
plot.prediction("MEAN_SM_5", seaso_mods)
plot.prediction("MEAN_ST_5", seaso_mods)  
plot.prediction("VAR_SM_5_tc", seaso_mods) 
plot.prediction("VAR_ST_5_tc", seaso_mods)
plot.prediction("PRED_SM_5", seaso_mods)
plot.prediction("PRED_ST_5", seaso_mods)
dev.off()



par(mfrow = c(1,1))
plot.new()
legend(0,0.5, legend = species, col = palette("Dark2")[1:6], lwd = 2)

#######################################

try.this = function(species){
  x0 = make.df(species)
  mean.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + (1|node), data = x0, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
  all.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = x0, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
  mod.list = list(mean.mod = mean.mod, all.mod = all.mod)
  return(mod.list)
  
}

mod.improve = function(type){
  as.df = make.df("ARTSCO")
  ar.df = make.df("ACOROS")
  dc.df = make.df("DESCAE")
  km.df = make.df("KOBMYO")
  tp.df = make.df("TRIPAR")
  bb.df = make.df("BISBIS")
  spp.list = list(as.df = as.df, ar.df = ar.df, dc.df = dc.df, km.df = km.df, tp.df = tp.df, bb.df = bb.df)
  mod.list = list()
  
  for(i in 1:6){
    x0 = spp.list[[i]]
    if(type == "means"){
      mean.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + (1|node), data = x0, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
      mod.list[[i]] = mean.mod

    } else if(type == "all"){
      all.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = x0, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
      mod.list[[i]] = all.mod

    } else if (type == "means_var"){
        means_var.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + (1|node), data = x0, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
        mod.list[[i]] = means_var.mod
    } else if (type == "means_pred"){
        means_pred.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + PRED_SM_5 + PRED_ST_5 + (1|node), data = x0, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
        mod.list[[i]] = means_pred.mod
      }
  }
  
  return(mod.list)
  
}

mean.mod.list = mod.improve(type = "means")
mean_var.mod.list = mod.improve(type = "means_var")
mean_pred.mod.list = mod.improve(type = "means_pred")
all.mod.list = mod.improve(type = "all")


png("output/improve_fig_means.png", height = 7, width = 4, res = 200, units = "in")
par(mfrow = c(2,1), mar = c(2,3,2,3))
plot.prediction("MEAN_SM_5", mean.mod.list)
plot.prediction("MEAN_ST_5", mean.mod.list)
dev.off()

png("output/improve_fig_all.png", height = 10, width = 8, res = 200, units = "in")
par(mfrow = c(3,2), mar = c(2,3,2,3), cex.axis = 1.3)
plot.prediction("MEAN_SM_5", all.mod.list)
plot.prediction("MEAN_ST_5", all.mod.list)
plot.prediction("VAR_SM_5_tc", all.mod.list)
plot.prediction("VAR_ST_5_tc", all.mod.list)
plot.prediction("PRED_SM_5", all.mod.list)
plot.prediction("PRED_ST_5", all.mod.list)
dev.off()

x0 = anova(all.mod.list[[6]], mean.mod.list[[6]])
summary(mean.mod.list[[6]])
summary(all.mod.list[[6]])


coef(summary(mean.mod.list[[1]]))

evaluate.model = function(mod.list1, mod.list2, mod.list3, mod.list4){
  tmp0 = data.frame(matrix(nrow = 6, ncol = 10))
  colnames(tmp0) = c("Species", "mv_p", "mv_aic", "mv_addsig", "mp_p", "mp_aic", "mp_addsig", "ma_p", "ma_aic", "ma_addsig")
  
  spp.list = c("ARTSCO", "ACOROS", "DESCAE", "KOBMYO", "TRIPAR", "BISBS")
  
  for(i in 1:6){
    x0 = anova(mod.list1[[i]], mod.list2[[i]])
    x1 = anova(mod.list1[[i]], mod.list3[[i]])
    x2 = anova(mod.list1[[i]], mod.list4[[i]])
    
    tmp0[i, "Species"] = spp.list[i]
    
    tmp0[i, "mv_p"] = x0[2, 8]
    tmp0[i, "mv_aic"] = x0[2,2] - x0[1,2]
    print(which(coef(summary(mod.list2[[i]]))[c(4,5), 4] < 0.05))
    tmp0[i, "mv_addsig"] = toString(which(coef(summary(mod.list2[[i]]))[c(4,5), 4] < 0.05))
    
    tmp0[i, "mp_p"] = x1[2, 8]
    tmp0[i, "mp_aic"] = x1[2,2] - x1[1,2]
    tmp0[i, "mp_addsig"] = toString(which(coef(summary(mod.list3[[i]]))[c(4,5), 4] < 0.05))
    
    tmp0[i, "ma_p"] = x2[2, 8]
    tmp0[i, "ma_aic"] = x2[2,2] - x2[1,2]
    tmp0[i, "ma_addsig"] = toString(which(coef(summary(mod.list4[[i]]))[c(4,5,6,7), 4] < 0.05))
    
  }
  
  return(tmp0)
  
}

modsss = evaluate.model(mean.mod.list, mean_var.mod.list, mean_pred.mod.list, all.mod.list)

c("ARTSCO", "ACOROS", "DESCAE", "KOBMYO", "TRIPAR", "BISBIS")
anova(mean.mod.list[[1]], mean_var.mod.list[[1]], mean_pred.mod.list[[1]], all.mod.list[[1]])
anova(mean.mod.list[[2]], mean_var.mod.list[[2]], mean_pred.mod.list[[2]], all.mod.list[[2]])
anova(mean.mod.list[[3]], mean_var.mod.list[[3]], mean_pred.mod.list[[3]], all.mod.list[[3]])
anova(mean.mod.list[[4]], mean_var.mod.list[[4]], mean_pred.mod.list[[4]], all.mod.list[[4]])
anova(mean.mod.list[[5]], mean_var.mod.list[[5]], mean_pred.mod.list[[5]], all.mod.list[[5]])
anova(mean.mod.list[[6]], mean_var.mod.list[[6]], mean_pred.mod.list[[6]], all.mod.list[[6]])



write.csv(modsss, "/Users/Will/Desktop/output_table2.csv", row.names = FALSE)

png("/Users/Will/Desktop/semi_Fig3.1_2.png", height = 10, width = 8, res = 200, units = "in")
par(mfrow = c(3,2), mar = c(2,3,2,3), cex.axis = 1.3)
plot.prediction2("MEAN_SM_5", mean.mod.list)
plot.prediction2("MEAN_ST_5", mean.mod.list)
plot.new()
plot.new()
plot.new()
plot.new()
dev.off()

png("/Users/Will/Desktop/semi_Fig3.2_2.png", height = 10, width = 8, res = 200, units = "in")
par(mfrow = c(3,2), mar = c(2,3,2,3), cex.axis = 1.3)
plot.prediction2("MEAN_SM_5", mean_var.mod.list)
plot.prediction2("MEAN_ST_5", mean_var.mod.list)
plot.prediction2("VAR_SM_5_tc", mean_var.mod.list)
plot.prediction2("VAR_ST_5_tc", mean_var.mod.list)
plot.new()
plot.new()
dev.off()

png("/Users/Will/Desktop/semi_Fig3.3_2.png", height = 10, width = 8, res = 200, units = "in")
par(mfrow = c(3,2), mar = c(2,3,2,3), cex.axis = 1.3)
plot.prediction2("MEAN_SM_5", mean_pred.mod.list)
plot.prediction2("MEAN_ST_5", mean_pred.mod.list)
plot.new()
plot.new()
plot.prediction2("PRED_SM_5", mean_pred.mod.list)
plot.prediction2("PRED_ST_5", mean_pred.mod.list)
dev.off()

png("/Users/Will/Desktop/semi_Fig3.4_2.png", height = 10, width = 8, res = 200, units = "in")
par(mfrow = c(3,2), mar = c(2,3,2,3), cex.axis = 1.3)
plot.prediction2("MEAN_SM_5", all.mod.list)
plot.prediction2("MEAN_ST_5", all.mod.list)
plot.prediction2("VAR_SM_5_tc", all.mod.list)
plot.prediction2("VAR_ST_5_tc", all.mod.list)
plot.prediction2("PRED_SM_5", all.mod.list)
plot.prediction2("PRED_ST_5", all.mod.list)
dev.off()

plot.prediction2 = function(term, mod.list) {
  
  as.df = make.df("ARTSCO")
  ar.df = make.df("ACOROS")
  dc.df = make.df("DESCAE")
  km.df = make.df("KOBMYO")
  tp.df = make.df("TRIPAR")
  bb.df = make.df("BISBIS")
  spp.list = list(as.df = as.df, ar.df = ar.df, dc.df = dc.df, km.df = km.df, tp.df = tp.df, bb.df = bb.df)
  
  
  plot.new()
  plot.window(xlim = c(range(as.df[, term])), ylim = c(-5, 5))
  box(bty = "l")
  axis(1, cex.lab = 2)
  axis(2)
  
  for(i in 1:length(mod.list)){
    tmp = mod.list[[i]]
    if(any(names(which(coef(summary(tmp))[, 4] < 0.05)) %in% term)){
      x2 = as.data.frame(ggpredict(tmp, terms = paste(term, "[sample=8]", sep = "")))
      print(x2)
      points(log(spp.list[[i]]$abundance + 0.01) ~ spp.list[[i]][, term], pch = 19, col = alpha(palette("Dark2")[i], 0.5))
      lines(x2[, 1], log(x2[, 2]), lwd = 2, col = palette("Dark2")[i])
    }   
  }
}


x1 = test[[1]]
names(which(coef(summary(x1[[1]]))[, 4] < 0.05))


any(names(coef(tmp)$node) %in% "MEAN_SM_5")


#########################################################


all.spp = c("ACHMIL", "ACOROS", "AGOGLA", "ALLGEY", "ANDSEP", "ANTMED", 
            "ANTROS", "AREFEN", "ARTPAT", "ARTSCO", "BISBIS", "BISVIV", 
            "CALLEP", "CALPUR", "CAMROT", "CAMUNI", "CARALB", "CARCAP", "CARELY", 
            "CARHET", "CARILL", "CARNOV", "CARROS", "CARRUP", 
            "CARSCO", "CARSIC", "CASOCC", "CERARV", "CHIJAM", "CRERUN", "DANINT", 
            "DESCAE", "DODPUL", "DRACRA", "ELYSCR", "ELYTRA", 
            "ERIGLA", "ERINAN", "ERIPIN", "ERISIM", "ERYCAP", "FESBRA", 
            "GENALG", "GENAMA", "HELMOR", "JUNI_sp", "JUNTRI", "KOBMYO", 
            "LEWPYG", "Lichen", "LIGTEN", "LLOSER", "LUZSPI", "MERLAN", 
            "MINOBT", "MINRUB",  "NOCFEN", "OREALP", "PACCAN", "PACCRO", 
            "PARPUL", "PEDGRO", "PEDPAR", "PEDRAC", "PENWHI", 
            "PHLPUL", "POAALP", "POAARC", "POAGLA", "POLVIS", "POTDIV","POTGLA", "POTRUB", "POTSUB", 
            "PRIPAR", "RHOINT", "RHORHO", "SALGLA", "SALIX", "SALPLA", 
            "SAXRHO", "SEDLAN", "SELDEN", "SIBPRO", "SILACA", "SOLMUL", 
            "SOLSIM", "STELON", "STEUMB", "SWEPER", "TETACA", "TETGRA", "TONLYA", 
            "TONPYG", "TRIDAS", "TRIPAR", "TRISPI")

# holder = as.data.frame(matrix(nrow = 0, ncol = 8))
# colnames(holder) = c("Species","(Intercept)", "MEAN_SM_5", "MEAN_ST_5", "PRED_SM_5", "PRED_ST_5", 
#                      "VAR_SM_5_tc", "VAR_ST_5_tc")

for(i in all.spp){
  tmp = make.df(i)
  holder = ex.dredge(i, tmp)
}

tmp_holder = holder


holder$group = 1
holder = holder[holder$Species != "Rock", ]
holder = holder[43:84,]


table(is.na(holder$MEAN_SM_5))
table(is.na(holder$MEAN_ST_5))
table(is.na(holder$VAR_SM_5_tc))
table(is.na(holder$VAR_ST_5_tc))
table(is.na(holder$PRED_SM_5))
table(is.na(holder$PRED_ST_5))

selected.daily = c(8, 12, 14, 12, 9, 21)
selected.seasonal = c(11, 14, 14, 8, 7, 18)

nope = c(35, 30, 28, 30, 33, 21)


png("/Users/Will/Desktop/Auto_fig1.png", height = 3, width = 5, res = 450, units = "in")
par(mfrow = c(1,2), mar = c(2,2,2,2))
#barplot(c(45,45,45,45,45,45), ylim = c(0,48), names.arg = c("MEAN_SM_5", "MEAN_ST_5", "VAR_SM_5_tc", "VAR_ST_5_tc", "PRED_SM_5", "PRED_ST_5"), cex.names = 0.8, col = "dark red", main = "Number of species where model terms were selected")
barplot(selected.daily, ylim = c(0,45), add = F, col = alpha("dark gray", 0.5), main = "Daily Variability & Predictability", cex.main = 0.6)
box(bty = "l", lwd = 1.5)
barplot(selected.seasonal, ylim = c(0,45), add = F, col = alpha("dark gray", 0.5), main = "Seasonal Variability & Predictability", cex.main = 0.6)
box(bty = "l", lwd = 1.5)
dev.off()
