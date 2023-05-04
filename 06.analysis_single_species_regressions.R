library(MASS)
library(lme4)

veg = read.csv("data/veg01.csv")
names(veg)
veg$year = as.character(veg$year)

env1 = read.csv("data/env01.csv")
env2 = read.csv("data/env02.csv")
env3 = read.csv("data/env03.csv")

env3.1 = as.data.frame(scale(env3[, 3:18]))
env3.1 = cbind(env3[ , c(1,2)], env3.1)


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


veg_tripar = veg_agg[veg_agg$species == "ACOROS", ]
veg_tripar = merge(veg_tripar, all_combinations, all = T, by = c("Year", "node"))
veg_tripar$species = "TRIPAR"
veg_tripar$abundance = ifelse(is.na(veg_tripar$abundance), 0, veg_tripar$abundance)

tripar_all = merge(veg_tripar, env3, by.x = c("Year", "node"), by.y = c("year", "node"), all.y = T)
tripar_all = merge(tripar_all, veg_tot, by = c("Year", "node"))
tripar_all$Year = as.numeric(tripar_all$Year)

#install.packages("MuMln")
library(MuMIn)


global.model2 = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = tripar_all, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(global.model1)

dredge(global.model1, beta = "n" )

test.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = tripar_all, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(test.mod)
ranef(test.mod)


#### extra ####
holder = as.data.frame(matrix(ncol = 8, nrow = 6))
colnames(holder) = c("SPECIES", "(Intercept)", "MEAN_SM_5", "MEAN_ST_5", "VAR_SM_5_tc", "VAR_ST_5_tc", "PRED_ST_5", "PRED_SM_5")

x = fixef(acoros.mod)
holder$SPECIES[1] = "Acoros"
holder$`(Intercept)`[1]  = 1.77834922
holder$MEAN_ST_5[1]    = 0.11587856
holder$MEAN_SM_5[1]    = 0.27560879
holder$VAR_ST_5_tc[1] = -0.09719104

as.data.frame(fixef(kobmyo.mod))

holder$SPECIES[2] = "Kobmyo"
holder$`(Intercept)`[2] =       -0.7805812
holder$MEAN_SM_5[2]  =       1.6712088
holder$MEAN_ST_5[2]     =      0.4786322
holder$VAR_SM_5_tc[2]    =    -0.2779713
holder$PRED_SM_5[2]       =    0.2780989
holder$PRED_ST_5[2]        =   0.2550884

as.data.frame(fixef(descae.mod))

holder$SPECIES[3] = "Descae"
holder$`(Intercept)`[3] =0.72522076
holder$MEAN_SM_5[3]    =     -0.18431834
holder$MEAN_ST_5[3]    =     -0.05842698
holder$VAR_SM_5_tc[3]   =    -0.30948294
holder$PRED_ST_5[3]    =     -0.10664494

as.data.frame(fixef(tripar.mod))

holder$SPECIES[4] = "Tripar"
holder$`(Intercept)`[4] =         1.5899501
holder$VAR_SM_5_tc[4] =         0.2961049
holder$PRED_SM_5[4] =           0.1753672


as.data.frame(fixef(bisbis.mod))

holder$SPECIES[5] = "Bisbis" 
holder$`(Intercept)`[5] =         1.2036472
holder$MEAN_ST_5[5] =         -0.1505450
holder$MEAN_SM_5[5]  =         0.5379067

as.data.frame(fixef(artsco.mod))

holder$SPECIES[6] = "Artsco"
holder$`(Intercept)`[6]  =       1.2864593
holder$MEAN_ST_5[6]   =       -0.2819675
holder$MEAN_SM_5[6]   =        0.3411765
holder$PRED_ST_5[6]   =       -0.2706665
holder$VAR_ST_5_tc[6] =        0.1700265

## ACOROS - BLUE
acoros.dat = tripar_all
acoros.mod = glmer(abundance ~ MEAN_ST_5 + MEAN_SM_5  + VAR_ST_5_tc + (1|node), data = acoros.dat, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(acoros.mod)

## KOBMYO - RED
kobmyo.dat = tripar_all
kobmyo.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = kobmyo.dat, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(kobmyo.mod)

## DESCAE - GREEN
descae.dat = tripar_all
descae.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + PRED_ST_5 + (1|node), data = descae.dat, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(descae.mod)


## TRIPAR - YELLOW
tripar.dat = tripar_all
tripar.mod = glmer(abundance ~  VAR_SM_5_tc + PRED_SM_5 + (1|node), data = tripar.dat, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(tripar.mod)

## BISBIS - ORANGE
bisbis.dat = tripar_all
bisbis.mod = glmer(abundance ~ MEAN_ST_5 + MEAN_SM_5 + (1|node), data = bisbis.dat, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(bisbis.mod)

## SILACA - BLACK
silaca.dat = tripar_all
silaca.mod = glmer(abundance ~ MEAN_SM_5 + PRED_SM_5 + VAR_SM_5_tc + VAR_ST_5_tc + (1|node), data = silaca.dat, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(silaca.mod)


## ARTSCO
artsco.dat = tripar_all
artsco.mod = glmer(abundance ~ MEAN_ST_5 + MEAN_SM_5 + PRED_ST_5+ VAR_ST_5_tc + (1|node), data = artsco.dat, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
summary(artsco.mod)

png("/Users/Will/Desktop/mod_slopes.png", height = 9, width = 6, res = 200, units = "in")
par(mfrow = c(3,2), mar = c(2,2,2,2))

library(RColorBrewer)
palette()

## MEAN_ST
plot.new()
plot.window(xlim = c(0, 1), ylim = c(-2,2))
box(bty = "l", lwd = 2)
axis(1)
axis(2)
for(i in 1:nrow(holder)){
  if(!is.na(holder$MEAN_ST_5[i])){
    abline(a=holder$`(Intercept)`[i], b = holder$MEAN_ST_5[i], lwd = 2, col = palette("Dark2")[i])
  }
}


## MEAN SM
plot.new()
plot.window(xlim = c(0, 1), ylim = c(-2,2))
box(bty = "l", lwd = 2)
axis(1)
axis(2)

for(i in 1:nrow(holder)){
  if(!is.na(holder$MEAN_SM_5[i])){
    abline(a=holder$`(Intercept)`[i], b = holder$MEAN_SM_5[i], lwd = 2, col = palette("Dark2")[i])
  }
}

## VAR_ST
plot.new()
plot.window(xlim = c(0, 1), ylim = c(-2,2))
box(bty = "l", lwd = 2)
axis(1)
axis(2)

for(i in 1:nrow(holder)){
  if(!is.na(holder$VAR_ST_5_tc[i])){
    abline(a=holder$`(Intercept)`[i], b = holder$VAR_ST_5_tc[i], lwd = 2, col = palette("Dark2")[i])
  }
}

## VAR_SM
plot.new()
plot.window(xlim = c(0, 1), ylim = c(-2,2))
box(bty = "l", lwd = 2)
axis(1)
axis(2)

for(i in 1:nrow(holder)){
  if(!is.na(holder$VAR_SM_5_tc[i])){
    abline(a=holder$`(Intercept)`[i], b = holder$VAR_SM_5_tc[i], lwd = 2, col = palette("Dark2")[i])
  }
}

## PRED_ST
plot.new()
plot.window(xlim = c(0, 1), ylim = c(-2,2))
box(bty = "l", lwd = 2)
axis(1)
axis(2)

for(i in 1:nrow(holder)){
  if(!is.na(holder$PRED_ST_5[i])){
    abline(a=holder$`(Intercept)`[i], b = holder$PRED_ST_5[i], lwd = 2, col = palette("Dark2")[i])
  }
}


## PRED_ST
plot.new()
plot.window(xlim = c(0, 1), ylim = c(-2,2))
box(bty = "l", lwd = 2)
axis(1)
axis(2)

for(i in 1:nrow(holder)){
  if(!is.na(holder$PRED_SM_5[i])){
    abline(a=holder$`(Intercept)`[i], b = holder$PRED_SM_5[i], lwd = 2, col = palette("Dark2")[i])
  }
}

dev.off()

par(mfrow = c(3,2))
plot(acoros.dat$abundance,(acoros.dat$abundance * acoros.dat$MEAN_ST_5), xlab = "abundance", ylab = "weighted abundance", main = "mean soil temperature")
plot(acoros.dat$abundance,(acoros.dat$abundance * acoros.dat$MEAN_SM_5), xlab = "abundance", ylab = "weighted abundance", main = "mean soil moisture")
plot(acoros.dat$abundance,(acoros.dat$abundance * acoros.dat$VAR_ST_5_tc), xlab = "abundance", ylab = "weighted abundance", main = "var soil temperature")
plot(acoros.dat$abundance,(acoros.dat$abundance * acoros.dat$VAR_SM_5_tc), xlab = "abundance", ylab = "weighted abundance", main = "var soil moisture")
plot(acoros.dat$abundance,(acoros.dat$abundance * acoros.dat$PRED_ST_5), xlab = "abundance", ylab = "weighted abundance", main = "pred soil temperature")
plot(acoros.dat$abundance,(acoros.dat$abundance * acoros.dat$PRED_SM_5), xlab = "abundance", ylab = "weighted abundance", main = "pred soil moisture")

par(mfrow = c(1,1))
plot.new()
legend(0,1, legend = holder$SPECIES, lty = 1, lwd = 2, col = palette("Dark2")[1:6])


aggregate(xx,FUN = mean)

## functions ####

get.coefs = function(mod.object, species){
  x0 = coef(mod.object)
  x1 = as.data.frame(x0$node)
  x1$species = species
  
  x2 = aggregate(x1[1:ncol(x1) - 1], by = list(x1$species), FUN = mean)
  return(x2)
  
}

holder = as.data.frame(matrix(nrow = 0, ncol = 8))
colnames(holder) = c("Species","(Intercept)", "MEAN_SM_5", "MEAN_ST_5", "PRED_SM_5", "PRED_ST_5", 
                     "VAR_SM_5_tc", "VAR_ST_5_tc")

ex.dredge = function(species, df){
  global.mod1 = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = df, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
  x0 = dredge(global.model1, beta = "n" )
  holder[nrow(holder) +1, 2:8] = x0[1, 1:7]
  holder[nrow(holder), "Species"] = species
  return(holder)
  
}



species = c("ARTSCO", "ACOROS", "DESCAE", "KOBMYO", "TRIPAR", "BISBIS")

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




library(ggplot2)

AS.df = make.df("ARTSCO")
AR.df = make.df("ACOROS")
DC.df = make.df("DESCAE")
KM.df = make.df("KOBMYO")
TP.df = make.df("TRIPAR")
BB.df = make.df("BISBIS")

par(mfrow = c(2,3))
plot(AS.df$abundance ~ AS.df$MEAN_ST_5, pch = 19, col = alpha(palette("Dark2")[1], 0.3), bty = "L", xlab = "Mean Soil Temperature (5 cm)", ylab = "abundance", ylim = c(0, 100))
abline(a = holder[holder$Species %in% "ARTSCO", "(Intercept)"], b = exp(holder[holder$Species %in% "ARTSCO", "MEAN_ST_5"]), col = palette("Dark2")[1])

points(AR.df$abundance ~ AR.df$MEAN_ST_5, pch = 19, col = alpha(palette("Dark2")[2], 0.3))
abline(a = exp(holder[holder$Species %in% "ACOROS", "(Intercept)"]), b = exp(holder[holder$Species %in% "ACOROS", "MEAN_ST_5"]), col = palette("Dark2")[2])

points(DC.df$abundance ~ DC.df$MEAN_ST_5, pch = 19, col = alpha(palette("Dark2")[3], 0.3))
abline(a = exp(holder[holder$Species %in% "DESCAE", "(Intercept)"]), b = exp(holder[holder$Species %in% "ACOROS", "MEAN_ST_5"]), col = palette("Dark2")[3])

points(KM.df$abundance ~ KM.df$MEAN_ST_5, pch = 19, col = alpha(palette("Dark2")[4], 0.3))
abline(a = exp(holder[holder$Species %in% "KOBMYO", "(Intercept)"]), b = exp(holder[holder$Species %in% "ACOROS", "MEAN_ST_5"]), col = palette("Dark2")[4])


points(KM.df$abundance ~ AS.df$MEAN_ST_5, pch = 19, col = alpha(palette("Dark2")[4], 0.3))
abline(a = holder[holder$Species %in% "KOBMYO", "(Intercept)"], b = exp(holder[holder$Species %in% "ACOROS", "MEAN_ST_5"]), col = palette("Dark2")[4])


range(tripar_all$MEAN_SM_5)
MEAN_SM_5 = seq(0.06, 0.32, length.out = 10)
range(tripar_all$MEAN_ST_5)
MEAN_ST_5 = seq(-3.17, 4.63, length.out = 10)
range(tripar_all$VAR_SM_5_tc)
VAR_SM_5 = seq(1, 15, length.out = 10)
range(tripar_all$VAR_ST_5_tc)
VAR_ST_5 = seq(2, 18, length.out = 10)
range(tripar_all$PRED_SM_5)
PRED_SM_5 = seq(0.38, 0.89, length.out = 10)
range(tripar_all$PRED_ST_5)
PRED_ST_5 = seq(0.82, 0.91, length.out = 10)

predict(global.model1, list(MEAN_SM_5 = MEAN_SM_5, MEAN_ST_5 = MEAN_ST_5, VAR_SM_5 = VAR_SM_5, VAR_ST_5 = VAR_ST_5, PRED_SM_5 = PRED_SM_5, PRED_ST_5 = PRED_ST_5), type = "response", times = 10)
predict(global.model1)


try.this = function(species){
  x0 = make.df(species)
  mean.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + (1|node), data = x0, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
  all.mod = glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + (1|node), data = x0, family = poisson(link = "log"), na.action = "na.fail", control = glmerControl(optimizer ="bobyqa"))
  mod.list = list(mean.mod = mean.mod, all.mod = all.mod)
  return(mod.list)
  
  
}



df = make.df("KOBMYO")
x = try.this("KOBMYO")
summary(x$mean.mod)
summary(x$all.mod)
anova(x$mean.mod, x$all.mod)

par(mfrow = c(3,2))
##
plot(abundance ~ MEAN_SM_5, data = df, xlim = c(0, 0.3))
get.coefs(x$mean.mod, species = "DESCAE")
#abline(a = get.coefs(x$mean.mod, species = "DESCAE")[1, 2], b = exp(get.coefs(x$mean.mod, species = "DESCAE")[1, 3]))

y0 = ggpredict(x$mean.mod)
y1= as.data.frame(y0$MEAN_SM_5)
lines(predicted ~ x, data = y1, lwd = 2, lty = 2)

y0 = ggpredict(x$all.mod)
y1= as.data.frame(y0$MEAN_SM_5)
lines(predicted ~ x, data = y0, lwd = 2, lty = 2)

##
plot(abundance ~ MEAN_ST_5, data = df)
get.coefs(x$mean.mod, species = "DESCAE")
#abline(a = get.coefs(x$mean.mod, species = "DESCAE")[1, 2], b = exp(get.coefs(x$mean.mod, species = "DESCAE")[1, 3]))

y0 = ggpredict(x$mean.mod)
y1= as.data.frame(y0$MEAN_ST_5)
lines(predicted ~ x, data = y0, lwd = 2, lty = 2)

y0 = ggpredict(x$all.mod)
y1= as.data.frame(y0$MEAN_ST_5)
lines(predicted ~ x, data = y1, lwd = 2, lty = 2)

##
plot(abundance ~ VAR_SM_5_tc, data = df)
get.coefs(x$all.mod, species = "DESCAE")
#abline(a = get.coefs(x$mean.mod, species = "DESCAE")[1, 2], b = exp(get.coefs(x$mean.mod, species = "DESCAE")[1, 3]))

y0 = ggpredict(x$all.mod)
y1= as.data.frame(y0$VAR_SM_5_tc)
lines(predicted ~ x, data = y1, lwd = 2, lty = 2)


range(tripar_all$MEAN_SM_5)
MEAN_SM_5 = seq(0.06, 0.32, length.out = 15)
range(tripar_all$MEAN_ST_5)
MEAN_ST_5 = seq(-3.17, 4.63, length.out = 15)
range(tripar_all$VAR_SM_5_tc)
VAR_SM_5_tc = seq(1, 15, length.out = 15)
range(tripar_all$VAR_ST_5_tc)
VAR_ST_5_tc = seq(2, 18, length.out = 15)
range(tripar_all$PRED_SM_5)
PRED_SM_5 = seq(0.38, 0.89, length.out = 15)
range(tripar_all$PRED_ST_5)
PRED_ST_5 = seq(0.82, 0.91, length.out = 15)

new_data = cbind(MEAN_SM_5, MEAN_ST_5, VAR_SM_5_tc, VAR_ST_5_tc, PRED_SM_5, PRED_ST_5)
new_data = as.data.frame(new_data)
new_data$node = sample(c(10, 11, 12, 13, 14, 16, 17, 19, 20, 6, 7, 8, 9), 1)
x = predict(x$all.mod, as.data.frame(new_data), type = "response")
plot(x)

