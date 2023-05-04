library(MASS)
library(lme4)
library(MuMIn)
library(ggeffects)
library(lmerTest)

### DAILY
env1 = read.csv("data/env01.csv")
colnames(env1) = c("node", "year", "MEAN_SM_5", "VAR_SM_5_tc", "PRED_SM_5", "MEAN_ST_5", 
                   "VAR_ST_5_tc", "PRED_ST_5")
env1 = env1[env1$node != 15, ]
env1 = env1[env1$node != 12 | env1$year != 2021, ]



### SEASONAL
env3 = read.csv("data/env03.csv")

env1 = env1[c(1:23, 25:51), ]
env3 = env3[c(1:23, 25:51), ]

env5 = read.csv("data/env05.csv")

env5 = merge(env1[c("node", "year", "MEAN_SM_5")], env5, by = c("node", "year"))


matrix2df = function(species, matrix){
  t1 = matrix[, colnames(matrix) %in% species]
  t2 = data.frame(rownames(matrix), t1)
  colnames(t2) = c("name", "abundance")
  t2$species = species
  t2$node = matrix(unlist(strsplit(t2$name, split = "-")), ncol = 2, byrow = T)[, 1]
  t2$year = matrix(unlist(strsplit(t2$name, split = "-")), ncol = 2, byrow = T)[, 2]
  t2 = t2[, c(4,5,3,2)]
  return(t2)
}

make.df = function(veg, env){
  all = merge(veg, env, by.x = c("year", "node"), by.y = c("year", "node"), all.y = T)
  all$year = as.numeric(all$year)
  return(all)
  
}

dredge.em = function(env, matrix, logistic = F){
  as.df = make.df(matrix2df("TRIPAR", matrix), env)
  ar.df = make.df(matrix2df("ACOROS", matrix), env)
  dc.df = make.df(matrix2df("CARSCO", matrix), env)
  km.df = make.df(matrix2df("KOBMYO", matrix), env)
  tp.df = make.df(matrix2df("DESCAE", matrix), env)
  bb.df = make.df(matrix2df("SEDLAN", matrix), env)
  spp.list = list(as.df = as.df, ar.df = ar.df, dc.df = dc.df, km.df = km.df, tp.df = tp.df, bb.df = bb.df)
  mod.list = list()
  
  
  
  for(i in 1:length(spp.list)){
    df = spp.list[[i]]
    if(logistic == T){
      global.mod = glm(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + year, family = binomial(link = "logit") ,data = df, na.action = na.fail)
    } else if (logistic == F){
      df = df[df$abundance != 0, ]
      global.mod = lmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + year + (1|node), data = df, na.action = na.fail, REML = FALSE)
    }
    x0 = dredge(global.mod)
    #x1 = summary(model.avg(x0, subset = delta < 4, fit = T))
    x1 = get.models(x0, subset = 1)
    mod.list[[i]] = x1
    
    
  }
  return(mod.list)
  
}

test


daily.mods = dredge.em(env1, matrix)
seasonal.mods = dredge.em(env3, matrix)
daily.mods.logistic = dredge.em(env1, matrix2, logistic = T)
seasonal.mods.logistic = dredge.em(env3, matrix2, logistic = T)



plot.prediction = function(term, mod.list, env, matrix, logistic = F) {
  
  as.df = make.df(matrix2df("TRIPAR", matrix), env)
  ar.df = make.df(matrix2df("ACOROS", matrix), env)
  dc.df = make.df(matrix2df("CARSCO", matrix), env)
  km.df = make.df(matrix2df("KOBMYO", matrix), env)
  tp.df = make.df(matrix2df("DESCAE", matrix), env)
  bb.df = make.df(matrix2df("SEDLAN", matrix), env)
  spp.list = list(as.df = as.df, ar.df = ar.df, dc.df = dc.df, km.df = km.df, tp.df = tp.df, bb.df = bb.df)

  if(logistic == F){
    ylim = c(0,100)
  } else if(logistic == T){
    ylim = c(0,1)
  }
  
  plot.new()
  plot.window(xlim = c(range(as.df[, term])), ylim = ylim)
  box(bty = "l")
  axis(1, cex.lab = 2)
  axis(2)
  
  for(i in 1:length(mod.list)){
    tmp = mod.list[[i]]
    tmp = tmp[[1]]
    
    if(logistic == T){
      names = names(coef(tmp))
    }else if(logistic == F){
      names = names(coef(tmp)$node)
    }
    
    if(any(names %in% term)){
      df = spp.list[[i]]
      if(logistic == F){
        df = df[df$abundance != 0,]
      }
      termm = paste(term, "[all]", sep = "")
      #termm = paste(term, "[sample=50]", sep = "")
      print(termm)
      x2 = as.data.frame(ggpredict(tmp, terms = termm))
      points(df$abundance ~ df[, term], pch = 21, cex = 1.4,col = "black" ,bg = alpha(cols[i], 0.6))
      polygon(c(x2[,1],rev(x2[,1])),c(x2[,4],rev(x2[,5])),border=NA, col = alpha(cols[i], 0.2))
      lines(x2[, 1], x2[, 2], lwd = 2, col = cols[i])
    } else{
      df = spp.list[[i]]
      if(logistic == F){
        df = df[df$abundance != 0,]
      }
      points(df$abundance ~ df[, term], cex = 1.4,pch = 4, col = alpha("gray", 0.6))
      # x2 = as.data.frame(ggpredict(tmp, terms = paste(term, "[sample=8]", sep = "")))
      # lines(x2[,1], x2[,2], lwd = 2, lty = 2, col = alpha("gray", 0.5))
    }  
  }
}

df=make.df(matrix2df("TRIPAR", matrix2), env1)
glmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + year + (1|node), family  = "binomial", data = df, na.action = na.fail)

cols = c("red", "orange", "yellow", "green", "blue", "purple")

png("/Users/Will/Desktop/newFig31.png", height = 7, width = 4.67, res = 400, units = "in")
par(mfrow = c(3,2), mar = c(2,2,2,2), cex.axis = 1.1)
plot.prediction("MEAN_SM_5", daily.mods, env1, matrix)
plot.prediction("MEAN_ST_5", daily.mods, env1, matrix)
plot.prediction("VAR_SM_5_tc", daily.mods, env1, matrix)
plot.prediction("VAR_ST_5_tc", daily.mods, env1, matrix)
plot.prediction("PRED_SM_5", daily.mods, env1, matrix)
plot.prediction("PRED_ST_5", daily.mods, env1, matrix)
dev.off()

png("/Users/Will/Desktop/newFig32.png", height = 7, width = 4.67, res = 400, units = "in")
par(mfrow = c(3,2), mar = c(2,2,2,2), cex.axis = 1.1)
plot.prediction("MEAN_SM_5", seasonal.mods, env3, matrix)
plot.prediction("MEAN_ST_5", seasonal.mods, env3, matrix)
plot.prediction("VAR_SM_5_tc", seasonal.mods, env3, matrix)
plot.prediction("VAR_ST_5_tc", seasonal.mods, env3, matrix)
plot.prediction("PRED_SM_5", seasonal.mods, env3, matrix)
plot.prediction("PRED_ST_5", seasonal.mods, env3, matrix)
dev.off()

png("/Users/Will/Desktop/newFig41.png", height = 7, width = 4.67, res = 400, units = "in")
par(mfrow = c(3,2), mar = c(2,2,2,2), cex.axis = 1.1)
plot.prediction("MEAN_SM_5", daily.mods.logistic, env1,matrix2, logistic = T)
plot.prediction("MEAN_ST_5", daily.mods.logistic, env1,matrix2, logistic = T)
plot.prediction("VAR_SM_5_tc", daily.mods.logistic, env1,matrix2, logistic = T) # change to [sample=25]
plot.prediction("VAR_ST_5_tc", daily.mods.logistic, env1,matrix2, logistic = T)
plot.prediction("PRED_SM_5", daily.mods.logistic, env1,matrix2, logistic = T)
plot.prediction("PRED_ST_5", daily.mods.logistic, env1,matrix2, logistic = T)
dev.off()

png("/Users/Will/Desktop/newFig42.png", height = 7, width = 4.67, res = 400, units = "in")
par(mfrow = c(3,2), mar = c(2,2,2,2), cex.axis = 1.1)
plot.prediction("MEAN_SM_5", seasonal.mods.logistic, env3,matrix2, logistic = T)
plot.prediction("MEAN_ST_5", seasonal.mods.logistic, env3,matrix2, logistic = T)
plot.prediction("VAR_SM_5_tc", seasonal.mods.logistic, env3,matrix2, logistic = T)
plot.prediction("VAR_ST_5_tc", seasonal.mods.logistic, env3,matrix2, logistic = T)
plot.prediction("PRED_SM_5", seasonal.mods.logistic, env3,matrix2, logistic = T)
plot.prediction("PRED_ST_5", seasonal.mods.logistic, env3,matrix2, logistic = T)
dev.off()



all.spp = c("ACHMIL", "ACOROS", "ALLGEY", "ANTMED", 
            "AREFEN", "ARTPAT", "ARTSCO", "BISBIS", "BISVIV", 
            "CALLEP", "CALPUR", "CAMROT", "CAMUNI", "CARALB", "CARELY", 
            "CARILL", "CARNOV", "CARROS", "CARRUP", 
            "CARSCO", "CARSIC", "CASOCC", "CERARV", "CHIJAM", "CRERUN", "DANINT", 
            "DESCAE", "DODPUL", "DRACRA", "ELYSCR", "ELYTRA", 
            "ERIGLA", "ERINAN", "ERIPIN", "ERISIM", "ERYCAP", "FESBRA", 
            "GENALG", "GENAMA", "HELMOR", "KOBMYO", 
            "LEWPYG", "LIGTEN", "LLOSER", "LUZSPI", "MERLAN", 
            "MINOBT", "MINRUB",  "NOCFEN", "OREALP", "PACCAN", "PACCRO", 
            "PARPUL", "PEDGRO", "PEDPAR", "PEDRAC", "PENWHI", 
            "PHLPUL", "POAALP", "POAARC", "POAGLA", "POLVIS", "POTDIV","POTGLA", "POTRUB", "POTSUB", 
            "PRIPAR", "RHOINT", "RHORHO", "SALGLA", "SALIX", "SALPLA", 
            "SAXRHO", "SEDLAN", "SELDEN", "SIBPRO", "SILACA", "SOLMUL", 
            "SOLSIM", "STELON", "STEUMB", "SWEPER", "TETACA", "TETGRA", "TONLYA", 
            "TONPYG", "TRIDAS", "TRIPAR", "TRISPI")

all.spp2 = c("ACHMIL", "ACOROS", "AREFEN", "ARTSCO", "BISBIS", "BISVIV", 
             "CALLEP", "CARRUP", "CARSCO", "CERARV", "DESCAE", "ERISIM", "FESBRA", 
             "GENALG", "KOBMYO", "LEWPYG", "LLOSER", "LUZSPI", "MINOBT", 
             "OREALP", "POTGLA", "Rock", "SEDLAN", "SELDEN", 
             "STELON", "TRIDAS", "TRIPAR", "TRISPI")


ex.dredge = function(species, df){
  holder = as.data.frame(matrix(nrow = 0, ncol = 8))
  colnames(holder) = c("Species","(Intercept)", "MEAN_SM_5", "MEAN_ST_5", "PRED_SM_5", "PRED_ST_5",
                        "VAR_SM_5_tc", "VAR_ST_5_tc")
  
  
  global.mod1 = lmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + year + (1|node), data = test01, na.action = na.fail, REML = FALSE)  
  x0 = dredge(global.mod1, beta = "n" )
  holder[nrow(holder) +1, 2:8] = x0[1, 1:7]
  holder[nrow(holder), "Species"] = species
  return(holder)
  
}


library(PerformanceAnalytics)
chart.Correlation(env3[, 3:8])                            
chart.Correlation(env5)

############################## dredging test ####
ex.dredge("ARTSCO", test01)


all.spp2

test01 = make.df(matrix2df("ARTSCO", matrix), env1)   
gm = lmer(abundance ~ MEAN_SM_5 + MEAN_ST_5 + VAR_SM_5_tc + VAR_ST_5_tc + PRED_SM_5 + PRED_ST_5 + year + (1|node), data = test01, na.action = na.fail, REML = FALSE)  
tt = dredge(gm, beta = "n")

png("/Users/Will/Desktop/ARTSCO_plots.png", height = 9, width = 6, res = 250, units = "in")
par(mfrow = c(3,2), mar = c(4,4, 2,2))
test01 = test01[test01$abundance != 0, ]
plot(test01$abundance ~ test01$MEAN_SM_5)
plot(test01$abundance ~ test01$MEAN_ST_5)
plot(test01$abundance ~ test01$VAR_SM_5_tc)
plot(test01$abundance ~ test01$VAR_ST_5_tc)
plot(test01$abundance ~ test01$PRED_SM_5)
plot(test01$abundance ~ test01$PRED_ST_5)
dev.off()
