env01 = read.csv("data/env01.csv")
env02 = read.csv("data/env02.csv")
env03 = read.csv("data/env03.csv")
env04 = read.csv("data/env04.csv")

env01 = env01[env01$plot != "15", ]

env01$NAME = paste(env01$plot, env01$year, sep = "-")
row.names(env01) = env01$NAME
env01_pca = env01[,c(3:20)]

env01_pca = env01_pca[c(1:11, 13:52), ]

env01_pca_sub1 = env01_pca[, c("soiltemp_30cm_avg", "soilmoisture_a_30cm_avg", "soiltemp_30cm_sd", "soilmoisture_a_30cm_sd", "st_acf_30cm_14", "sm_acf_30cm_14")]
env01_pca_sub2 = env01_pca[, c("soiltemp_5cm_avg", "soilmoisture_a_5cm_avg", "soiltemp_5cm_sd", "soilmoisture_a_5cm_sd", "st_acf_5cm_14", "sm_acf_5cm_14")]
env01_pca_sub3 = env01_pca[, c("soiltemp_30cm_avg", "soilmoisture_a_30cm_avg", "soiltemp_30cm_dr", "soilmoisture_a_30cm_dr", "st_acf_30cm_14", "sm_acf_30cm_14")]
env01_pca_sub4 = env01_pca[, c("soiltemp_5cm_avg", "soilmoisture_a_5cm_avg", "soiltemp_5cm_dr", "soilmoisture_a_5cm_dr", "st_acf_5cm_14", "sm_acf_5cm_14")]

library(factoextra)
png("output/env_pca_a.png", height = 5, width = 5, res = 200, units = "in")
pca_sub1 = prcomp(env01_pca_sub1, scale. = T)
fviz_pca_var(pca_sub1, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969")  # Individuals color
dev.off()

png("output/env_pca_b.png", height = 5, width = 5, res = 200, units = "in")
pca_sub2 = prcomp(env01_pca_sub2, scale. = T)
fviz_pca_var(pca_sub2, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969")  # Individuals color
dev.off()

png("output/env_pca_c.png", height = 5, width = 5, res = 200, units = "in")
pca_sub3 = prcomp(env01_pca_sub3, scale. = T)
fviz_pca_var(pca_sub3, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969")  # Individuals color
dev.off()

png("output/env_pca_d.png", height = 5, width = 5, res = 200, units = "in")
pca_sub4 = prcomp(env01_pca_sub4, scale. = T)
fviz_pca_var(pca_sub4, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969")  # Individuals color
dev.off()


env02$NAME = env02$plot
row.names(env02) = env02$NAME
env02_pca = env02[,c(3:19)]
env02_pca = env02_pca[c(1:9, 11:14), ]

env02_pca_sub1 = env02_pca[, c("soiltemp_30cm_avg", "soilmoisture_a_30cm_avg", "soiltemp_30cm_sd", "soilmoisture_a_30cm_sd", "st_acf_30cm_14", "sm_acf_30cm_14")]
env02_pca_sub2 = env02_pca[, c("soiltemp_5cm_avg", "soilmoisture_a_5cm_avg", "soiltemp_5cm_sd", "soilmoisture_a_5cm_sd", "st_acf_5cm_14", "sm_acf_5cm_14")]
env02_pca_sub3 = env02_pca[, c("soiltemp_30cm_avg", "soilmoisture_a_30cm_avg", "soiltemp_30cm_dr", "soilmoisture_a_30cm_dr", "st_acf_30cm_14", "sm_acf_30cm_14")]
env02_pca_sub4 = env02_pca[, c("soiltemp_5cm_avg", "soilmoisture_a_5cm_avg", "soiltemp_5cm_dr", "soilmoisture_a_5cm_dr", "st_acf_5cm_14", "sm_acf_5cm_14")]


png("output/env_pca_ay_a.png", height = 5, width = 5, res = 200, units = "in")
pca_sub1 = prcomp(env02_pca_sub1, scale. = T)
fviz_pca_biplot(pca_sub1, repel = TRUE,
             col.var = "#2E9FDF", # Variables color
             col.ind = "#696969")  # Individuals color
dev.off()

png("output/env_pca_ay_b.png", height = 5, width = 5, res = 200, units = "in")
pca_sub2 = prcomp(env02_pca_sub2, scale. = T)
fviz_pca_biplot(pca_sub2, repel = TRUE,
             col.var = "#2E9FDF", # Variables color
             col.ind = "#696969")  # Individuals color
dev.off()

png("output/env_pca_ay_c.png", height = 5, width = 5, res = 200, units = "in")
pca_sub3 = prcomp(env02_pca_sub3, scale. = T)
fviz_pca_biplot(pca_sub3, repel = TRUE,
             col.var = "#2E9FDF", # Variables color
             col.ind = "#696969")  # Individuals color
dev.off()

png("output/env_pca_d_ay.png", height = 5, width = 5, res = 200, units = "in")
pca_sub4 = prcomp(env02_pca_sub4, scale. = T)
fviz_pca_var(pca_sub4, repel = TRUE,
             col.var = "#2E9FDF", # Variables color
             col.ind = "#696969")  # Individuals color
dev.off()


png("output/var_corr.png", height = 5, width = 10, res = 300, units = "in")
par(mfrow = c(1,2))
plot(env01_pca$soiltemp_30cm_avg ~ env01_pca$soiltemp_30cm_sd, pch = 22, col = "black", bty = "l", xlim = c(0, 1.4),
     ylab = "Daily average soil temperature (C)", xlab = ("Daily average soil temperature variation"))
points(env01_pca$soiltemp_30cm_avg ~ env01_pca$soiltemp_30cm_dr, pch = 23, col = "red")
legend(x=1, y =-1, legend = c("SD", "range"), pch = c(22, 23), cex = 1, col = c("black", "red"), bty = 'n')

plot(env02_pca$soiltemp_30cm_avg ~ env02_pca$soiltemp_30cm_sd, pch = 22, col = "black", bty = "l", xlim = c(0, 1.4),
     ylab = "Daily average soil temperature (C)", xlab = ("Daily average soil temperature variation"))
points(env02_pca$soiltemp_30cm_avg ~ env02_pca$soiltemp_30cm_dr, pch = 23, col = "red")
legend(x=1, y =1, legend = c("SD", "range"), pch = c(22, 23), cex = 1, col = c("black", "red"), bty = 'n')
dev.off()

png("output/pred_variation.png", height = 5, width = 10, res = 300, units = "in")
par(mfrow = c(1,2))
hist(env01_pca$st_acf_5cm_14, breaks = c(0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 
                                            0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 
                                            0.94, 0.96, 0.98, 1), col = rgb(0.7, 0.1, 0.1, 0.3), main = "", xlab = "Predictability (5cm)")
hist(env01_pca$sm_acf_5cm_14, breaks = c(0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 
                                         0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 
                                         0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 
                                         0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1), add =T, col = rgb(0.1, 0.1, 0.7, 0.3))

hist(env01_pca$st_acf_30cm_14, breaks = c(0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 
                                         0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 
                                         0.94, 0.96, 0.98, 1), col = rgb(0.7, 0.1, 0.1, 0.3), main = "", xlab = "Predictability (30cm)")
hist(env01_pca$sm_acf_30cm_14, breaks = c(0.0, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 
                                          0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 
                                          0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 
                                          0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 
                                          0.98, 1), add =T, col = rgb(0.1, 0.1, 0.7, 0.3))
dev.off()


row.names(env03) = paste(env03$node, env03$year, sep = "-")
row.names(env04) = env04$Group.1

env03_sub1 = env03[,c("MEAN_SM_5", "MEAN_ST_5", 
                      "VAR_SM_5_tc", "VAR_ST_5_tc", 
                      "PRED_ST_5", "PRED_SM_5")]
env03_sub2 = env03[,c("MEAN_SM_5", "MEAN_ST_5", 
                      "VAR_SM_5_cv", "VAR_ST_5_cv", 
                      "PRED_ST_5", "PRED_SM_5")]


env03_sub3 = env03[,c("MEAN_SM_30", "MEAN_ST_30", 
                      "VAR_SM_30_tc", "VAR_ST_30_tc", 
                      "PRED_ST_30", "PRED_SM_30")]
env03_sub4 = env03[,c("MEAN_SM_30", "MEAN_ST_30", 
                      "VAR_SM_30_cv", "VAR_ST_30_cv", 
                      "PRED_ST_30", "PRED_SM_30")]


pca_sub1 = prcomp(env03_sub1, scale. = T)
fviz_pca_var(pca_sub1, repel = TRUE,
             col.var = "#2E9FDF", # Variables color
             col.ind = "#696969")

pca_sub2 = prcomp(env03_sub2, scale. = T)
fviz_pca_var(pca_sub2, repel = TRUE,
             col.var = "#2E9FDF", # Variables color
             col.ind = "#696969")

par(mfrow = c(1,2))
pca_sub3 = prcomp(env03_sub3, scale. = T)
fviz_pca_var(pca_sub3, repel = TRUE,
             col.var = "#2E9FDF", # Variables color
             col.ind = "#696969")

pca_sub4 = prcomp(env03_sub4, scale. = T)
fviz_pca_var(pca_sub4, repel = TRUE,
             col.var = "#2E9FDF", # Variables color
             col.ind = "#696969")



env04_sub1 = env04[,c("MEAN_SM_5", "MEAN_ST_5", 
                      "VAR_SM_5_tc", "VAR_ST_5_tc", 
                      "PRED_ST_5", "PRED_SM_5")]
env04_sub2 = env04[,c("MEAN_SM_5", "MEAN_ST_5", 
                      "VAR_SM_5_cv", "VAR_ST_5_cv", 
                      "PRED_ST_5", "PRED_SM_5")]

env04_sub3 = env04[,c("MEAN_SM_30", "MEAN_ST_30", 
                      "VAR_SM_30_tc", "VAR_ST_30_tc", 
                      "PRED_ST_30", "PRED_SM_30")]
env04_sub4 = env04[,c("MEAN_SM_30", "MEAN_ST_30", 
                      "VAR_SM_30_cv", "VAR_ST_30_cv", 
                      "PRED_ST_30", "PRED_SM_30")]



pca_sub5 = prcomp(env04_sub1, scale. = T)
fviz_pca_biplot(pca_sub5, repel = TRUE,
             col.var = "red", # Variables color
             col.ind = "#696969")

pca_sub6 = prcomp(env04_sub2, scale. = T)
fviz_pca_biplot(pca_sub6, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969")

pca_sub7 = prcomp(env04_sub3, scale. = T)
fviz_pca_biplot(pca_sub7, repel = TRUE,
                col.var = "red", # Variables color
                col.ind = "#696969")

pca_sub8 = prcomp(env04_sub4, scale. = T)
fviz_pca_biplot(pca_sub8, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969")


pairs(env04)
