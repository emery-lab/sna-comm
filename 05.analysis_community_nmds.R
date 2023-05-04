set.seed(11)

library("vegan")
## Vegetation data 
data = read.csv("data/veg01.csv")

table(data$year)
str(data)

data$year = as.character(data$year)

years = unique(data$year)
years = as.character(years)

data = data[data$year != 2017, ]

table(data$multi.hit_or_top.hit)
data$multi.hit_or_top.hit = ifelse(data$multi.hit_or_top.hit == "top-hit-of-multihit-protocol", 
                                   "top-hit", data$multi.hit_or_top.hit)
data$multi.hit_or_top.hit = ifelse(data$multi.hit_or_top.hit == "non-top-hit-of-multihit-protocol", 
                                   "multi-hit", data$multi.hit_or_top.hit)

data = data[data$multi.hit_or_top.hit != "present", ]


xyear = list()
xyear_p = list()
matrix.xyear = list()
matrix.xyear_p = list()

j = 1

str(data)

for(i in years){
  data.xyear = data[data$year == i, ]
  data.xyear = data.xyear[data.xyear$multi.hit_or_top.hit == "top-hit", ]
  data.xyear = data.xyear[, c("plot_name", "abundance", "species")]
  
  data.xyear_p = data.xyear
  data.xyear_p$abundance = 1
  
  xyear[[j]] = data.xyear
  xyear_p[[j]] = data.xyear_p
  matrix.xyear[[j]] = picante::sample2matrix(data.xyear)
  matrix.xyear_p[[j]] = picante::sample2matrix(data.xyear_p)
  
  j = j + 1
}

for(i in 1:4){
  matrix.xyear[[i]] = matrix.xyear[[i]][, -which(names(matrix.xyear[[i]]) %in% c("Rock", "Litter", "Soil", "Moss"))]
  matrix.xyear_p[[i]] = matrix.xyear_p[[i]][, -which(names(matrix.xyear_p[[i]]) %in% c("Rock", "Litter", "Soil", "Moss"))]
  
  matrix.xyear[[i]] = decostand(matrix.xyear[[i]], method = "hellinger")
  
}


nmds = list()
nmds_p = list()

for(i in 1:4){
  nmds[[i]] = metaMDS(matrix.xyear[[i]], k = 2)
  nmds_p[[i]] = metaMDS(matrix.xyear_p[[i]], k = 2)
}


ordiplot(nmds[[1]], type = "n")
orditorp(nmds[[1]], display = "sites")


data$name = paste(data$sensor_node, sep = "-")
data.sub = data[, c("name", "abundance", "species")]
matrix = picante::sample2matrix(data.sub)
matrix = matrix[, -which(names(matrix) %in% c("Rock", "Litter", "Soil", "Moss", "Unknown_seedling"))]

full.nmds = metaMDS(matrix, distance = "bray", k = 2)


sites = full.nmds$points
species = full.nmds$species
par(mfrow = c(1,1))
plot(sites[,1], sites[,2], xlim = c(-2,2), ylim = c(-2,2), pch = 21, col = "black", bg = cols, cex = 1.2, bty = "n", xlab = "MDS1", ylab = "MDS2")
abline(v = 0, lty = 2, col = "gray")
abline(h = 0, lty = 2, col = "gray")
text(sites[,1] + 0.10, sites[,2], labels = row.names(sites))

c = sqrt(species[,1]^2 + species[,2]^2)

species2 = cbind(species, c)
species2 = species2[species2[,3] > 1.25, ]
species2 = species2[c(2:6,13), ]
segments(0,0, species2[,1], species2[,2], col = "blue", lty = 2)
text(species2[,1], species2[,2], labels = rownames(species2))

matrix[, rownames(species2)]

matrix[14:17, ]
matrix[c(3,6,11), ]
matrix[c(5,8,9,12), ]
matrix[c(2,4), ]

cols = c("black", "green", "blue", "green", "orange", "blue", "black", "orange", "orange", "black", "blue", "orange", "black", "red", "red", "red", "red")

