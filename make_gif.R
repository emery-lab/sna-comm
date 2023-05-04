library(magick)

data  = read.csv("data/sn02.csv")
colnames(data)
data$date = as.Date(data$date)

p = ggplot(data) +
  geom_point (aes(x=node, y=soiltemp_5cm_avg, colour = soiltemp_5cm_avg, size = 4)) +
  scale_colour_gradient(low="blue", high="red")+ 
  theme_bw() +
  transition_time(date) +
  ease_aes('linear')

animate(p, duration = 40)

rbPal <- colorRampPalette(c('brown','blue'))
data$Col <- rbPal(10)[as.numeric(cut(data$soilmoisture_a_5cm_avg,breaks = 10))]
dates[730]

dates = unique(data$date)
for(i in 366:730){
  
  png(paste("gif/", dates[i], ".gif", sep = ""), height = 5, width = 5, res = 180, units = "in")
  par(mar = c(4,4,2,2))
  plot(x=data[data$date %in% dates[i], "node"],
       y=data[data$date %in% dates[i], "soilmoisture_a_5cm_avg"], 
       ylim = c(-0.1, 1), pch = 19, cex = 2, 
       col = data[data$date %in% dates[i], "Col"],
       xlab = "Sensor Node", ylab = "Average Daily Soil Temperature", bty = "l",
       main = dates[i])
  abline(h = 0.13, lty =2 , col = "gray")
  for(k in 1:10){
    points(data[data$date %in% dates[i-k], "node"], data[data$date %in% dates[i-k], "soilmoisture_a_5cm_avg"], pch = 19, cex = 2, col = alpha(data[data$date %in% dates[i-k], "Col"], 1-(0.1*k)))
  }
  dev.off()
  
}

imgs <- list.files("gif", full.names = TRUE)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 20)
image_write(image = img_animated,
            path = "test2.gif")

dates = seq.Date(as.Date("1-1-2019", format = "%m-%d-%Y"), 
                 as.Date("12-31-2019", format = "%m-%d-%Y"), 
                 length.out = 365)

par(mar = c(0,0,0,0))

for(i in 1:365){
  png(paste("gif/", dates[i], ".gif", sep = ""), height = 4, width = 6, res = 300, units = "in")
  
  plot.new()
  plot.window(xlim = c(0,1), ylim = c(0,1))
  
  polygon(c(0,1,1,0), c(0.45, 0.45,0.55, 0.55))
  polygon(c(0,i/365, i/365, 0), c(0.45, 0.45, 0.55, 0.55), col = "gray")
  
  #segments(0,0.5, 1,0.5, lwd = 3)
  #segments(0,0.45, 0, 0.55, lwd = 3)
  #segments(1,0.45, 1, 0.55, lwd = 3)
  text(x=0.5, y = 0.65, dates[i], cex = 1.7)
  
  #segments(x0 = i/365, y0 = 0.45, x1 = i/365, y1 = 0.55, lwd = 4)
  dev.off()
  
}

imgs <- list.files("gif", full.names = TRUE)
img_list <- lapply(imgs, image_read)
img_joined <- image_join(img_list)
img_animated <- image_animate(img_joined, fps = 20)
image_write(image = img_animated,
            path = "timeline.gif")
