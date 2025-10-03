################################################################################
# Project Colors
################################################################################
#
# Lina Aragón
# linamaragonb@gmail.com
# September 24, 2025
#
# Description
#
################################################################################

#Colors for species
library(viridis)
viridis::viridis(6,option = "C", direction = -1, begin = 0.1, end = 0.9)
# "#FCCE25FF" "#F79143FF" "#DD5E66FF" "#B6308BFF" "#8204A7FF" "#42049EFF"
# Create a blank plot with colored rectangles
my_colors <- c("#FCCE25FF", "#F79143FF", "#B6308BFF", "#8204A7FF", "#42049EFF")
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 5), axes = FALSE, xlab = "", ylab = "")
for(i in 1:5) {
  rect(0, i-1, 1, i, col = my_colors[i], border = NA)
  text(0.5, i-0.5, my_colors[i], col = "white", cex=1.5)
}


#Colors for biologically hierarchical levels
viridis::viridis(10, option = "F", direction = -1, begin = 0.1, end = 0.9)
[1] "#F7C5A5FF" "#ED4F3EFF" "#931C5BFF" "#261433FF"
[1] "#F7C5A5FF" "#F4835BFF" "#E33641FF" "#AA185AFF" "#641F54FF" "#261433FF"
[1] "#F7C5A5FF" "#F5966DFF" "#F16244FF" "#DD2C45FF" "#B41658FF" "#821E5AFF" "#521E4DFF" "#261433FF"
[1] "#F7C5A5FF" "#F6A178FF" "#F47A54FF" "#ED4F3EFF" "#D92847FF" "#B91657FF" "#931C5BFF" "#6D1F56FF" "#481C48FF" "#261433FF"

my_colors <- c("#F7C5A5FF", "#ED4F3EFF", "#931C5BFF", "#261433FF")
par(mar = c(0, 0, 0, 0))
plot(1, type = "n", xlim = c(0, 1), ylim = c(0, 4), axes = FALSE, xlab = "", ylab = "")
for(i in 1:5) {
  rect(0, i-1, 1, i, col = my_colors[i], border = NA)
  text(0.5, i-0.5, my_colors[i], col = "white", cex=1.5)
}


