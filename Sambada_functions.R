# Function to add a scalebar to a base-graphics plot
Scalebar = function(bar.text, xpos, ypos, lbar) { # text legend, xpos: position on x axis (in % of x axis), ypos: same for y axis, lbar = length of scale
  
  # Get plot coordinates
  pc = par("usr") 
  
  y.bar = rep(pc[3]+((pc[4]-pc[3])*ypos),2)
  x.bar = c(pc[1]+((pc[2]-pc[1])*xpos), pc[1]+((pc[2]-pc[1])*xpos)+lbar  )
  
  lines(x.bar, y.bar)
  text(mean(x.bar), mean(y.bar), labels=bar.text, pos=1)
}
