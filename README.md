# Spatial 2D "density" generator

This application is hosted at (hbwaddel.shinyapps.io/GenerateDensity)

For some of my movement modeling work, I wanted a quick way to generate a 2d spatial "density" for the purposes of simulations.

This application generates a function defined over a 2-dimensional square by overlaying multiple bivariate normal distributions on top of each other.

You can generate the density randomly, or fill in the x and y coordinates of the peaks, lists of the horizontal and vertical spreads, and the domain. The "resolution" box determines the resolution of the plot that is generated.

The density also generates the R code that you need to get the density at a point in R. After pasting the R code into your own program, you can generate the density at a point (x,y) with density.func(x,y).
