rm(list=ls(all.names=TRUE))

#### Install MCSim ####

# install.packages("devtools")
# install.packages("tidyverse")
# install.packages("fields")
# install.packages("vegan")
# install.packages("betapart")
# devtools::install_github('sokole/MCSim')

library(MCSim)
library(tidyverse)
library(fields)
library(vegan)
library(betapart)

####Model Initialization####

#Set a seed
set.seed(1234)

#We will simulate 10 species over 50 discrete sites. We will first define the sites and their
#characteristics. 
nsites <- 50 ##number of sites

# Let's generate a data frame with coordinates for nsites sites,
#randomly distributed in a 10 x 10 environment
xy.coordinates <- data.frame(
  x = runif(nsites)*10,
  y = runif(nsites)*10
)
xy.coordinates <-
  xy.coordinates[order(xy.coordinates$x,xy.coordinates$y),]
print(xy.coordinates)
plot(y ~ x, xy.coordinates)

##The following creates a list with information on sites that is
#later needed by the function running the model
my.landscape <- MCSim::fn.make.landscape(
  site.coords = xy.coordinates, ## site coordinates
  m = rep(0.05,nsites),#rep(0.01,10), ##immigration rate for each site from the species pool
  Ef = sort(runif(nsites)), ## Environmental variable value for each
  #site - since the sites are ordered on the x-y axise, using sort()
  #introduces an environmental gradient along the x axis – remove sort
  #for a more random environment
  JM = 1000000 ## metacommunity size (sensu Hubbell 2001) – i.e. species pool size
)

print(my.landscape$site.info)

#We can plot the sites again, colouring them according to their environment:
  plot(y ~ x, my.landscape$site.coords, pch=16,
       col=tim.colors(n=64)[ceiling(my.landscape$site.info$Ef*64)])
  
#We will then define the 10 species, including their niches along the environmental axis. Each
#niche is modelled as a Normal distribution, which is typical for Hutchinsonian niches, as we
#saw in class.
# niche positions, niche breadths, and relative abundances for 10 species
niche.positions <- runif(10) ##niche optima are distributed randomly along the environmental axis
niche.breadths <- rep(0.1,10) ##all species have the same niche width of 0.1 (keep in mind the environment varies between [0,1])
regional.rel.abund <- rep(1/10,10) ##each species will be
#initialised with the same number of individual in each site

#####Plot niches#####
# -- function for plotting bell curves (the curve of a normal distribution)
fn_norm_curve <- function(sigma=1, mu = 0,...) {
  curve(
    (1/sigma * sqrt(2 * pi)) * exp((-1 *(x - mu)^2) / (2 *
                                                         sigma^2)), ...) #formula for bell curve
}
# -- Initialize plot
plot(1,1,
     xlim = c(-0.2,1.2),
     ylim = c(0, (1/niche.breadths[1]* sqrt(2 * pi))),
     type = 'n',
     xlab = 'Environmental gradient',
     ylab = 'Prob. dens.',
     main = 'Niche positions and niche widths')
mypal <- tim.colors(length(niche.positions))
# -- loop to plot each species' habitat preference
for (i.spp in 1:length(niche.positions)){
  fn_norm_curve(
    mu = niche.positions[i.spp],
    sigma = niche.breadths[i.spp],
    add = TRUE,
    col = mypal[i.spp],
    n=1000)
}
# -- plot sites along the x-axis
rug(my.landscape$site.info$Ef)
legend(x = 'topleft',
       legend = 1:10,
       lty = 1,
       col = mypal,
       cex = .75,
       ncol = 4)

####Run The Model####

#We will run simulations over 100 timesteps. If you play with the model (i.e. if you increase
#simulation length and look at the different community patterns), you will see this is not
#always enough for the model to reach equilibrium depending on the parameters, but this is
#enough for this practical while keeping simulation times reasonable.
class(niche.breadths)
n.timestep <- 100
sim.result <- MCSim::fn.metaSIM(
  landscape = my.landscape,
  trait.Ef = niche.positions,
  trait.Ef.sd = niche.breadths,
  gamma.abund = regional.rel.abund,
  W.r = 2, ##slope of exponential dispersal kernel – the higher
  #the value, the more restricted the dispersal. 0 means infinite dispersal
  nu = 0, ##probability that a novel species will appear during a recruitment event
  n.timestep = n.timestep,
  sim.ID = "my_test_sim",
  output.dir.path = "my_sim_output_directory")
head(sim.result$J.long)

##### Analyse the model outputs

#If you look at sim.result$J.long, which is the main output of the function, it contains
#the species composition of all sites over all the simulated timesteps. You therefore need to
#extract the part of the data frame corresponding to the last time step. Once this is done, you
#need to transform this data frame into a species-by-site data frame, as usual. This can be
#done with the following code:

MC.fin <-
  sim.result$J.long[which(sim.result$J.long$timestep==n.timestep),]
MC.fin.mat <- MC.fin %>%
  pivot_wider(names_from=spp,values_from=c(count)) ##this data frame contains abundance data
MC.fin.mat <- MC.fin.mat[,-c(1,2)]
MC.fin.mat.pa <- as.data.frame((MC.fin.mat>0)*1) ##we convert the abundance data frame into a presence-absence one to compute some indices afterwards

#Using these data frames, compute gamma diversity, alpha diversity (richness-based), the
#Shannon index, incidence-based beta diversity from the presence-absence data (using
#Sørensen and Simpson indices), and the abundance-based beta diversity using the Bray-
#Curtis index (whis is the equivalent of Sørensen for abundance data):

  Gamma <- length(which(colSums(MC.fin.mat)>0))
#Alpha
#rowSums(MC.fin.mat.pa)
Alpha <- mean(rowSums(MC.fin.mat.pa))
#Shannon
#diversity(MC.fin.mat)
Shannon <- mean(diversity(MC.fin.mat))
#Beta
#beta.pair(MC.fin.mat.pa)
#beta.pair.abund(MC.fin.mat)$beta.bray
beta <- lapply(beta.pair(MC.fin.mat.pa),mean)
beta.Sim <- beta$beta.sim
beta.Sor <- beta$beta.sor
beta.ab <- mean(beta.pair.abund(MC.fin.mat)$beta.bray)

####Check for Convergence####

#Run the following code at the end of the previous one. What does this code do, and can you
#consider that the simulation has reached equilibrium?

MC.fin.mat.pa <- array(NA,c(n.timestep,50,10))
MC.fin.mat <- array(NA,c(n.timestep,50,10))
for(i in 1:n.timestep){
  MC.fin <- sim.result$J.long[which(sim.result$J.long$timestep==i),]
  MC.fin.mat.temp <- MC.fin %>%
    pivot_wider(names_from=spp,values_from=c(count)) ##this data frame contains abundance data
  MC.fin.mat[i,,] <- as.matrix(MC.fin.mat.temp[,-c(1,2)])
  MC.fin.mat.pa[i,,] <- as.matrix(as.data.frame((MC.fin.mat.temp[,-
                                                                   c(1,2)]>0)*1)) ##we convert the abundance data frame into a presence-absence one to compute some indices afterwards
}
oc.time <- apply(MC.fin.mat.pa,c(1),colSums)
plot(1:n.timestep,oc.time[1,],type="l",col=tim.colors(n=10)[1],ylim=
       c(0,max(oc.time)))
for(i in 2:10){
  lines(1:n.timestep,oc.time[i,],col=tim.colors(n=10)[i])
}
ab.time <- apply(MC.fin.mat,c(1),colSums)
plot(1:n.timestep,ab.time[1,],type="l",col=tim.colors(n=10)[1],ylim=
       c(0,max(ab.time)))
for(i in 2:10){
  lines(1:n.timestep,ab.time[i,],col=tim.colors(n=10)[i])
}

#### Run Multiple Simulations for Different Paramater Values ####

#We want to explore the influence of dispersal ability and niche width on the model outputs
#and the community patterns described above. You will then run your script within two for
#loops, to increment these two parameters. These for loop will go through the following two
#vectors that you need to initialise

niche.width <- c(0.01,0.1,2,5,1000) ##vector of niche widths
disp.ab <- c(0,2,5,10,1000,5e6) ##vector of dispersal abilities – as mentioned above, the higher the value, the more restricted the dispersal. 0 means infinite dispersal

#Parameters niche.breadths and W.r will take these values in the code above.

#Store the values of gamma diversity, alpha diversity, Shannon index and the different beta
#diversity in 5 x 6 matrices. Store also MC.fin.mat in a list of list, as it may help you to better
#understand how the community patterns change.

Gamma <- Alpha <- Shannon <- beta.Sim <- beta.Sor <- beta.ab <-
  matrix(NA,length(niche.width),length(disp.ab))
MC.fin.mat.store <- list()

#You now just need to add the two for() loops to your code, and the indices to the patterns,
#and run the model. It should take about 7 minutes, or you can load the
#Simulation_results.RData file instead.

# Number of time steps
n.timestep <- 100

# Outer loop: niche width
for (i in seq_along(niche.width)) {
  
  # Inner loop: dispersal ability
  for (j in seq_along(disp.ab)) {
    
    cat("Running sim with niche width =", niche.width[i],
        "and dispersal =", disp.ab[j], "\n")
    
    # Define species traits for this iteration
    niche.positions <- runif(10)
    niche.breadths <- rep(niche.width[i], 10)
    regional.rel.abund <- matrix(rep(1/10, 10), ncol = 1)
    
    # Run the simulation
    sim.result <- try(
      MCSim::fn.metaSIM(
        landscape = my.landscape,
        trait.Ef = niche.positions,
        trait.Ef.sd = niche.breadths,
        gamma.abund = regional.rel.abund,
        W.r = disp.ab[j],
        nu = 0,
        n.timestep = n.timestep,
        sim.ID = paste0("sim_", i, "_", j),
        output.dir.path = "my_sim_output_directory"
      ),
      silent = TRUE
    )
    
    # Skip to next iteration if error occurs
    if (inherits(sim.result, "try-error")) {
      cat("Simulation failed for i =", i, "j =", j, "\n")
      next
    }
    
    # Extract the last timestep community matrix
    MC.fin <- sim.result$J.long[which(sim.result$J.long$timestep == n.timestep), ]
    MC.fin.mat <- MC.fin %>%
      pivot_wider(names_from = spp, values_from = count)
    MC.fin.mat <- MC.fin.mat[, -c(1, 2)]
    MC.fin.mat.pa <- as.data.frame((MC.fin.mat > 0) * 1)
    
    # Compute metrics
    Gamma[i, j] <- length(which(colSums(MC.fin.mat) > 0))
    Alpha[i, j] <- mean(rowSums(MC.fin.mat.pa))
    Shannon[i, j] <- mean(vegan::diversity(MC.fin.mat))
    
    beta_vals <- lapply(betapart::beta.pair(MC.fin.mat.pa), mean)
    beta.Sim[i, j] <- beta_vals$beta.sim
    beta.Sor[i, j] <- beta_vals$beta.sor
    beta.ab[i, j] <- mean(betapart::beta.pair.abund(MC.fin.mat)$beta.bray)
    
    # Store final community for later inspection
    MC.fin.mat.store[[paste(i, j, sep = "_")]] <- MC.fin.mat
  }
}

cat("All simulations completed.\n")

#####Visualize Matrices as Heatmaps#####

##Plot alpha diversity
image(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(Alpha),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  col = heat.colors(20)
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Alpha diversity")

# Add color bar
image.plot(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(Alpha),
  col = heat.colors(20),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  legend.lab = "Value"
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Alpha diversity")

##Plot Shannon diversity
image(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(Shannon),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  col = heat.colors(20)
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Shannon diversity")

# Add color bar
image.plot(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(Shannon),
  col = heat.colors(20),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  legend.lab = "Value"
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Shannon diversity")

##Plot Sorensen diversity
image(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(beta.Sor),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  col = heat.colors(20)
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Sorensen diversity")

# Add color bar
image.plot(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(beta.Sor),
  col = heat.colors(20),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  legend.lab = "Value"
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Sorensen diversity")

##Plot Simpson diversity
image(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(beta.Sim),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  col = heat.colors(20)
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Simpson diversity")

# Add color bar
image.plot(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(beta.Sim),
  col = heat.colors(20),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  legend.lab = "Value"
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Simpson diversity")

##Plot Bray-Curtis diversity
image(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(beta.ab),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  col = heat.colors(20)
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Bray-Curtis diversity")

# Add color bar
image.plot(
  x = 1:length(disp.ab),
  y = 1:length(niche.width),
  z = t(beta.ab),
  col = heat.colors(20),
  axes = FALSE,
  xlab = "Dispersal ability",
  ylab = "Niche width",
  legend.lab = "Value"
)
axis(1, at = 1:length(disp.ab), labels = disp.ab)
axis(2, at = 1:length(niche.width), labels = rev(niche.width))
title("Bray-Curtis diversity")