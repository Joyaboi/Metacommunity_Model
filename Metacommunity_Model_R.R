#### Install MCSim ####

install.packages("devtools")
install.packages("tidyverse")
install.packages("fields")
install.packages("vegan")
install.packages("betapart")
devtools::install_github('sokole/MCSim')

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