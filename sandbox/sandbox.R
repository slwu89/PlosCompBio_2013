rm(list=ls())
dir <- "/Users/slwu89/Desktop/git/PlosCompBio_2013/sandbox/"

library(stats4)
source(paste0(dir,"functions.R"))

################################################################################
# parametersCommon.R
################################################################################

T = 400 # number of time steps to run
Nf = 50 # number of blood-feeding habitats
Nl = 50 # number of larval habitats
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.9
survF = 0.9

sigma = 3 # mosquito incubation period
tau = 3 # host incubation period
rho_max = 40 # rho is host recovery and these lines define its distribution
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 3, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 3, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.55 # mosquito-to-host transmission efficiency

fac = matrix(seq(0.9, 1, length.out = xi), byrow = T, nrow = Nl, ncol = xi)
betaParam1 = 14; betaParam2 = 6
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alphaEven = matrix(betaMean, nrow = Nl, ncol = xi) * fac
alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = xi) * fac
alpha = alphaEven
larvalLocations = sample(1 : Nf, Nl, replace = TRUE)

bites = read.table(paste0(dir,"/bitingDeBenedictis2003.txt"), sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))

hSize = 5.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of hosts across feeding stations, with two-host minimum
Nh = sum(Hh) # total number of hosts


################################################################################
# Set up the landscape
################################################################################

pF = points.poisson(Nf)
XFcoord = pF$x # blood-feeding habitat coordinate X
YFcoord = pF$y # blood-feeding habitat coordinate Y
XLcoord = rnorm(Nl, XFcoord[larvalLocations], .01) # larval habitat coordinate X
YLcoord = rnorm(Nl, YFcoord[larvalLocations], .01) # larval habitat coordinate Y

L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
# just needs to have rows that sum to 1*proportion survive mvmt
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'well')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'well')
# just needs to have row that sum to 1-(proportion of time not at risk)
H = makeHrossMacdonald(Nf, Nh, Hh)
# just needs to have columns that sum to 1
U = makeUrossMacdonald(H)
Nh = dim(H)[1]


################################################################################
# variables.R
################################################################################

# initialize variables

M_larvae = matrix(0, nrow = Nl, ncol = xi)
MM = array(0, dim = c(Nl, xi, T))
# debug(estimateLarvalEquilib)
M_larvae = round(estimateLarvalEquilib())
# undebug(estimateLarvalEquilib)

Sml = round(estimateMosquitoEquilibLarval())
Smf = round(estimateMosquitoEquilibFeeding())
Eml = matrix(0, nrow = Nl, ncol = sigma + 1)
Emf = matrix(0, nrow = Nf, ncol = sigma + 1)
Iml = rep(0, Nl)
Imf = rep(0, Nf)

Sh = rep(0, Nh)
Eh = matrix(0, nrow = Nf, ncol = tau + 1)
Ih = rep(0, Nh)
Rh = rep(0, Nh)

SSml = matrix(0, nrow = Nl, ncol = T)
SSmf = matrix(0, nrow = Nf, ncol = T)
EEml = array(0, dim = c(Nl, sigma + 1, T))
EEmf = array(0, dim = c(Nf, sigma + 1, T))
IIml = matrix(0, nrow = Nl, ncol = T)
IImf = matrix(0, nrow = Nf, ncol = T)

SSh = matrix(0, nrow = Nh, ncol = T)
EEh = array(0, dim = c(Nf, tau + 1, T))
IIh = matrix(0, nrow = Nh, ncol = T)
RRh = matrix(0, nrow = Nh, ncol = T)
