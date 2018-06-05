rm(list=ls())
dir <- "/Users/slwu89/Desktop/git/PlosCompBio_2013/sandbox/"

library(stats4)
source(paste0(dir,"functions.R"))

################################################################################
# parametersCommon.R
################################################################################

T = 600 # number of time steps to run
Nf = 50 # number of blood-feeding habitats
Nl = 50 # number of larval habitats
v = 50 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.95
survF = 0.95

sigma = 3 # mosquito incubation period
tau = 3 # host incubation period
rho_max = 40 # rho is host recovery and these lines define its distribution
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 3, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 3, .2)))
c = 0.5 # host-to-mosquito transmission efficiency
b = 0.5 # mosquito-to-host transmission efficiency

fac = matrix(seq(0.95, 1, length.out = xi), byrow = T, nrow = Nl, ncol = xi)
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
M_larvae = round(estimateLarvalEquilib())

Sml = round(estimateMosquitoEquilibLarval())
Smf = round(estimateMosquitoEquilibFeeding())
Eml = matrix(0, nrow = Nl, ncol = sigma)
Emf = matrix(0, nrow = Nf, ncol = sigma)
Iml = rep(0, Nl)
Imf = rep(0, Nf)

SSml = matrix(0, nrow = Nl, ncol = T)
SSmf = matrix(0, nrow = Nf, ncol = T)
EEml = array(0, dim = c(Nl, sigma + 1, T))
EEmf = array(0, dim = c(Nf, sigma + 1, T))
IIml = matrix(0, nrow = Nl, ncol = T)
IImf = matrix(0, nrow = Nf, ncol = T)

# humans
Sh = rep(1, Nh)
Eh = matrix(0, nrow = Nh, ncol = tau)
Ih = matrix(0, nrow = Nh, ncol = rho_max)
Rh = rep(0, Nh)

SSh = matrix(0, nrow = Nh, ncol = T)
EEh = array(0, dim = c(Nf, tau + 1, T))
IIh = matrix(0, nrow = Nh, ncol = T)
RRh = matrix(0, nrow = Nh, ncol = T)


################################################################################
# simulateEpidemic.R
################################################################################

hh = 1

Sh[hh] = 0
Ih[hh, 1] = 1
EIR_cum = rep(0, T + 1)
EIR_cum[1] = Nh - 1
Ih_sum = rep(0, T + 1)
Ih_sum[1] = 1

# loop over time
# for(tt in 2 : (T + 1)){
for(tt in 2 : 10){

  # mosquito movement from blood feeding habitats to larval habitats
  # Sml: susceptible mosquitoes at egg laying sites
  # Sml = Emerging Larvae + Incoming Migration of susceptible mosquitoes from Feeding Sites
  Sml = M_larvae[,xi] + rmultinomMatrix(Smf, L)

  # Eml: exposed mosquitoes at egg laying sites
  # Eml = Incoming migration of exposed mosquitoes from feeding sites for all stages of incubation (1,...,sigma)
  for(ii in 1 : sigma){
    Eml[, ii] = matrix(rmultinomMatrix(Emf[, ii], L), nrow = Nl, ncol = 1)
  }

  # Iml: infectious mosquitoes at egg laying sites
  # Iml = Incoming migration of infectious mosquitoes from feeding sites
  Iml = rmultinomMatrix(Imf, L)

  # egg laying and advancement through larval stages

  # M_larvae: incoming L1's in site i are Poisson distributed according to v * Sml_{i}
  # (this replaces that column of the larvae matrix, this is okay because we already accounted for the emerging ones when we set Sml above)
  M_larvae = cbind(rpois(Nl, v * (Sml + rowSums(Eml) + Iml)), M_larvae[, c(1 : xi)[-xi]])
  # binomial density dependent survival (apply it to the entire matrix)
  M_larvae = matrix(rbinom(Nl * xi, M_larvae, (M_larvae + 1) ^ (matrix(rep(alpha, xi), Nl, xi) - 1)), Nl, xi)

  # mosquito movement from larval habitats to blood feeding habitats
  # Smf: susceptible mosquitoes at feeding sites
  # Smf = Incoming migration of susceptible mosquitoes from egg laying sites
  Smf = rmultinomMatrix(Sml, F)

  # Emf: Exposed mosquitoes at feeding sites
  # Emf = migration from exposed mosquitoes at egg laying sites in incubation stages (1,...,sigma-1)
  #       going into exposed mosquitoes at feeding sites in incubation stages (2,...,sigma)
  for(ii in 1 : (sigma - 1))
    Emf[, ii + 1] = matrix(rmultinomMatrix(Eml[, ii], F), nrow = Nf, ncol = 1)

  # Imf: infectious mosquitoes at feeding sites
  # Imf = Incoming migration of infectious mosquitoes from egg laying sites +
  #       incoming migration/advancement of exposed mosquitoes at egg laying sites
  Imf = rmultinomMatrix(Iml, F) + rmultinomMatrix(Eml[, sigma], F)

  # transmission from infectious hosts to susceptible mosquitoes
  # This only occurs at feeding sites (humans do not go to egg laying sites)
  # this is Smf (susceptible mosquitoes at feeding sites) biting on infected hosts and geting infected
  sapply(1 : Nf, function(ii){
    sum(rmultinom(1, Smf[ii], U[, ii])[which(rowSums(Ih) > 0)])
  })

  if(sum(Ih) > 1){
    # for each feeding site (i) do:
    #   1.  partition the bites from susceptible mosquitoes in site i onto all people site i following
    #       multinomial distribution with probability given by U[,i], then sum how many bites fell onto infectious humans (that's the output)
    #   2.  the number of mosquitoes getting infected (S -> E)  is binomial(n=that number of bites,p = c (transmission efficiency))
    Emf[, 1] = rbinom(Nf, sapply(1 : Nf, function(ii) sum(rmultinom(1, Smf[ii], U[, ii])[which(rowSums(Ih) > 0)])), c)
  }
  if(sum(Ih) == 1){
    Emf[, 1] = rbinom(Nf, sapply(1 : Nf, function(ii) sum(rmultinom(1, Smf[ii], U[, ii])[which(rowSums(Ih) > 0)])), c)
  }
  if(sum(Ih) == 0){
    Emf[, 1] = 0
  }
  Smf = Smf - Emf[, 1]

  # infectious host recovery
  IhOld = Ih
  for(ii in seq(rho_max - 1, 1, -1))
    Ih[, ii + 1] = rbinom(Nh, Ih[, ii], 1 - rho_fail[ii])
  Rh = Rh + rowSums(IhOld) - rowSums(Ih[, 2 : rho_max])
  rm(IhOld)

  # host progression through the pathogen incubation period
  Ih[, 1] = Eh[, tau]
  for(ii in seq(tau - 1, 1, -1))
    Eh[, ii + 1] = Eh[, ii]

  # transmission from infectious mosquitoes to susceptible hosts
  secBites_h = rowSums(sapply(1 : Nf, function(ii) rmultinom(1, Imf[ii], U[, ii])))
  Eh[, 1] = rbinom(Nh, Sh, 1 - (1 - b) ^ secBites_h)
  Sh = Sh - Eh[, 1]

  EIR_cum[tt] = Nh - sum(Sh)
  Ih_sum[tt] = length(which(rowSums(Ih) > 0))
}
