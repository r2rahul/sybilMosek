context("Test Optimization Methods and Associated Functions")
library(sybilMosek)
test_lp1 <- function(){
  lp <- optObj(solver = "sybilMosek", method = "mosek")
  lp <- initProb(lp)
  cm <- Matrix(c(0.5, 2, 1, 1), nrow = 2)
  loadLPprob(lp, nCols = 2, nRows = 2, mat = cm,
   lb = c(0, 0), ub = rep(1000, 2), obj = c(1, 1),
   rlb = c(0, 0), rub = c(4.5, 9), rtype = c("U", "U"),
   lpdir = "max")
  status <- solveLp(lp)
  return(list(status = status, lp = lp))
}

test_that("lp optimization works", {
  expect_equal(test_lp1()$status, 0)
  expect_equal(getObjVal(test_lp1()$lp), 6)
  expect_equal(getFluxDist(test_lp1()$lp), c(3, 3))
  expect_equal(getRedCosts(test_lp1()$lp), c(-8.243784e-12, -8.243784e-12))
  expect_equal(getSolStat(test_lp1()$lp), 0)
  expect_equal(getNumRows(test_lp1()$lp), 2)
  expect_equal(getNumCols(test_lp1()$lp), 2)
  expect_equal(getColsLowBnds(test_lp1()$lp), c(0, 0))
  expect_equal(getColsUppBnds(test_lp1()$lp), c(1000, 1000))
  expect_equal(getObjCoefs(test_lp1()$lp), c(1, 1))
})

test_mip1 <- function(){
  lp <- optObj(solver = "sybilMosek", method = "mosek", pType = "mip")
  lp <- initProb(lp)
  cm <- Matrix(rbind(c(50, 31),
    c( 3, -2)), sparse=TRUE)
  loadLPprob(lp, nCols = 2, nRows = 2, mat = cm,
   lb = c(  0,   0), ub = c(Inf, Inf), obj =  c(1, 0.64),
   rlb = c(-Inf,  -4), rub = c( 250, Inf), ctype = c("I", "I"),
   rtype = c("U", "L"),
   lpdir = "max")
  status <- solveLp(lp)
  return(list(status = status, lp = lp))
}
test_that("MIP Optimization Works", {
  expect_equal(getObjVal(test_mip1()$lp), 5)
  expect_equal(getFluxDist(test_mip1()$lp), c(5.000000e+00, 1.015061e-15))
})


test_mip2 <- function(){
  lp <- optObj(solver = "sybilMosek", method = "mosek", pType = "mip")
  lp <- initProb(lp)
  cm <- Matrix(c(1, 10), nrow=1, sparse=TRUE)
  loadLPprob(lp, nCols = 2, nRows = 1, mat = cm,
             lb = rep(0,0), ub = c(10, 10), obj =  c(1, 1),
             rlb = c(20),
             rtype = c("U"),
             lpdir = "max", ctype = c("B", "B"))
  status <- solveLp(lp)
  return(list(status = status, lp = lp))
}
test_that("MIP With Binary Variable Optimization Works", {
  expect_equal(getObjVal(test_mip2()$lp), 2)
  expect_equal(getFluxDist(test_mip2()$lp), c(1, 1))
})

test_qp1 <- function(){
  lp <- optObj(solver = "sybilMosek", method = "mosek", pType = "qp")
  lp <- initProb(lp)
  cm <- Matrix(c(1, 1, 1), nrow=1, sparse=TRUE)
  loadLPprob(lp, nCols = 3, nRows = 1, mat = cm,
             lb = rep(0,3), ub = rep(Inf,3), obj =  c(0, -1, 0),
             rlb = c(1), rub = c(Inf),
             rtype = c("L"),
             lpdir = "min")
  # Specify the quadratic objective matrix in triplet form.
  i <- c(1,  3,   2,  3)
  j <- c(1,  1,   2,  3)
  v <- c(2, -1, 0.2,  2)
  mat <- Matrix(rep(0, 9), ncol = 3, sparse = TRUE)
  for (k in seq(along = i)){
    mat[i[k], j[k]] <- v[k]
  }
  loadQobj(lp, mat)
  status <- solveLp(lp)
  return(list(status = status, lp = lp))
}
test_that("Quadratic Optimization Works", {
  expect_equal(getObjVal(test_qp1()$lp), -2.5)
  expect_equal(getFluxDist(test_qp1()$lp), c(0.0001586161, 4.9999999539, 0.0001586388))
})
#################################################
#Now test COBRA Methods
test_fba1 <- function(){
  mp  <- system.file(package = "sybil", "extdata")
  mod <- readTSVmod(prefix = "Ec_core", fpath = mp, quoteChar = "\"")
  SYBIL_SETTINGS("SOLVER", "sybilMosek", loadPackage = FALSE)
  #Flux Balance Analysis
  optL <- optimizeProb(mod)
  #Gene Deletion
  optg <- oneGeneDel(mod, algorithm = "fba", fld = "all")
  #Robustness Analysis
  #optrob <- robAna(mod, ctrlreact = "EX_o2(e)")
  #Phenotypic Phase Plane Analysis
  Ec_core_wo_glc <- changeUptake(mod, off = "glc_D[e]")
  optphpp <- phpp(Ec_core_wo_glc,
              ctrlreact = c("EX_succ(e)", "EX_o2(e)"),
              redCosts = TRUE,
              numP = 25,
              verboseMode = 1)
  #Flux Variability Analysis
  optvar <- fluxVar(mod, percentage = 80, verboseMode = 0)
  #MOMA method testing Quadratic Optimization Framework
  fba <- optimizeProb(mod, algorithm = "fba")
  mtf <- optimizeProb(mod, algorithm = "mtf", wtobj = mod_obj(fba))
  ko <- optimizeProb(mod, gene = "b2276", lb = 0, ub = 0,
                     algorithm = "lmoma", wtflux = getFluxDist(mtf))
  
  out <- list(obj = round(lp_obj(optL), 2),
              flux = round(getFluxDist(optL), 2),
              ogene = round(lp_obj(optg), 2),
              oflux = round(getFluxDist(optg), 2),
              #orob = round(lp_obj(optrob), 2),
              #robflux = round(getFluxDist(optrob), 2),
              ophpp = round(lp_obj(optphpp), 2),
              phppflux = round(getFluxDist(optphpp), 2),
              ovar = round(lp_obj(optvar), 2),
              varflux = round(getFluxDist(optvar), 2),
              ovar = round(lp_obj(ko), 2),
              varflux = round(getFluxDist(ko), 2),
              ovar = round(lp_obj(mtf), 2),
              varflux = round(getFluxDist(mtf), 2)
  )
  return(out)
}

test_fba2 <- function(){
  mp  <- system.file(package = "sybil", "extdata")
  mod <- readTSVmod(prefix = "Ec_core", fpath = mp, quoteChar = "\"")
  SYBIL_SETTINGS("SOLVER", "glpkAPI", loadPackage = FALSE)
  SYBIL_SETTINGS("METHOD", "interior")
  #Flux Balance Analysis
  optL <- optimizeProb(mod)
  #Gene Deletion
  optg <- oneGeneDel(mod, algorithm = "fba", fld = "all")
  #Robustness Analysis
  optrob <- robAna(mod, ctrlreact = "EX_o2(e)")
  #Phenotypic Phase Plabe Analysis
  Ec_core_wo_glc <- changeUptake(mod, off = "glc_D[e]")
  optphpp <- phpp(Ec_core_wo_glc,
                  ctrlreact = c("EX_succ(e)", "EX_o2(e)"),
                  redCosts = TRUE,
                  numP = 25,
                  verboseMode = 1)
  #Flux Variability Analysis
  optvar <- fluxVar(mod, percentage = 80, verboseMode = 0)
  #MOMA method testing Quadratic Optimization Framework
  fba <- optimizeProb(mod, algorithm = "fba")
  mtf <- optimizeProb(mod, algorithm = "mtf", wtobj = mod_obj(fba))
  ko <- optimizeProb(mod, gene = "b2276", lb = 0, ub = 0,
                     algorithm = "lmoma", wtflux = getFluxDist(mtf))
  
  out <- list(obj = round(lp_obj(optL), 2),
              flux = round(getFluxDist(optL), 2),
              ogene = round(lp_obj(optg), 2),
              oflux = round(getFluxDist(optg), 2),
              orob = round(lp_obj(optrob), 2),
              robflux = round(getFluxDist(optrob), 2),
              ophpp = round(lp_obj(optphpp), 2),
              phppflux = round(getFluxDist(optphpp), 2),
              ovar = round(lp_obj(optvar), 2),
              varflux = round(getFluxDist(optvar), 2),
              kovar = round(lp_obj(ko), 2),
              koflux = round(getFluxDist(ko), 2),
              mtfvar = round(lp_obj(mtf), 2),
              mtfflux = round(getFluxDist(mtf), 2)
  )
  return(out)
}

outmosek <- test_fba1()
outglpk <- test_fba2()
test_that("cobra methods", {
  expect_equal(outmosek$obj, outglpk$obj)
  expect_equal(outmosek$flux, outglpk$flux)
  expect_equal(outmosek$ogene, outglpk$ogene)
  expect_equal(outmosek$oflux, outglpk$oflux)
  #expect_equal(outmosek$orob, outglpk$orob)
  #expect_equal(outmosek$robflux, outglpk$robflux)
  expect_equal(outmosek$ophpp, outglpk$ophpp)
  expect_equal(outmosek$phppflux, outglpk$phppflux)
  expect_equal(outmosek$ovar, outglpk$ovar)
  expect_equal(outmosek$varflux, outglpk$varflux)
  expect_equal(outmosek$kovar, outglpk$kovar)
  expect_equal(outmosek$koflux, outglpk$koflux)
  expect_equal(outmosek$mtfvar, outglpk$mtfvar)
  expect_equal(outmosek$mtfflux, outglpk$mtfflux)
})