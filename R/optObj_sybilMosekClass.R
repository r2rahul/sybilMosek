#  optObj_sybilMosekClass.R
# Mosek support for sybil.
#
#  Copyright (C) 2020 Department of Biochemistry, McGill University.
#  All right reserved.
#
#  This file is part of sybilMosek.
#
#  sybilMosek is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SybilMosek is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sybilMosek.  If not, see <http://www.gnu.org/licenses/>.



#------------------------------------------------------------------------------#
#               definition of the class optObj_sybilMosek                     #
#------------------------------------------------------------------------------#

#' An S4 class to represent a optObj for SybilMosek
#
#' @exportClass optObj_sybilMosek
setClass(Class = "optObj_sybilMosek",
         slots = c(msk = "character"),
         contains = "optObj")


#------------------------------------------------------------------------------#
#                                  methods                                     #
#------------------------------------------------------------------------------#
#' An S4 class to remove a optObj class.
#' @exportClass optObj_sybilMosek
setMethod("delProb", signature(lp = "optObj_sybilMosek"),

    function(lp, ...) {

        .MSKenv[[lp@msk]] <- NULL

    }
)


#------------------------------------------------------------------------------#
#' An S4 class to initialize the optObj_sybilMosek class.
#'
#' @slot lp
#'
setMethod("initProb", signature(lp = "optObj_sybilMosek"),

    function(lp, to = FALSE, nrows = 0, ncols = 0 ,...) {

        prob        <- vector(mode = "list", length = 5)
        names(prob) <- c("A", "c", "sense", "bc", "bx")

        repeat{
            pn <- paste(sample(letters, 7), collapse = "")
            if ( (is.null(.MSKenv)) || (! pn %in% ls(.MSKenv)) ) {
                break
            }
        }

        lp@msk <- pn

        if (isTRUE(to)) {
            toflag <- 2
        }
        else {
            toflag <- 2
        }

        .MSKenv[[lp@msk]] <- list(lp   = prob,
                                  opts = list(soldetail = toflag, verbose = 1),
                                  sol  = vector(mode = "list", length = 0))

        .MSKenv[[lp@msk]][["lp"]][["A"]]           <- Matrix(0,
                                                             nrow = nrows,
                                                             ncol = ncols,
                                                             sparse = TRUE)
        .MSKenv[[lp@msk]][["lp"]][["c"]]         <- numeric(ncols)
        .MSKenv[[lp@msk]][["lp"]][["sense"]]       <- "max"
        .MSKenv[[lp@msk]][["lp"]][["bc"]]         <- rbind(numeric(nrows),numeric(nrows))
        .MSKenv[[lp@msk]][["lp"]][["bx"]]          <- rbind(numeric(ncols),numeric(ncols))
        return(lp)
    }
)


#------------------------------------------------------------------------------#

setMethod("backupProb", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        repeat{
            pn <- paste(sample(letters, 7), collapse = "")
            if ( (is.null(.MSKenv)) || (! pn %in% ls(.MSKenv)) ) {
                break
            }
        }

        out <- new("optObj_sybilMosek", lp@solver, lp@method, lp@probType)

        out@msk <- pn

        .MSKenv[[out@msk]] <- .MSKenv[[lp@msk]]

        return(out)
    }
)


#------------------------------------------------------------------------------#


setMethod("setSolverParm", signature(lp = "optObj_sybilMosek"),

    function(lp, solverParm) {
        #https://docs.mosek.com/9.2/rmosek/parameters.html#doc-sparam-list

        if ( !(is.list(solverParm)) ) {
            stop(sQuote(solverParm), " must be list")
        }

        if (any(is.na(solverParm))) {
            stop(sQuote(solverParm), " contains NA values")
        }
        
        #Assumes SolverParam as list
        message("Please provide params as list in format \n
                list(dparam = list(MSK_PARAMS), iparam = list(MSK_PARAMS), \n
                sparam = list(MSK_PARAMS), opts = list(MSK_OPTIONS))\n
                you can provide either of dparam, iparam, sparam, or opts\n
                leave blank if not required like you can only provide\n
                list(opts = list(MSK_OPTIONS)) similarly for others")
        message("For parameter MSK_PARAMS list please check url: \n 
                https://docs.mosek.com/9.2/rmosek/parameters.html#doc-sparam-list \n
                For MSK_Options check help menu ?mosek, like verbose etc.")
        
        if(length(solverParm[['dparam']]) > 0){
            .MSKenv[[lp@msk]][["lp"]][["dparam"]] <- solverParm[['dparam']]
        }
        
        if(length(solverParm[['iparam']]) > 0){
            .MSKenv[[lp@msk]][["lp"]][["iparam"]] <- solverParm[['iparam']]
        }
        
        if(length(solverParm[['sparam']]) > 0){
         .MSKenv[[lp@msk]][["lp"]][["sparam"]] <- solverParm[['sparam']]
        }
        ##################
        #optional params
        if(length(solverParm[['opts']]) > 0){ 
            .MSKenv[[lp@msk]][["opts"]] <- solverParm[['opts']] 
        }
        
        
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolverParm", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        out <- .MSKenv[[lp@msk]][["opts"]]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("setObjDir", signature(lp = "optObj_sybilMosek",
                                 lpdir = "character"),

    function(lp, lpdir) {

        .MSKenv[[lp@msk]][["lp"]][["sense"]] <- ifelse(lpdir == "max",
                                                               "max", "min")

    }
)


#------------------------------------------------------------------------------#

setMethod("getObjDir", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        out <- .MSKenv[[lp@msk]][["lp"]][["sense"]]

        return(out)

    }
)


#------------------------------------------------------------------------------#

setMethod("getNumRows", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        out <- nrow(.MSKenv[[lp@msk]][["lp"]][["A"]])

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumCols", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        out <- ncol(.MSKenv[[lp@msk]][["lp"]][["A"]])

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("addRowsToProb", signature(lp = "optObj_sybilMosek"),

    # i: vector containing the new row indices (must be ascending)
    # cind: list, containing the column indices of the new nz elements
    # nzval: list, containing the new nz elements
    #
    # i, type, lb, ub, cind and nzval must have the same length
    #
    # type can be one of the following:
    # "F" = free variable                -INF <  x <  INF
    # "L" = variable with lower bound      lb <= x <  INF
    # "U" = variable with upper bound    -INF <  x <= ub
    # "D" = double-bounded variable        lb <= x <= ub
    # "E" = fixed variable                 lb  = x  = ub

    function(lp, i, type, lb, ub, cind, nzval, rnames = NULL) {
        
        stopifnot(length(lb) == length(ub))
        #Translate F to Mosek
        f <- which(type == "F")
        lb[f] <- -Inf
        ub[f] <- Inf
        #Translate L type
        l <- which(type == "L")
        ub[l] <- Inf
        #Translate U type
        u <- which(type == "U")
        lb[u] <- -Inf
        

        nc  <- ncol(.MSKenv[[lp@msk]][["lp"]][["A"]])
        nr  <- nrow(.MSKenv[[lp@msk]][["lp"]][["A"]])
        mat <- Matrix(0, nrow = length(i), ncol = nc)

        .MSKenv[[lp@msk]][["lp"]][["A"]] <- rbind(
                                       .MSKenv[[lp@msk]][["lp"]][["A"]], mat)

        bc <- rbind(
            append(
                .MSKenv[[lp@msk]][["lp"]][["bc"]][1, ], rep(0, length(i))),
            append(
                .MSKenv[[lp@msk]][["lp"]][["bc"]][2, ], rep(0, length(i)))
        )
        
        .MSKenv[[lp@msk]][["lp"]][["bc"]] <- bc

        for (k in seq(along = i)) {
            
            .MSKenv[[lp@msk]][["lp"]][["A"]][(nr+k), cind[[k]]] <- nzval[[k]]

            .MSKenv[[lp@msk]][["lp"]][["bc"]][1,(nr+k)] <- lb[k]

			.MSKenv[[lp@msk]][["lp"]][["bc"]][2,(nr+k)] <- ub[k]
        }


    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBnds", signature(lp = "optObj_sybilMosek"),

    function(lp, j, lb, ub) {

        .MSKenv[[lp@msk]][["lp"]][["bx"]][1, j] <- lb
        .MSKenv[[lp@msk]][["lp"]][["bx"]][2, j] <- ub
    }
)


#------------------------------------------------------------------------------#

setMethod("changeColsBndsObjCoefs", signature(lp = "optObj_sybilMosek"),

    function(lp, j, lb, ub, obj_coef) {
        .MSKenv[[lp@msk]][["lp"]][["bx"]][1, j]  <- lb
        .MSKenv[[lp@msk]][["lp"]][["bx"]][2, j]  <- ub
        .MSKenv[[lp@msk]][["lp"]][["c"]][j] <- obj_coef

    }
)


#------------------------------------------------------------------------------#

setMethod("getColsLowBnds", signature(lp = "optObj_sybilMosek"),

    function(lp, j) {

        out <- .MSKenv[[lp@msk]][["lp"]][["bx"]][1,j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColsUppBnds", signature(lp = "optObj_sybilMosek"),

    function(lp, j) {

        out <- .MSKenv[[lp@msk]][["lp"]][["bx"]][2,j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeRowsBnds", signature(lp = "optObj_sybilMosek"),

    function(lp, i, lb, ub) {
        .MSKenv[[lp@msk]][["lp"]][["bc"]][1,i] <- lb
        .MSKenv[[lp@msk]][["lp"]][["bc"]][2,i] <- ub
    }
)


#------------------------------------------------------------------------------#

setMethod("setRhsZero", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        nr  <- nrow(.MSKenv[[lp@msk]][["lp"]][["A"]])
        #Assuming Sv=0 constraint
        .MSKenv[[lp@msk]][["lp"]][["bc"]][1, ]   <- rep(0, nr)
        .MSKenv[[lp@msk]][["lp"]][["bc"]][2, ]   <- rep(0, nr)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeObjCoefs", signature(lp = "optObj_sybilMosek"),

    function(lp, j, obj_coef) {

        .MSKenv[[lp@msk]][["lp"]][["c"]][j] <- obj_coef
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjCoefs", signature(lp = "optObj_sybilMosek"),

    function(lp, j) {

        out <- .MSKenv[[lp@msk]][["lp"]][["c"]][j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("changeMatrixRow", signature(lp = "optObj_sybilMosek"),

    function(lp, i, j, val) {

        .MSKenv[[lp@msk]][["lp"]][["A"]][i, j] <- val

    }
)


#------------------------------------------------------------------------------#

setMethod("loadLPprob", signature(lp = "optObj_sybilMosek"),

    function(lp, nCols, nRows, mat, ub, lb, obj, rlb, rtype,
             lpdir = "max", rub = NULL, ctype = NULL,
             cnames = NULL, rnames = NULL, pname = NULL) {

        stopifnot(is(mat, "Matrix"))

        # optimization direction
        setObjDir(lp, lpdir = lpdir)
        
        #Translate variable bounds for the binary variable
        # variable type
        if(!is.null(ctype)){
            vartype <- which((ctype == "I") | (ctype == "B"))
            .MSKenv[[lp@msk]][["lp"]][["intsub"]] <- vartype
            bflag <- which(ctype == "B")
            lb[bflag] <- 0
            ub[bflag] <- 1
        }

        # constraint matrix
        .MSKenv[[lp@msk]][["lp"]][["A"]] <- as(mat, "CsparseMatrix")

        # column (variable) bounds and objective function
        .MSKenv[[lp@msk]][["lp"]][["bx"]]  <- rbind(lb,ub)
        .MSKenv[[lp@msk]][["lp"]][["c"]] <- obj
        
        #Linear Constraints on matrix A
        # type can be one of the following:
        # "F" = free variable                -INF <  x <  INF
        # "L" = variable with lower bound      lb <= x <  INF
        # "U" = variable with upper bound    -INF <  x <= ub
        # "D" = double-bounded variable        lb <= x <= ub
        # "E" = fixed variable                 lb  = x  = ub
        if(is.null(rub)){
            #message("No upper row bound provided assuming L type lb <= x <  INF")
            rub <- rlb
        }else{ 
            stopifnot(length(rlb) == length(rub))
        }
        #Translate F to Mosek
        f <- which(rtype == "F")
        rlb[f] <- -Inf
        rub[f] <- Inf
        #Translate L type
        l <- which(rtype == "L")
        rub[l] <- Inf
        #Translate U type
        u <- which(rtype == "U")
        rlb[u] <- -Inf
        .MSKenv[[lp@msk]][["lp"]][["bc"]] <- rbind(rlb, rub)
    }
)


#------------------------------------------------------------------------------#

setMethod("loadQobj", signature(lp = "optObj_sybilMosek", mat = "Matrix"),

    function(lp, mat) {
        #browser()
        mat <- tril(mat)
        ind <- which(as.matrix(!(mat == 0)), arr.ind = TRUE)
        .MSKenv[[lp@msk]][["lp"]][["qobj"]][["i"]] <- ind[, 1]
        .MSKenv[[lp@msk]][["lp"]][["qobj"]][["j"]] <- ind[, 2]
        .MSKenv[[lp@msk]][["lp"]][["qobj"]][["v"]] <- mat@x

    }
)

#------------------------------------------------------------------------------#

setMethod("loadQobj", signature(lp = "optObj_sybilMosek", mat = "numeric"),
          
          function(lp, mat) {
              
              ind <- 1:length(mat)
              .MSKenv[[lp@msk]][["lp"]][["qobj"]][["i"]] <- ind
              .MSKenv[[lp@msk]][["lp"]][["qobj"]][["j"]] <- ind
              .MSKenv[[lp@msk]][["lp"]][["qobj"]][["v"]] <- mat
              
          }
)

#------------------------------------------------------------------------------#

setMethod("scaleProb", signature(lp = "optObj_sybilMosek"),

    function(lp, opt) {
        
        if(length(.MSKenv[[lp@msk]][["iparam"]]) > 0 ){
            .MSKenv[[lp@msk]][["iparam"]][["SIM_SCALING_METHOD"]] <- "MSK_SCALING_METHOD_FREE"
            .MSKenv[[lp@msk]][["iparam"]][["INTPNT_SCALING"]] <- "MSK_SCALING_AGGRESSIVE"
        } else{
            .MSKenv[[lp@msk]][["iparam"]] <- list(
                SIM_SCALING_METHOD = "MSK_SCALING_METHOD_FREE",
                INTPNT_SCALING = "MSK_SCALING_AGGRESSIVE")
        }

    }
)


#------------------------------------------------------------------------------#

setMethod("solveLp", signature(lp = "optObj_sybilMosek"),

    function(lp) {
        
        #For backward Compatability
        sol_param <- list()
        if (length(.MSKenv[[lp@msk]][["lp"]][["dparam"]]) > 0) {
            sol_param[["dparam"]]  <- c(sol_param , .MSKenv[[lp@msk]][["lp"]][["dparam"]])
        } 
        if(length(.MSKenv[[lp@msk]][["lp"]][["iparam"]]) > 0){
            sol_param[["iparam"]]  <- c(sol_param , .MSKenv[[lp@msk]][["lp"]][["iparam"]])
        }
        if(length(.MSKenv[[lp@msk]][["lp"]][["sparam"]]) > 0){
            sol_param[["sparam"]]  <- c(sol_param , .MSKenv[[lp@msk]][["lp"]][["sparam"]])
        }
        #Set opts field for the solver
        if(length(.MSKenv[[lp@msk]][["opts"]]) > 0){
            opt <- .MSKenv[[lp@msk]][["opts"]]
            if(opt$soldetail <= 1.5){
                message("Option soldetail should be great than 1, setting it to 2")
                opt$soldetail <- 2
            }
        }else{
            #Set Solution details param and information to print
            opt <- list(soldetail = 2, verbose = 1)
        }
            
        #For backward Compatability
        .MSKenv[[lp@msk]][["opts"]] <- sol_param 
        rm(sol_param)
        
        
        .MSKenv[[lp@msk]][["sol"]] <- Rmosek::mosek(.MSKenv[[lp@msk]][["lp"]], opts = opt)

        .MSKenv[[lp@msk]][["sol"]][["status"]] <- .MSKenv[[lp@msk]][["sol"]][["response"]][["code"]]
        #msg <- .MSKenv[[lp@msk]][["sol"]][["response"]][["code"]]
        
        if(lp@probType == "mip"){
            .MSKenv[[lp@msk]][["sol"]][["sol"]][["itr"]] <- 
                .MSKenv[[lp@msk]][["sol"]][["sol"]][["int"]]
        }
        
        objval <- .MSKenv[[lp@msk]][["sol"]][["sol"]][["itr"]][["pobjval"]]
        if (is.null(objval)) {
            out <- 1
            .MSKenv[[lp@msk]][["sol"]][["sol"]][["itr"]][["pobjval"]] <- 0
        }
        else {
            out <- 0
        }
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getObjVal", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        out <- .MSKenv[[lp@msk]][["sol"]][["sol"]][["itr"]][["pobjval"]]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getRedCosts", signature(lp = "optObj_sybilMosek"),

    function(lp) {
         slx <- .MSKenv[[lp@msk]][["sol"]][["sol"]][["itr"]][["slx"]]
         sux <- .MSKenv[[lp@msk]][["sol"]][["sol"]][["itr"]][["sux"]]
         out <- slx - sux
         out <- ifelse(abs(out) < 1e-6, 0, out)
        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getSolStat", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        out <- .MSKenv[[lp@msk]][["sol"]][["status"]] 

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getFluxDist", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        fld <- .MSKenv[[lp@msk]][["sol"]][["sol"]][["itr"]][["xx"]]
        if (is.null(fld)) {
            out <- rep(0, getNumCols(lp))
        }
        else {
            out <- fld
        }

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getColPrim", signature(lp = "optObj_sybilMosek"),

    function(lp, j) {

        out <- .MSKenv[[lp@msk]][["sol"]][["sol"]][["itr"]][["xx"]][j]

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("getNumNnz", signature(lp = "optObj_sybilMosek"),

    function(lp) {

        out <- nnzero(.MSKenv[[lp@msk]][["lp"]][["A"]])

        return(out)
    }
)


#------------------------------------------------------------------------------#

setMethod("writeProb", signature(lp = "optObj_sybilMosek"),

    function(lp, fname, ff = "lp") {

        if (length(.MSKenv[[lp@msk]][["opts"]]) > 0) {
            gp <- .MSKenv[[lp@msk]][["opts"]]
        }
        else {
            gp <- vector(mode = "list", length = 0)

        }

        if (grepl(".", fname, fixed = TRUE)) {
            fn <- fname
        }
        else {
            fn <- paste(fname, ff, sep = ".")
        }

        gp[["writeafter"]] <- fn
        Rmosek::mosek(.MSKenv[[lp@msk]][["lp"]], opts = gp)


        return(TRUE)
    }
)
