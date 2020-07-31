# mosekcompatibility.R Mosek support for sybil.  Copyright (C) 2014 Rahul,
# Department of Biochemistry, McGill University.  All right reserved.
# Email: rahul.rahul@mcgill.ca This file is part of sybilMosek.  sybilMosek
# is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.  SybilMosek is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.  You should have received a copy of the
# GNU General Public License along with sybilMosek.  If not, see
# <http://www.gnu.org/licenses/>.


MSK_STATUS <- list(
    "MSK_RES_OK: No error occurred." = 0
    )

#' @import Rmosek
#' @importFrom sybil optObj
checkSolutionStatus <- function(stat) {
    out <- which(stat != 0)
    return(out)
}

#' @import Rmosek
#' @importFrom sybil optObj
getReturnString <- function(code) {
    out <- ifelse(code == 0, "Solution was successful",
                  "Something went wrong in optimization")
    return(out)
}

#' @import Rmosek
#' @importFrom sybil optObj
getStatusString <- function(code) {
    url <- "https://docs.mosek.com/7.0/capi/Response_codes.html"
    if( code == 0 ){
        out <- names(MSK_STATUS)[1]
    }else{
        out <- paste0("Code ", as.character(code), ": ",
                      "For code details Please check: ", url)
    }
    return(out)
}
