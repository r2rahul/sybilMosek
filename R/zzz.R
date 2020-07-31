# initSybilMosek.R Mosek support for sybil.  Copyright (C) 2014 Rahul,
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


.MSKenv <- new.env()

.onLoad <- function(lib, pkg) {
    
    sybil::addSolver(solver = "sybilMosek", 
            method = "mosek", 
            probType = list(c("lp", "mip", "qp")))
} 
