#   gauss.py
#   Supporting material software to the article 
#    "Probing single cell fermentation flux and intercellular exchange networks via
#    pH-microenvironment sensing and inverse modeling"
#    
#    Copyright (C) 2022 V. Onesto et al.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
                                                                         


import numpy as np
import quadprog
def quadprog_solve_qp(P, q, G=None, h=None, A=None, b=None):
    qp_G = .5 * (P + P.T) 
    qp_a = -q
    if A is not None:
        qp_C = -numpy.vstack([A, G]).T
        qp_b = -numpy.hstack([b, h])
        meq = A.shape[0]
    else:  
        qp_C = -G.T
        qp_b = -h
        meq = 0
    return quadprog.solve_qp(qp_G, qp_a, qp_C, qp_b, meq)[0]
A = np.loadtxt("A.dat")                   # INPUT: inverse distance probes-cell matrix   
PH = np.loadtxt("ph.dat")
ERR= np.loadtxt("err.dat")
C = np.power(10,-PH)
CERR = np.log(10)*C*ERR
Ap = A.transpose()/CERR
M = Ap.transpose()
P = np.dot(M.transpose(),M) 
gamma2 = 1e9                               # Lagrange Multiplier see SI 3.2.1
deltaI = np.identity(np.sqrt(P.size).astype(int))
deltaIs = np.sqrt(P.size).astype(int)
deltaI[deltaIs-1,deltaIs-1]=0
P = P + gamma2*deltaI
b= C/CERR
q = -np.dot(M.transpose(), b)
G = -M
en2 = M[0].size
sizef = b.size
h = np.zeros(sizef).reshape((sizef,))
flux=quadprog_solve_qp(P, q, G, h)
factor=800000                            
np.savetxt("flux.dat",flux*factor)
COV = np.linalg.inv(0.5*P)-np.outer(flux,flux) 
np.savetxt("cov.dat",factor*factor*COV)
