# Imports
import numpy as np

from CFD_1_PotentialFlows_StaticCylinder.Mesher import Mesher
from CFD_1_PotentialFlows_StaticCylinder.Fields import Fields
from CFD_1_PotentialFlows_StaticCylinder.DiscretCoeffs import DiscretCoeffs
"""
Class SolverGaussSeidelPsi.
- In charge of: solving the potential flow equations, using the stream function (Psi) formulation.
Physical units in International System.
"""


class SolverGaussSeidelPsi:
    def __init__(self, a_epsilon:float = 1e-6, a_the_mesh:Mesher = None, a_the_fields:Fields = None,
                 a_the_coeffs:DiscretCoeffs = None):
        self.epsilon = a_epsilon
        self.msh = a_the_mesh
        self.fields = a_the_fields
        self.coeffs = a_the_coeffs
        self.solution_error = []
        self.solution_iterations = []
    
    def solve_psi(self):
        iter_counter = 0
        iter_error = 1
        # Initial stream function fields, from which we iterate to solve for the stream function.
        psi0 = self.fields.Psi.copy()

        vc_x = self.msh.vc_x
        vc_y = self.msh.vc_y
        while iter_error > self.epsilon:
            # For each row and column, compute Psi.
            # Boundary condition nodes are not solved (top & bottom walls, and inlet).
            for row in range(1,vc_y+1):
                for col in range(1,vc_x+1):
                    if self.msh.fluid_nodes[row,col]:
                        psiE = self.fields.Psi[row,col+1]
                        psiW = self.fields.Psi[row,col-1]
                        psiN = self.fields.Psi[row+1,col]
                        psiS = self.fields.Psi[row-1,col]
                        self.fields.Psi[row,col] = (
                                                    self.coeffs.aE[row-1,col-1]*psiE +
                                                    self.coeffs.aW[row-1,col-1]*psiW +
                                                    self.coeffs.aN[row-1,col-1]*psiN +
                                                    self.coeffs.aS[row-1,col-1]*psiS ) / self.coeffs.aP[row-1,col-1]

            # Neumann bc at the outlet.
            self.fields.Psi[:,self.msh.N-1] = self.fields.Psi[:,self.msh.N-2]

            # Error computation.
            sum_psi0 = np.sum(psi0)
            sum_psi = np.sum(self.fields.Psi)
            iter_error = abs(sum_psi-sum_psi0)

            if not iter_counter % 100:
                print('(Iter: ', iter_counter, ', Error: ',iter_error,')')
            self.solution_iterations.append(iter_counter)
            self.solution_error.append(iter_error)

            psi0 = self.fields.Psi.copy()
            iter_counter = iter_counter + 1
            