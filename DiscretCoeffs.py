# Imports
import numpy as np

from CFD_1_PotentialFlows_StaticCylinder.PhysicProp import PhysicProp
from CFD_1_PotentialFlows_StaticCylinder.Mesher import Mesher
from CFD_1_PotentialFlows_StaticCylinder.Fields import Fields

"""
Class DiscretCoeffs.
- In charge of: initialize and store the stream function discretization coefficients for all the domain.
Physical units in International System.
"""


class DiscretCoeffs:
    def __init__(self, a_the_mesh:Mesher = None, a_physical_properties: PhysicProp = None, a_the_fields:Fields = None):
        self.msh = a_the_mesh
        self.pprop = a_physical_properties
        self.fields = a_the_fields

        self.aE = None
        self.aW = None
        self.aN = None
        self.aS = None
        self.aP = None

    def build_coeffs(self):
        dx = self.msh.dx
        dy = self.msh.dy
        d_margin_x = self.msh.d_margin_x
        d_margin_y = self.msh.d_margin_y
        vc_x = self.msh.vc_x
        vc_y = self.msh.vc_y
        rho0 = self.pprop.rho0

        self.aE = np.zeros((vc_y, vc_x))
        self.aW = np.zeros((vc_y, vc_x))
        self.aN = np.zeros((vc_y, vc_x))
        self.aS = np.zeros((vc_y, vc_x))
        self.aP = np.zeros((vc_y, vc_x))

        rho_rel = np.zeros((self.msh.M, self.msh.N))
        for row in range(self.msh.M):
            for col in range(self.msh.N):
                # Solid nodes will result in inf density, when dividing by zero.
                try:
                    rho_rel[row,col] = rho0 / self.fields.rho[row,col]
                except RuntimeWarning:
                    # Ignore zero division warning.
                    pass

        for row in range(vc_y):
            for col in range(vc_x):
                #   *     *         *             *         *
                #   *    (i,j)     (i,j+1)        *
                #   *   (i+1,j)    (i+1,j+1)   (i+1,j+2)    *
                #   *     *        (i+2,j+2)      *         *
                #   *     *         *             *         *

                # If bottom-right node is fluid.
                if self.msh.fluid_nodes[row+1,col+1]:
                    # If P-right node is fluid.
                    if self.msh.fluid_nodes[row,col+1]:
                        rho_rel_s = dy / (d_margin_y / (rho_rel[row+1,col+1]) + d_margin_y / rho_rel[row,col+1])
                    else:
                        rho_rel_s = rho_rel[row+1,col+1]*2
                    # If P-bottom node is fluid.
                    if self.msh.fluid_nodes[row+1,col]:
                        rho_rel_w = dx / (d_margin_x / (rho_rel[row+1,col+1]) + d_margin_x / rho_rel[row+1,col])
                    else:
                        rho_rel_w = rho_rel[row+1,col+1]*2
                    # If P-bottomright-right node is fluid.
                    if self.msh.fluid_nodes[row+1,col+2]:
                        rho_rel_e = dx / (d_margin_x / rho_rel[row+1,col+1] + d_margin_x / rho_rel[row+1,col+2])
                    else:
                        rho_rel_e = rho_rel[row+1,col+1]*2
                    # If P-bottomright-bottom node is fluid.
                    if self.msh.fluid_nodes[row+2,col+1]:
                        rho_rel_n = dy / (d_margin_y / (rho_rel[row+1,col+1]) + d_margin_y / rho_rel[row+2,col+1])
                    else:
                        rho_rel_n = rho_rel[row+1,col+1]*2
                else:
                    # Bottom-right node is solid.
                    # If P-right node is fluid.
                    if self.msh.fluid_nodes[row, col + 1]:
                        rho_rel_s = rho_rel[row, col+1]*2
                    else:
                        rho_rel_s = 0
                    # If P-bottom node is fluid.
                    if self.msh.fluid_nodes[row + 1, col]:
                        rho_rel_w = rho_rel[row+1, col]*2
                    else:
                        rho_rel_w = 0
                    # If P-bottomright-right node is fluid.
                    if self.msh.fluid_nodes[row + 1, col + 2]:
                        rho_rel_e = rho_rel[row+1, col+2]*2
                    else:
                        rho_rel_e = 0
                    # If P-bottomright-bottom node is fluid.
                    if self.msh.fluid_nodes[row + 2, col + 1]:
                        rho_rel_n = rho_rel[row+2, col+1]*2
                    else:
                        rho_rel_n = 0
                self.aE[row,col] = rho_rel_e * dy/dx
                self.aW[row,col] = rho_rel_w * dy/dx
                self.aN[row,col]= rho_rel_n * dx/dy
                self.aS[row,col] = rho_rel_s * dx/dy

        self.aP = self.aE + self.aW + self.aN + self.aS