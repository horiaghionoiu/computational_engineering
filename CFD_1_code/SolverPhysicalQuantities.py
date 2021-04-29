# Imports
import numpy as np

from Mesher import Mesher
from Fields import Fields
from DiscretCoeffs import DiscretCoeffs
from PhysicProp import PhysicProp

"""
Class SolverPhysicalQuantities.
- In charge of: evaluating the velocity, pressure, temperature, and density fields,
given a mesh and a Psi stream function.
Physical units in International System.
"""


class SolverPhysicalQuantities:
    def __init__(self, a_the_mesh:Mesher = None,
                 a_the_fields:Fields = None,
                 a_the_coeffs:DiscretCoeffs = None,
                 a_physical_properties:PhysicProp = None):
        self.msh = a_the_mesh
        self.fields = a_the_fields
        self.coeffs = a_the_coeffs
        self.pprop = a_physical_properties

        self.temperature = None

    def solve_velocity(self):
        vc_x = self.msh.vc_x
        vc_y = self.msh.vc_y
        dx = self.msh.dx
        dy = self.msh.dy
        for row in range(vc_y):
            for col in range(vc_x):
                if self.msh.fluid_nodes[row+1,col+1]:
                    psiE = self.fields.Psi[row+1,col+2]
                    psiW = self.fields.Psi[row+1,col]
                    psiN = self.fields.Psi[row+2,col+1]
                    psiS = self.fields.Psi[row,col+1]
                    psiP = self.fields.Psi[row+1,col+1]

                    vye = self.coeffs.aE[row,col]/dy*(psiP-psiE)
                    vyw = self.coeffs.aW[row,col]/dy*(psiW-psiP)
                    vxn = self.coeffs.aN[row,col]/dx*(psiP-psiN)
                    vxs = self.coeffs.aS[row,col]/dx*(psiS-psiP)

                    vx = (vxn + vxs) / 2
                    vy = (vye + vyw) / 2
                    self.fields.vy[row+1,col+1] = vy
                    self.fields.vx[row+1,col+1] = vx
                    self.fields.v[row+1,col+1] = np.sqrt(vx**2 + vy**2)
                else:
                    self.fields.vy[row+1,col+1] = 0
                    self.fields.vx[row+1,col+1] = 0
                    self.fields.v[row+1,col+1] = 0

        # Neumann bc at the outlet.
        self.fields.vy[:, self.msh.N - 1] = self.fields.vy[:, self.msh.N - 2]
        self.fields.vx[:, self.msh.N - 1] = self.fields.vx[:, self.msh.N - 2]
        self.fields.v[:, self.msh.N - 1] = self.fields.v[:, self.msh.N - 2]

    def solve_pressure(self):
        T_relative = self.fields.T / self.fields.T0
        self.fields.P = self.fields.P0 * T_relative ** (self.pprop.gamma / (self.pprop.gamma - 1))

        v_relative = self.fields.v / self.pprop.u0
        self.fields.Cp = 1 - v_relative**2

    def solve_temperature(self):
        self.fields.T = self.fields.T0 + (self.fields.v0**2 - self.fields.v**2)/2/self.pprop.cp
        self.fields.T = self.fields.T * self.msh.fluid_nodes  # Remove temperature for solid nodes.
        self.fields.T[self.fields.T==0] = np.nan  # Better to convert zeros to nans, for plotting.

    def solve_rho(self):
        self.fields.rho = self.fields.P / (self.pprop.R * self.fields.T)
