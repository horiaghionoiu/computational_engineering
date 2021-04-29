# Imports
from Mesher import Mesher
from Fields import Fields
import numpy as np

"""
Class DiscretCoeffs.
- In charge of: initialize and store the discretization coefficients for all the domain.
Physical units in International System.
"""


class DiscretCoeffs:
    def __init__(self, a_user_input: dict, a_the_mesh:Mesher = None, a_the_fields:Fields = None):
        self.usr = a_user_input  # The user inputs.
        self.msh = a_the_mesh  # The mesh.
        self.fields = a_the_fields  # The fields.

        self.aE = None
        self.aW = None
        self.aN = None
        self.aS = None
        self.aP = None

    def build_coeffs(self):
        # Build the discretization coefficients.
        self.aE = np.zeros((self.msh.M, self.msh.N))
        self.aW = np.zeros((self.msh.M, self.msh.N))
        self.aN = np.zeros((self.msh.M, self.msh.N))
        self.aS = np.zeros((self.msh.M, self.msh.N))
        self.aP = np.zeros((self.msh.M, self.msh.N))

        numerical_scheme = self.usr.get('numerical_scheme', 'CDS')
        for row in range(self.msh.M):
            for col in range(self.msh.N):
                self.aE[row,col] = self.fields.D_e[row,col] * \
                          DiscretCoeffs.evaluate_convective_scheme(self.fields.P_e[row,col], numerical_scheme) + \
                          max(-self.fields.F_e[row,col],0)
                self.aW[row,col] = self.fields.D_w[row,col] * \
                          DiscretCoeffs.evaluate_convective_scheme(self.fields.P_w[row,col], numerical_scheme) + \
                          max(self.fields.F_w[row,col],0)
                self.aN[row,col] = self.fields.D_n[row,col] * \
                          DiscretCoeffs.evaluate_convective_scheme(self.fields.P_n[row,col], numerical_scheme) + \
                          max(-self.fields.F_n[row,col],0)
                self.aS[row,col] = self.fields.D_s[row,col] * \
                          DiscretCoeffs.evaluate_convective_scheme(self.fields.P_s[row,col], numerical_scheme) + \
                          max(self.fields.F_s[row,col],0)

        # Depending on the problem selected, some boundary conditions apply.
        the_problem = self.usr.get('problem', 'Smith-Hutton')
        if the_problem == 'Diagonal flow':
            # The row below the top wall...
            for col in range(1,self.msh.N-1):
                self.aN[-2,col] = 0

            # The column at the left of the right wall...
            for row in range(1, self.msh.M - 1):
                self.aE[row,-2] = 0
        else:
            # The row above the bottom wall (outlet)...
            for col in range(int(self.msh.cv_y / 2 + 1), self.msh.cv_x):
                self.aS[1,col] = 0

        # Compute the aP discretization coefficient.
        for row in range(self.msh.M):
            for col in range(self.msh.N):
                self.aP[row,col] =  self.aE[row,col] + self.aW[row,col] + self.aN[row,col] + self.aS[row,col]

    @staticmethod
    def evaluate_convective_scheme(a_peclet_number, a_convective_scheme):
        if a_convective_scheme == 'CDS':
            return 1 - 0.5 * abs(a_peclet_number)
        elif a_convective_scheme == 'UDS':
            return 1
        elif a_convective_scheme == 'EDS':
            return abs(a_peclet_number) / (np.exp(abs(a_peclet_number)) - 1)
        else:
            print('DiscretCoeffs.evaluate_convective_scheme."Wrong convective scheme selected."')
            return a_peclet_number