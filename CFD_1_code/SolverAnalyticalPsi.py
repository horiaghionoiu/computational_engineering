# Imports
import numpy as np

from Mesher import Mesher
from Fields import Fields
from PhysicProp import PhysicProp
"""
Class SolverAnalyticalPsi.
- In charge of: building the analytical flow solution.
Physical units in International System.
"""


class SolverAnalyticalPsi:
    def __init__(self, a_the_mesh:Mesher = None, a_the_fields:Fields = None, a_physical_properties:PhysicProp = None):
        self.msh = a_the_mesh
        self.fields = a_the_fields
        self.pprop = a_physical_properties
    
    def solve_psi(self):
        self.fields.Psi_analytical = np.zeros((self.msh.M, self.msh.N))
        rows = self.msh.M
        cols = self.msh.N

        for row in range(rows):
            for col in range(cols):
                if self.msh.fluid_nodes[row, col]:
                    self.fields.Psi_analytical[row, col] =\
                        self.compute_analytical_psi(self.msh,
                                                    self.pprop.u0,
                                                    self.msh.the_nodes_x_coords[row,col],
                                                    self.msh.the_nodes_y_coords[row,col])

    @staticmethod
    def compute_analytical_psi(a_msh: Mesher, a_u0, a_node_x, a_node_y):
        r = np.sqrt((a_node_x-a_msh.cylinder_x)**2 + (a_node_y-a_msh.cylinder_y)**2)
        theta = np.arctan2(a_node_y,a_node_x)
        return a_u0 * np.sin(theta) * (r - a_msh.cylinder_R**2 / r)

