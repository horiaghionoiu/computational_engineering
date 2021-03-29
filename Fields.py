# Imports
import numpy as np
from CFD_1_PotentialFlows_StaticCylinder.Mesher import Mesher
from CFD_1_PotentialFlows_StaticCylinder.PhysicProp import PhysicProp
"""
Class Fields.
- In charge of: initialize and store the velocity, pressure, temperature, density, and stream function fields. 
These are stored for each node of the domain.
Physical units in International System.
"""


class Fields:
    def __init__(self, a_the_mesh:Mesher = None, a_physical_properties: PhysicProp = None):
        self.msh = a_the_mesh
        self.pprop = a_physical_properties

        # self.v0 = None  # Initial velocity (x)
        # self.v = None  # Total velocity.
        # self.vx = None  # x velocity.
        # self.vy = None  # y velocity.
        self.P0 = None  # Initial pressure.
        self.P = None  # Pressure.
        # self.Cp = None   # Pressure coefficient.
        self.T0 = None  # Initial temperature.
        self.T = None  # Temperature.
        self.rho0 = None  # Initial density.
        self.rho = None  # Density.
        self.Psi = None  # Stream function.
        self.Psi_analytical = None  # Stream function (analytical solution).

    def init_v(self):
        u0 = self.pprop.u0
        N = self.msh.N
        M = self.msh.M
        self.v0 = u0 * np.ones((M,N))
        self.v = u0 * np.ones((M,N))
        self.vx = u0 * np.ones((M,N))
        self.vy = np.zeros((M,N))

    def build_P(self):
        P0 = self.pprop.P0
        # N = self.msh.N
        # M = self.msh.M

        self.P0 = P0 * self.msh.fluid_nodes
        self.P = P0 * self.msh.fluid_nodes
        # self.Cp = np.zeros((M,N))

    def build_T(self):
        T0 = self.pprop.T0
        # N = self.msh.N
        # M = self.msh.M

        self.T0 = T0 * self.msh.fluid_nodes
        self.T = T0 * self.msh.fluid_nodes

    def build_rho(self):
        rho0 = self.pprop.rho0
        # N = self.msh.N
        # M = self.msh.M

        self.rho0 = rho0 * self.msh.fluid_nodes
        self.rho = rho0 * self.msh.fluid_nodes

    def build_Psi(self):
        self.Psi = np.zeros((self.msh.M, self.msh.N))

        u0 = self.pprop.u0
        # Dirichlet condition at the inlet.
        for i in range(self.msh.M):
            self.Psi[i,0] = u0 * self.msh.the_nodes_y_coords[i,0]

        # Dirichlet condition at the top wall.
        self.Psi[0,:] = u0 * (self.msh.yf)

        # Dirichlet condition at the bottom wall.
        self.Psi[-1,:] = u0 * self.msh.the_nodes_y_coords[-1,:]

        # No slip condition inside the solid body.
        self.deactivate_solid_nodes(p_mesh=self.msh, p_field=self.Psi)

    @staticmethod
    def deactivate_solid_nodes(p_mesh=None, p_field=None):
        """
        This method sets to zero all the nodes corresponding to the solid.
        :param a_field: a numpy nd array.
        """
        if not p_field.shape == p_mesh.fluid_nodes.shape:
            print('Fields.deactivate_solid_nodes."Wrong arrays shape"')
            return

        for row in range(p_mesh.M):
            for col in range(p_mesh.N):
                if not p_mesh.fluid_nodes[row,col]:
                    p_field[row,col] = 0  # No slip condition.
