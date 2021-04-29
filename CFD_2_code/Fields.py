# Imports
from Mesher import Mesher
import numpy as np

"""
Class Fields.
- In charge of: initialize and store the problem fields. 
Physical units in International System.
"""


class Fields:
    def __init__(self, a_user_input: dict, a_the_mesh: Mesher):
        self.usr = a_user_input  # The user inputs.
        self.msh = a_the_mesh  # The mesh.

        # Velocity field.
        the_problem = self.usr.get('problem', 'Smith-Hutton')
        if the_problem == 'Diagonal flow':
            # In diagonal flow problem, the velocity field is the same for all nodes.
            self.vx = None
            self.vy = None
        else:
            # In Smith-Hutton problem, the velocity field depends on each node's coords.
            # The velocity field is known (class notes) and is computed at each cell face.
            self.vx = None  # At nodes.
            self.vy = None  # At nodes.
            self.vx_e = None  # At cell faces.
            self.vx_w = None  # At cell faces.
            self.vy_n = None  # At cell faces.
            self.vy_s = None  # At cell faces.

        # Phi.
        self.Phi = None

        # Convection strength.
        self.F_e = None
        self.F_w = None
        self.F_n = None
        self.F_s = None

        # Diffusion strength.
        self.D_e = None
        self.D_w = None
        self.D_n = None
        self.D_s = None

        # Peclet number.
        self.P_e = None
        self.P_w = None
        self.P_n = None
        self.P_s = None

    def build_velocity(self):
        # Smith-Hutton or Diagonal flow problem.
        the_problem = self.usr.get('problem', 'Smith-Hutton')
        if the_problem == 'Diagonal flow':
            alpha = self.usr.get('alpha_diagonal_flow', 45) * (np.pi/180)
            v0 = self.usr.get('v0_diagonal_flow', 10)
            self.vx = v0 * np.cos(alpha)
            self.vy = v0 * np.sin(alpha)
        else:
            self.vx = np.zeros((self.msh.M, self.msh.N))
            self.vy = np.zeros((self.msh.M, self.msh.N))
            self.vx_e = np.zeros((self.msh.M, self.msh.N))
            self.vx_w = np.zeros((self.msh.M, self.msh.N))
            self.vy_n = np.zeros((self.msh.M, self.msh.N))
            self.vy_s = np.zeros((self.msh.M, self.msh.N))

            for row in range(self.msh.M):
                for col in range(self.msh.N):
                    x = self.msh.cells_x_coords[row,col]
                    y = self.msh.cells_y_coords[row,col]
                    self.vx[row,col] = + 2 * y * (1 - x**2)
                    self.vy[row,col] = - 2 * x * (1 - y**2)

                    self.vx_e[row,col] = 2 * self.msh.y_e[row,col] * (1 - self.msh.x_e[row,col]**2)
                    self.vx_w[row,col] = 2 * self.msh.y_w[row,col] * (1 - self.msh.x_w[row,col]**2)
                    self.vy_n[row,col] = - 2 * self.msh.x_n[row,col] * (1 - self.msh.y_n[row,col]**2)
                    self.vy_s[row,col] = - 2 * self.msh.x_s[row,col] * (1 - self.msh.y_s[row,col]**2)

    def build_phi(self):
        # Build phi and apply the boundary conditions, depending on the problem chosen. (Smith-Hutton/Diagonal flow)
        self.Phi = np.zeros((self.msh.M, self.msh.N))

        the_problem = self.usr.get('problem', 'Smith-Hutton')
        if the_problem == 'Diagonal flow':
            # Boundary conditions for Diagonal flow case, extracted from class notes.
            # x = 0 (left wall) & y=1 (top wall).
            self.Phi[:, 0] = 1
            self.Phi[-1, :] = 1
            # y = 0 (lower wall) & x=1 (right wall).
            self.Phi[0, :] = 0
            self.Phi[:, -1] = 0
            # x = 0, y = 0 (left bottom corner).
            self.Phi[0, 0] = 0.5
        else:
            # Boundary conditions for Smith-Hutton case, extracted from class notes.
            # Alpha is 10, extracted from class notes.
            alpha = 10

            # Get the inlet indices, cutremente.
            x_coords = self.msh.cells_x_coords[0,:]
            inlet_ix = 0
            for i in range(len(x_coords)):
                if x_coords[i] <= 0:
                    inlet_ix = i

            # Inlet -1<=x<=0 y=0 (lower left wall).
            for col in range(inlet_ix):
                self.Phi[0,col] = 1 + np.tanh(alpha * (2 * self.msh.cells_x_coords[0,col] + 1))
            # x = 0 (left wall).
            self.Phi[:, 0] = 1 - np.tanh(alpha)
            # x = 1 (right wall).
            self.Phi[:, -1] = 1 - np.tanh(alpha)
            # y = 1 (top wall).
            self.Phi[-1,:] = 1 - np.tanh(alpha)

    def build_convection_strength(self):
        self.F_e = np.zeros((self.msh.M, self.msh.N))
        self.F_w = np.zeros((self.msh.M, self.msh.N))
        self.F_n = np.zeros((self.msh.M, self.msh.N))
        self.F_s = np.zeros((self.msh.M, self.msh.N))

        the_problem = self.usr.get('problem', 'Smith-Hutton')
        rho = self.usr.get('rho', 1e9)
        for row in range(self.msh.M):
            for col in range(self.msh.N):
                if the_problem == 'Diagonal flow':
                    self.F_e[row,col] = rho * self.msh.dy * self.vx
                    self.F_w[row,col] = rho * self.msh.dy * self.vx
                    self.F_n[row,col] = rho * self.msh.dx * self.vy
                    self.F_s[row,col] = rho * self.msh.dx * self.vy
                else:
                    self.F_e[row,col] = rho * self.msh.dy * self.vx_e[row,col]
                    self.F_w[row,col] = rho * self.msh.dy * self.vx_w[row,col]
                    self.F_n[row,col] = rho * self.msh.dx * self.vy_n[row,col]
                    self.F_s[row,col] = rho * self.msh.dx * self.vy_s[row,col]

    def build_diffusion_strength(self):
        self.D_e = np.zeros((self.msh.M, self.msh.N))
        self.D_w = np.zeros((self.msh.M, self.msh.N))
        self.D_n = np.zeros((self.msh.M, self.msh.N))
        self.D_s = np.zeros((self.msh.M, self.msh.N))

        # The conductance.
        gamma = self.usr.get('gamma', 1)

        for row in range(self.msh.M):
            for col in range(self.msh.N):
                self.D_e[row,col] = gamma * self.msh.dy / self.msh.dx
                self.D_w[row,col] = gamma * self.msh.dy / self.msh.dx
                self.D_n[row,col] = gamma * self.msh.dx / self.msh.dy
                self.D_s[row,col] = gamma * self.msh.dx / self.msh.dy

    def build_peclet_number(self):
        # self.P_e = np.zeros((self.msh.M, self.msh.N))
        # self.P_w = np.zeros((self.msh.M, self.msh.N))
        # self.P_n = np.zeros((self.msh.M, self.msh.N))
        # self.P_s = np.zeros((self.msh.M, self.msh.N))

        # for row in range(self.msh.M):
        #     for col in range(self.msh.N):
        #         pass
        self.P_e = self.F_e / self.D_e
        self.P_w = self.F_w / self.D_w
        self.P_n = self.F_n / self.D_n
        self.P_s = self.F_s / self.D_s
