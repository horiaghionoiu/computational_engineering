# Imports
import numpy as np

"""
Class Mesher.
- In charge of: building the mesh.
Physical units in International System.
"""


class Mesher:
    def __init__(self, a_user_input: dict):
        self.usr = a_user_input  # The user inputs.

        self.cv_x = None  # Control volumes in x direction.
        self.cv_y = None  # Control volumes in y direction.
        self.x0 = None  # Mesh x origin.
        self.y0 = None  # Mesh y origin.
        self.xf = None  # Mesh x end.
        self.yf = None  # Mesh y end.

        self.L = None  # Mesh length.
        self.H = None  # Mesh height.
        self.dx = None  # Distance between nodes in x direction.
        self.dy = None  # Distance between nodes in y direction.
        self.M = None  # Total nodes in x direction.
        self.N = None  # Total nodes in y direction.

        # The cells coords.
        self.cells_x_coords = None
        self.cells_y_coords = None

        #  The cell faces coords.
        self.x_w = None
        self.x_e = None
        self.x_n = None
        self.x_s = None

        self.y_w = None
        self.y_e = None
        self.y_n = None
        self.y_s = None

    def build_mesh(self):
        """
        Build the mesh.
        """
        self.load_user_inputs()
        self.compute_mesh()

    def load_user_inputs(self):
        the_problem = self.usr.get('problem', 'Smith-Hutton')

        # The domain size.
        if the_problem == 'Diagonal flow':
            self.x0 = 0
            self.y0 = 0
            self.xf = 1
            self.yf = 1
        else:
            self.x0 = -1
            self.y0 = 0
            self.xf = 1
            self.yf = 1

        self.L = self.xf - self.x0
        self.H = self.yf - self.y0

        # The control volumes.
        self.cv_y = self.usr.get('cv_y', 100)
        self.cv_x = self.usr.get('cv_x', 100)

        # The nodes distance.
        self.dx = self.L / self.cv_x
        self.dy = self.H / self.cv_y

        # Nodes matrix size.
        self.M = self.cv_y + 1
        self.N = self.cv_x + 1

    def compute_mesh(self):
        # Compute the cells coords.
        self.cells_x_coords = np.zeros((self.M, self.N))
        self.cells_y_coords = np.zeros((self.M, self.N))

        for row in range(self.M):
            for col in range(self.N):
                self.cells_x_coords[row, col] = self.x0 + col * self.dx
                self.cells_y_coords[row, col] = self.y0 + row * self.dy

        # Compute the cell faces coords.
        self.x_w = np.zeros((self.M, self.N))
        self.x_e = np.zeros((self.M, self.N))
        self.x_n = np.zeros((self.M, self.N))
        self.x_s = np.zeros((self.M, self.N))

        self.y_w = np.zeros((self.M, self.N))
        self.y_e = np.zeros((self.M, self.N))
        self.y_n = np.zeros((self.M, self.N))
        self.y_s = np.zeros((self.M, self.N))

        for row in range(self.M):
            for col in range(self.N):
                self.x_w[row,col] = self.cells_x_coords[row,col] - self.dx / 2
                self.x_e[row,col] = self.cells_x_coords[row,col] + self.dx / 2
                self.x_n[row,col] = self.cells_x_coords[row,col]
                self.x_s[row,col] = self.cells_x_coords[row,col]
                self.y_w[row,col] = self.cells_y_coords[row,col]
                self.y_e[row,col] = self.cells_y_coords[row,col]
                self.y_n[row,col] = self.cells_y_coords[row,col] + self.dy / 2
                self.y_s[row,col] = self.cells_y_coords[row,col] - self.dy / 2