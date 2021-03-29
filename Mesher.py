# Imports
import numpy as np

"""
Class Mesher.
- In charge of: building the mesh.
This case, is an evolution of the channel flow. Here, a solid cylinder is introduced
in the middle of the mesh.
Physical units in International System.
"""


class Mesher:
    def __init__(self):
        # The inputs.

        # Control volumes in x,y.
        self.vc_x = None
        self.vc_y = None
        # Cylinder.
        self.cylinder_x = 0
        self.cylinder_y = 0
        self.cylinder_R = None
        # Mesh size.
        # a_n_Radius_west
        # a_n_Radius_east
        # a_n_Radius_north
        # a_n_Radius_south
        self.x0 = None
        self.y0 = None
        self.xf = None
        self.yf = None
        # Total nodes in x,y.
        self.M = None
        self.N = None

        # Calculated domain parameters.

        # VC size.
        self.dx = None
        self.dy = None
        # VC margins size.
        self.d_margin_x = None
        self.d_margin_y = None
        # Inner domain size.
        self.x0_inner = None
        self.y0_inner = None
        self.xf_inner = None
        self.yf_inner = None
        # The nodes.
        self.the_nodes_x_coords = None
        self.the_nodes_y_coords = None
        self.the_inner_nodes_x_coords = None
        self.the_inner_nodes_y_coords = None
        self.fluid_nodes = None

    def set_control_volumes(self, a_cv_x=50, a_cv_y=30):
        """
        Sets the inner mesh control volumes, it doesn't include the boundary nodes.
        :param a_cv_x: number of x control volumes
        :param a_cv_y: number of y control volumes
        """
        self.vc_x = a_cv_x
        self.vc_y = a_cv_y

        # M would be the number of matrix rows.
        # N would be the number of matrix columns.
        self.M = self.vc_y + 2
        self.N = self.vc_x + 2  # Two nodes extra (begin+end).

    def set_cylinder(self, a_R=1):
        """
        The cylinder is set in the origin.
        :param a_R: radius.
        """
        self.cylinder_R = a_R

    def set_domain_size(self, a_n_Radius_west=-10, a_n_Radius_east=10, a_n_Radius_north=5, a_n_Radius_south=-5):
        """
        Sets the mesh size, as a function of the cylinder radius.
        :param a_n_Radius_west: number of cylinder radius length upwards the cylinder.
        :param a_n_Radius_east: number of cylinder radius length downwards the cylinder.
        :param a_n_Radius_north: number of cylinder radius length above the cylinder.
        :param a_n_Radius_south: number of cylinder radius length below the cylinder.
        """
        self.x0 = self.cylinder_R * a_n_Radius_west
        self.y0 = self.cylinder_R * a_n_Radius_south
        self.xf = self.cylinder_R * a_n_Radius_east
        self.yf = self.cylinder_R * a_n_Radius_north

    def build_mesh(self):
        """
        Build the mesh.
        The outputs of building the mesh are:
            - self.the_nodes_x_coords: the X coords of the mesh nodes.
            - self.the_nodes_y_coords: the Y coords of the mesh nodes.
            - self.dx: the x distance between nodes.
            - self.dy: the y distance between nodes.
            - self.d_margin_x: the x distance between nodes, at the left and right boundaries.
            - self.d_margin_y: the y distance between nodes, at the top and bottom boundaries.
        """
        # VC size.
        self.dx = (self.xf - self.x0 ) / self.vc_x
        self.dy = (self.yf - self.y0 ) / self.vc_y

        # Distance between VC center node and end node.
        self.d_margin_x = self.dx/2
        self.d_margin_y = self.dy/2

        # Inner nodes.
        self.x0_inner = self.x0 + self.d_margin_x
        self.y0_inner = self.y0 + self.d_margin_y
        self.xf_inner = self.xf - self.d_margin_x
        self.yf_inner = self.yf - self.d_margin_y

        # Building the X coords array.
        self.the_inner_nodes_x_coords = np.linspace(self.x0_inner, self.xf_inner, self.vc_x)
        extra_node_begin = self.x0
        extra_node_end = self.xf
        the_nodes_x_coords = np.r_[extra_node_begin, self.the_inner_nodes_x_coords, extra_node_end]
        self.the_nodes_x_coords = np.tile(the_nodes_x_coords, (self.M, 1))

        # Building the Y coords array.
        self.the_inner_nodes_y_coords = np.flip(np.linspace(self.yf_inner, self.y0_inner, self.vc_y))
        extra_node_begin = self.y0
        extra_node_end = self.yf
        the_nodes_y_coords = np.r_[extra_node_begin, self.the_inner_nodes_y_coords, extra_node_end]
        the_nodes_y_coords = np.flip(the_nodes_y_coords)
        the_nodes_y_coords = the_nodes_y_coords[..., None]  # Promote it to column vector aka 2D array.
        self.the_nodes_y_coords = np.tile(the_nodes_y_coords, (1, self.N))

        # Building the fluid nodes.
        self.build_fluid_nodes()

    def build_fluid_nodes(self):
        """
        Builds a matrix of zeros and ones,
        being ones the fluid nodes,
        and zeros the solid (aka cylinder) nodes.
        """
        self.fluid_nodes = np.ones((self.M,self.N))

        for row in range(self.M):
            for col in range(self.N):
                x = self.the_nodes_x_coords[row,col]
                y = self.the_nodes_y_coords[row,col]
                distance_from_node_to_origin = np.sqrt(x**2 + y**2)

                if distance_from_node_to_origin < self.cylinder_R:
                    self.fluid_nodes[row,col] = 0
                # print('(',row,',',col,')')