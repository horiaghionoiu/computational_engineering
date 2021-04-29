# Imports
from Mesher import Mesher
from Fields import Fields
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

"""
Class Plotter.
- In charge of: plotting interesting quantities.
"""


class Plotter:
    def __init__(self, a_user_input,
                 a_the_mesh = None,
                 a_the_fields = None,
                 a_reference_solution={}):
        self.usr = a_user_input  # The user inputs.
        self.msh = a_the_mesh  # The mesh.
        self.fields = a_the_fields  # The fields.
        self.reference_solution = a_reference_solution  # The reference solution, from class notes.
    # @staticmethod
    # def map_fluid_nodes(a_mesh:Mesher, nodes_x:np.array, nodes_y:np.array):
    #
    #     fluid_nodes = np.ones((nodes_x.shape[0],nodes_y.shape[1]))
    #
    #     for row in range(nodes_x.shape[0]):
    #         for col in range(nodes_x.shape[1]):
    #             x = nodes_x[row,col]
    #             y = nodes_y[row,col]
    #             distance_from_node_to_origin = np.sqrt(x**2 + y**2)
    #
    #             if distance_from_node_to_origin < a_mesh.cylinder_R:
    #                 fluid_nodes[row,col] = 0.000000000000001
    #             # print('(',row,',',col,')')
    #     return fluid_nodes.copy()

    def plot_mesh(self):
        # Flatten the nodes matrixes.
        x = np.concatenate(self.msh.cells_x_coords).ravel().tolist()
        y = np.concatenate(self.msh.cells_y_coords).ravel().tolist()

        # Get the boundary nodes.
        x_boundary_nodes = list(self.msh.cells_x_coords[0,:]) + \
                           list(self.msh.cells_x_coords[-1,:]) + \
                           list(self.msh.cells_x_coords[:,0]) + \
                           list(self.msh.cells_x_coords[:,-1])
        y_boundary_nodes = list(self.msh.cells_y_coords[0,:]) + \
                           list(self.msh.cells_y_coords[-1,:]) + \
                           list(self.msh.cells_y_coords[:,0]) + \
                           list(self.msh.cells_y_coords[:,-1])
        # Plot the nodes.
        fig = go.Figure(data=go.Scatter(x=x, y=y, mode='markers', name='Internal nodes'))
        fig.add_scatter(x=x_boundary_nodes, y=y_boundary_nodes, mode='markers', name='Boundary nodes')

        # Figure config.
        separator = '-'
        space_m = '\:'
        # Text display hack in case of Smith-Hutton problem.
        problem = self.usr.get('problem','')
        if problem == 'Smith-Hutton':
            problem = 'Smith \: Hutton'
        title = 'Mesh' + separator + problem + separator + str(self.usr.get('cv_y','')) + 'x' + \
        str(self.usr.get('cv_x','')) + ' nodes' + separator + '\\Delta_x' + space_m + "{:.4f}".format(self.msh.dx) + \
                space_m + '\\Delta_y' + \
                space_m + "{:.4f}".format(self.msh.dy)
        title = '$' + title + '$'
        fig.update_layout( title={
                                'text': title,
                                'x':0.05,
                                'xanchor': 'left',
                                'yanchor': 'top'},
                          xaxis_title="x (m)",
                          yaxis_title="y (m)",
                          font={'size':15})
        fig.show()
        # pio.write_image(fig, files_output_path + self.output_files_name + 'Mesh.png', width=2864, height=1356)

    def plot_quiver(self, p_meshgrid_spacing=0.02):
        # Arrange data.
        x, y = np.meshgrid(np.arange(self.msh.x0, self.msh.xf, p_meshgrid_spacing),
                           np.arange(self.msh.y0, self.msh.yf, p_meshgrid_spacing))

        the_problem = self.usr.get('problem', 'Smith-Hutton')
        if the_problem == 'Diagonal flow':
            u = np.ones(x.shape) * self.fields.vx
            v = np.ones(x.shape) * self.fields.vy

            the_scale = 0.0025
            the_arrow_scale = 0.5
        else:
            # Prepare the mesh and velocity field for the quiver plot function.
            n_nodes = self.msh.M * self.msh.N

            # Rearrange the nodes x,y coords in a single array, and flatten the vx,vy velocities.
            the_points = np.zeros((n_nodes, 2))
            u = np.zeros((n_nodes,))
            v = np.zeros((n_nodes,))
            i = 0
            for row in range(self.msh.M):
                for col in range(self.msh.N):
                    the_points[i,0] = self.msh.cells_x_coords[row,col]
                    the_points[i,1] = self.msh.cells_y_coords[row,col]
                    u[i] = self.fields.vx[row,col]
                    v[i] = self.fields.vy[row,col]
                    i = i + 1

            x_ = np.linspace(self.msh.x0, self.msh.xf, 50)
            y_ = np.linspace(self.msh.y0, self.msh.yf, 50)
            x, y = np.meshgrid(x_, y_)
            u = griddata(the_points, u, (x, y), method='cubic')
            v = griddata(the_points, v, (x, y), method='cubic')

            the_scale = 0.02
            the_arrow_scale = 0.3


        # Build the figure.
        fig = ff.create_quiver(x, y, u, v, name='Velocity field', scale=the_scale,
                               arrow_scale=the_arrow_scale)

        # Plot origin.-
        fig.add_trace(go.Scatter(x=[0], y=[0],
                                 name='Origin',
                                 mode='markers'#,
                                 # marker_size=15
                                 ))
        # Figure config.
        separator = '-'
        space_m = '\:'
        problem = self.usr.get('problem','')
        if problem == 'Smith-Hutton':
            problem = 'Smith \: Hutton'
        title = 'Quiver' + separator + problem
        title = '$' + title + '$'
        fig.update_layout( title={
                                'text': title,
                                'x':0.05,
                                'xanchor': 'left',
                                'yanchor': 'top'},
                          xaxis_title="x (m)",
                          yaxis_title="y (m)",
                          font={'size':15})
        fig.show()

    def plot_phi_vs_x(self):
        # Plot reference solution.
        x_reference_solution = self.reference_solution['x']
        y_reference_solution = None

        rho_gamma = self.usr.get('rho_gamma',10)
        if rho_gamma == 1e1:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e1', None)
        elif rho_gamma == 1e3:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e3', None)
        elif rho_gamma == 1e6:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e6', None)

        fig = go.Figure(data=go.Scatter(x=x_reference_solution, y=y_reference_solution, line=dict(color='firebrick'), name='Reference solution'))

        # Plot numerical solution.
        x = self.msh.cells_x_coords[0,:].tolist()[self.get_outlet_indices()::]
        y = self.fields.Phi[0,:].tolist()[self.get_outlet_indices()::]
        fig.add_scatter(x=x, y=y, line=dict(color='royalblue'),name='Numerical solution')

        # Figure config.
        separator = '-'
        space_m = '\:'
        # Text display hack in case of Smith-Hutton problem.
        problem = self.usr.get('problem','')
        if problem == 'Smith-Hutton':
            problem = 'Smith \: Hutton'
        title = '\\phi' + space_m + 'value' + space_m + '@' + space_m + 'outlet' + separator + problem + separator + str(self.usr.get('cv_y','')) + 'x' + \
        str(self.usr.get('cv_x','')) + ' nodes' + separator + '\\Delta_x' + space_m + "{:.4f}".format(self.msh.dx) + \
                space_m + '\\Delta_y' + space_m + "{:.4f}".format(self.msh.dy) + \
                separator + '\\rho / \\Gamma = ' + str(rho_gamma)
        title = '$' + title + '$'
        fig.update_layout( title={
                                'text': title,
                                'x':0.05,
                                'xanchor': 'left',
                                'yanchor': 'top'},
                          xaxis_title="x (m)",
                          yaxis_title="$\\phi$",
                          font={'size':15})

        fig.show()

    def plot_phi_contours(self):
        field_to_plot = self.fields.Phi
        title1 = "$\\phi$"
        rho_gamma = self.usr['rho']
        space_m = '\:'
        separator = '-'
        title2 = '\\phi' + space_m + ' value' + space_m + ' contours' + separator + \
                '\\rho / \\Gamma = ' + str(rho_gamma)
        title2 = '$' + title2 + '$'
        fig = go.Figure(
            layout=go.Layout(
                title=go.layout.Title(text=title2)
            ))
        fig.add_trace(go.Contour(
            z=field_to_plot,
            colorscale="aggrnyl",
            colorbar=dict(
                title=title1
            )
        ))

        # fig.update_layout()
        fig.update_xaxes(title_text="cv_x")
        fig.update_yaxes(title_text="cv_y")
        fig.update_layout(font={'size': 15})

        fig.show()

    def get_outlet_indices(self):
        # Get the outlet first index, cutremente.
        x_coords = self.msh.cells_x_coords[0, :]
        outlet_ix = 0
        for i in range(len(x_coords)):
            if x_coords[i] >= 0:
                outlet_ix = i
                break
        return outlet_ix

    def get_outlet_indices2(self, ix):
        # Get the outlet first index, cutremente.
        msh = self.msh[ix]
        x_coords = msh.cells_x_coords[0, :]
        outlet_ix = 0
        for i in range(len(x_coords)):
            if x_coords[i] >= 0:
                outlet_ix = i
                break
        return outlet_ix

    def plot_meshes(self):

        stop = True

    def plot_quivers(self, a_selected_mesh=0):
        pass

    def plot_phi_vs_x_mesh_size_influence_study(self):
        # Plot reference solution.
        x_reference_solution = self.reference_solution['x']
        y_reference_solution = None

        rho_gamma = self.usr[0].get('rho')
        if rho_gamma == 1e1:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e1', None)
        elif rho_gamma == 1e3:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e3', None)
        elif rho_gamma == 1e6:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e6', None)

        fig = go.Figure(data=go.Scatter(x=x_reference_solution, y=y_reference_solution, line=dict(color='firebrick'),
                                        name='Reference solution'))

        # Plot numerical solution.
        i = 0
        for i in range(3):
            msh = self.msh[i]
            field = self.fields[i]
            x = msh.cells_x_coords[0, :].tolist()[self.get_outlet_indices2(i)::]
            y = field.Phi[0, :].tolist()[self.get_outlet_indices2(i)::]
            fig.add_scatter(x=x, y=y, name='Numerical solution (mesh ' + str(i) + ')')
            i = i + 1
            # , line = dict(color='royalblue')

        # Figure config.
        separator = '-'
        space_m = '\:'
        # Text display hack in case of Smith-Hutton problem.
        problem = self.usr[0].get('problem', '')
        if problem == 'Smith-Hutton':
            problem = 'Smith \: Hutton'
        title = '\\phi' + space_m + 'value' + space_m + '@' + space_m + 'outlet' + separator + problem + \
                separator + '\\rho / \\Gamma = ' + str(rho_gamma)
        title = '$' + title + '$'
        fig.update_layout(title={
            'text': title,
            'x': 0.05,
            'xanchor': 'left',
            'yanchor': 'top'},
            xaxis_title="x (m)",
            yaxis_title="$\\phi$",
            font={'size': 15})

        fig.show()

    def plot_phi_vs_x_scheme_study(self):
        # Plot reference solution.
        x_reference_solution = self.reference_solution['x']
        y_reference_solution = None

        rho_gamma = self.usr[0].get('rho')
        if rho_gamma == 1e1:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e1', None)
        elif rho_gamma == 1e3:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e3', None)
        elif rho_gamma == 1e6:
            y_reference_solution = self.reference_solution.get('rho_gamma_1e6', None)

        fig = go.Figure(data=go.Scatter(x=x_reference_solution, y=y_reference_solution, line=dict(color='firebrick'),
                                        name='Reference solution'))

        # Plot numerical solution.
        n = len(self.fields)
        i = 0
        for i in range(n):
            msh = self.msh[i]
            field = self.fields[i]
            x = msh.cells_x_coords[0, :].tolist()[self.get_outlet_indices2(i)::]
            y = field.Phi[0, :].tolist()[self.get_outlet_indices2(i)::]
            fig.add_scatter(x=x, y=y, name='Numerical solution (' + self.usr[i]['numerical_scheme'] + ')')
            i = i + 1
            # , line = dict(color='royalblue')

        # Figure config.
        separator = '-'
        space_m = '\:'
        # Text display hack in case of Smith-Hutton problem.
        problem = self.usr[0].get('problem', '')
        if problem == 'Smith-Hutton':
            problem = 'Smith \: Hutton'
        title = '\\phi' + space_m + 'value' + space_m + '@' + space_m + 'outlet' + separator + problem + \
                separator + '\\rho / \\Gamma = ' + str(rho_gamma)
        title = '$' + title + '$'
        fig.update_layout(title={
            'text': title,
            'x': 0.05,
            'xanchor': 'left',
            'yanchor': 'top'},
            xaxis_title="x (m)",
            yaxis_title="$\\phi$",
            font={'size': 15})

        fig.show()
