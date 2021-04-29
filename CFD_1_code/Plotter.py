# Imports
import plotly.express as px
import plotly.figure_factory as ff
import plotly.graph_objects as go
import plotly.io as pio
import plotly.express as px

import numpy as np
import pandas as pd
from scipy.interpolate import griddata

from Mesher import Mesher
from Fields import Fields
from SolverGaussSeidelPsi import SolverGaussSeidelPsi
from CFDUtils import files_output_path

"""
Class Plotter.
- In charge of: plotting interesting quantities.
"""


class Plotter:
    def __init__(self, a_the_mesh:Mesher = None, a_the_fields:Fields = None, a_output_files_name=''):
        self.msh = a_the_mesh
        self.fields = a_the_fields
        self.output_files_name = a_output_files_name

    @staticmethod
    def map_fluid_nodes(a_mesh:Mesher, nodes_x:np.array, nodes_y:np.array):

        fluid_nodes = np.ones((nodes_x.shape[0],nodes_y.shape[1]))

        for row in range(nodes_x.shape[0]):
            for col in range(nodes_x.shape[1]):
                x = nodes_x[row,col]
                y = nodes_y[row,col]
                distance_from_node_to_origin = np.sqrt(x**2 + y**2)

                if distance_from_node_to_origin < a_mesh.cylinder_R:
                    fluid_nodes[row,col] = 0.000000000000001
                # print('(',row,',',col,')')
        return fluid_nodes.copy()

    def plot_mesh(self):
        # Flatten the nodes matrixes.
        x = np.concatenate(self.msh.the_nodes_x_coords).ravel().tolist()
        y = np.concatenate(self.msh.the_nodes_y_coords).ravel().tolist()

        # Get the boundary nodes.
        x_boundary_nodes = list(self.msh.the_nodes_x_coords[0,:]) + \
                           list(self.msh.the_nodes_x_coords[-1,:]) + \
                           list(self.msh.the_nodes_x_coords[:,0]) + \
                           list(self.msh.the_nodes_x_coords[:,-1])
        y_boundary_nodes = list(self.msh.the_nodes_y_coords[0,:]) + \
                           list(self.msh.the_nodes_y_coords[-1,:]) + \
                           list(self.msh.the_nodes_y_coords[:,0]) + \
                           list(self.msh.the_nodes_y_coords[:,-1])
        # Plot the nodes.
        fig = go.Figure(data=go.Scatter(x=x, y=y, mode='markers', name='Internal nodes'))
        fig.add_scatter(x=x_boundary_nodes, y=y_boundary_nodes, mode='markers', name='Boundary nodes')

        # Plot the cylinder.
        theta = np.linspace(0, 2 * np.pi, 100)
        fig.add_scatter(x=self.msh.cylinder_R*np.cos(theta),y=self.msh.cylinder_R*np.sin(theta), fill='toself', fillcolor='violet',line_color='violet', name='Cylinder')
        # fig.layout.update({'title': 'Mesh'})
        fig.update_layout( title={
                                'text': "Mesh",
                                # 'x':0.5,
                                'xanchor': 'center',
                                'yanchor': 'top'},
                          xaxis_title="x (m)",
                          yaxis_title="y (m)",
                          font={'size':28})
        fig.show()
        pio.write_image(fig, files_output_path + self.output_files_name + 'Mesh.png', width=2864, height=1356)

    def plot_stream_lines(self):
        # Prepare the mesh and vx,vy velocity values for the streamlines plot function.
        nodes = self.msh.M * self.msh.N


        # Rearrange the nodes x,y coords in a single array, and flatten the vx,vy velocities.
        the_points = np.zeros((nodes, 2))
        the_vx_values = np.zeros((nodes,))
        the_vy_values = np.zeros((nodes,))
        i = 0
        for row in range(self.msh.M):
            for col in range(self.msh.N):
                the_points[i,0] = self.msh.the_nodes_x_coords[row,col]
                the_points[i,1] = self.msh.the_nodes_y_coords[row,col]
                the_vx_values[i] = self.fields.vx[row,col]
                the_vy_values[i] = self.fields.vy[row,col]
                i = i + 1

        x_ = np.linspace(self.msh.x0, self.msh.xf, 1000)
        y_ = np.linspace(self.msh.y0, self.msh.yf, 1000)
        x, y = np.meshgrid(x_, y_)
        VX = griddata(the_points, the_vx_values, (x, y), method='cubic')
        VY = griddata(the_points, the_vy_values, (x, y), method='cubic')

        # Last step, the original solid nodes, should remain solid in the mapped VX,VY velocities
        # for the streamlines plot.
        FLUID_NODES = self.map_fluid_nodes(self.msh, x, y)
        VX = VX * FLUID_NODES
        VY = VY * FLUID_NODES
        fig = ff.create_streamline(x_,y_, VX, VY, arrow_scale=.1, density=2,name='Streamlines')
        theta = np.linspace(0, 2 * np.pi, 100)
        fig.add_scatter(x=self.msh.cylinder_R*np.cos(theta),y=self.msh.cylinder_R*np.sin(theta), fill='toself', fillcolor='violet', line_color='violet', name='Cylinder')
        fig.layout.update({'title': 'Streamlines'},
                          font={'size':28})
        fig.show()
        pio.write_image(fig, files_output_path + self.output_files_name + 'Streamlines.png', width=2864, height=1356)

    def plot_pressure(self, p_cp=False, p_relative=False):
        if p_cp:
            # Plot the pressure coefficient.
            field_to_plot = self.fields.Cp
            title = '-'
        else:
            # Plot the pressure.
            if p_relative:
                field_to_plot = self.fields.P - self.fields.P0
            else:
                field_to_plot = self.fields.P
            title = 'P-P0 (Pa)'

        fig = go.Figure(
            layout=go.Layout(
                title=go.layout.Title(text="Pressure plot")
            ))
        fig.add_trace(go.Contour(
            z=field_to_plot,
            colorscale="aggrnyl",
            colorbar=dict(
                title=title
            )
        ))

        # fig.update_layout()
        fig.update_xaxes(title_text="VCx")
        fig.update_yaxes(title_text="VCy")
        fig.update_layout(font={'size':28})

        fig.show()
        pio.write_image(fig, files_output_path + self.output_files_name + 'Pressure.png', width=2864, height=1356)

    def plot_temperature(self, p_relative=False):
        # ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
        #  'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
        #  'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
        #  'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
        #  'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
        #  'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
        #  'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
        #  'orrd', 'oryel', 'oxy', 'peach', 'phase', 'picnic', 'pinkyl',
        #  'piyg', 'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn',
        #  'puor', 'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu',
        #  'rdgy', 'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar',
        #  'spectral', 'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn',
        #  'tealrose', 'tempo', 'temps', 'thermal', 'tropic', 'turbid',
        #  'turbo', 'twilight', 'viridis', 'ylgn', 'ylgnbu', 'ylorbr',
        #  'ylorrd']

        if p_relative:
            field_to_plot = self.fields.T - self.fields.T0
            title = 'T-T0 (K)'
        else:
            field_to_plot = self.fields.T
            title = 'K'

        fig = go.Figure(
            layout=go.Layout(
                title=go.layout.Title(text="Temperature plot")
            ))
        fig.add_trace(go.Contour(
            z=field_to_plot,
            colorscale="aggrnyl",
            colorbar=dict(
                title=title
            )
        ))

        # fig.update_layout()
        fig.update_xaxes(title_text="VCx")
        fig.update_yaxes(title_text="VCy")
        fig.update_layout(font={'size':28})

        fig.show()
        pio.write_image(fig, files_output_path + self.output_files_name + 'Temperature.png', width=2864, height=1356)

    def plot_rho(self, p_relative=False):

        if p_relative:
            the_field = self.fields.rho - self.fields.rho0
        else:
            the_field = self.fields.rho

        fig = go.Figure(
            layout=go.Layout(
                title=go.layout.Title(text="Density plot")
            ))
        fig.add_trace(go.Contour(
            z=the_field,
            colorscale="aggrnyl",
            colorbar=dict(
                title='rho-rho0 (kg/m3)'
            )
        ))

        # fig.update_layout()
        fig.update_xaxes(title_text="VCx")
        fig.update_yaxes(title_text="VCy")
        fig.update_layout(font={'size':28})

        fig.show()
        pio.write_image(fig, files_output_path + self.output_files_name + 'Density.png', width=2864, height=1356)

    def plot_error_vs_iterations(self, a_solver: SolverGaussSeidelPsi):
        df = pd.DataFrame(data={'Iteration':a_solver.solution_iterations, 'Error':a_solver.solution_error})
        fig = px.line(df, x="Iteration", y="Error", title='Solution convergence')
        fig.update_xaxes(title_text="Iteration")
        fig.update_yaxes(title_text="Error")
        fig.update_layout(font={'size':28})

        fig.show()
        pio.write_image(fig, files_output_path + self.output_files_name + 'Error_vs_Iteration.png', width=2864, height=1356)


    def plot_numerical_vs_analytical_contour(self):
        fig = go.Figure(
            layout=go.Layout(
                title=go.layout.Title(text="Numerical vs Analytical Stream function (Psi) contour plot")
            ))
        fig.add_trace(go.Contour(
            z=abs(self.fields.Psi - self.fields.Psi_analytical),
            colorscale="aggrnyl",
            colorbar=dict(
                title='abs(Psi-Psi_ana)'
            )
        ))

        # fig.update_layout()
        fig.update_xaxes(title_text="VCx")
        fig.update_yaxes(title_text="VCy")
        fig.update_layout(font={'size':28})

        fig.show()
        pio.write_image(fig, files_output_path + self.output_files_name + 'Numerical_vs_Analytical_Psi.png', width=2864, height=1356)