from Mesher import Mesher
from Fields import Fields
from DiscretCoeffs import DiscretCoeffs
from SolverGaussSeidelPhi import SolverGaussSeidelPhi
from Plotter import Plotter
from CFDUtils import timer, tic, toc, load_reference_solution
import sys
import numpy as np


"""
MAIN script

Main code for the 2nd assignment: Convection-Diffusion.
This code runs solves the Diagonal flow or the Smith-Hutton problem, 
taking into account the user input throught user_input dict.

There is a separate script called NumericalStudyGenerator, which does the
mesh convergence, mesh size influence, and numerical scheme influence studies,
for the Smith-Hutton case only.

Physical units in International System.
"""

# User inputs (international units, unless otherwise stated).
user_input = dict()
user_input['problem'] = 'Diagonal flow'  # 'Smith-Hutton' | 'Diagonal flow'
user_input['alpha_diagonal_flow'] = 45
user_input['v0_diagonal_flow'] = 10
user_input['numerical_scheme'] = 'UDS'  # 'CDS' | 'UDS' | 'EDS'
user_input['rho'] = 1e9
user_input['gamma'] = 1
user_input['rho_gamma'] = user_input['rho'] / user_input['gamma']
user_input['cv_y'] = 19
user_input['cv_x'] = 2 * user_input['cv_y']  # For Smith-Hutton problem, to keep the mesh constant, as the domain is twice as long as high.
user_input['epsilon'] = 1e-6  # Convergence criteria.
user_input['output_files_name'] = 'CaseStudy_'

tic('Preprocessing')
# Build the mesh.
mesh = Mesher(user_input)
mesh.build_mesh()

# Initialize fields.
fields = Fields(user_input, mesh)
fields.build_velocity()
fields.build_phi()
fields.build_convection_strength()
fields.build_diffusion_strength()
fields.build_peclet_number()

# Initialize the discretization coefficients: aP, aE, aW, aN, aS.
coeffs = DiscretCoeffs(user_input, mesh, fields)
coeffs.build_coeffs()
toc('Preprocessing')

# Solve.
tic('Solving')
solver = SolverGaussSeidelPhi(user_input, mesh, fields, coeffs)
solver.solve_phi()
toc('Solving')

# Check solution.
if solver.solution_divergence:
    print('Solution has diverged. Try with a more stable numerical scheme. Aborting.')
    sys.exit()

# Reference solution.
reference_solution = load_reference_solution()

plotter =  Plotter(user_input, mesh, fields, reference_solution)
plotter.plot_mesh()
plotter.plot_quiver()
plotter.plot_phi_vs_x()
plotter.plot_phi_contours()
toc('Postprocessing')

# print('Elapsed time:', timer)
# print('abs(Psi - Psi_ana):', np.sum(fields.Psi) - np.sum(fields.Psi_analytical))
