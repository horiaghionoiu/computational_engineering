from PhysicProp import PhysicProp
from Mesher import Mesher
from Fields import Fields
from DiscretCoeffs import DiscretCoeffs
from SolverGaussSeidelPsi import SolverGaussSeidelPsi
from SolverPhysicalQuantities import SolverPhysicalQuantities
from Plotter import Plotter
from SolverAnalyticalPsi import SolverAnalyticalPsi
from CFDUtils import timer, tic, toc
import numpy as np

"""
Main code for the 1st assignment: Potential flow around a static cylinder.
Physical units in International System.
"""

# Simulation settings.
epsilon = 1e-6  # Convergence criteria.
phys_prop = PhysicProp('air')

tic('Preprocessing')
# Build the mesh.
mesh = Mesher()
mesh.set_control_volumes(a_cv_y=78, a_cv_x=158)
mesh.set_cylinder(a_R=1)
output_files_name = 'CaseStudy_'
n_lengths = 2
mesh.set_domain_size(a_n_Radius_west=-10*n_lengths,
                     a_n_Radius_east=10*n_lengths,
                     a_n_Radius_north=5*n_lengths,
                     a_n_Radius_south=-5*n_lengths)
mesh.build_mesh()

# Initialize velocity, pressure, temperature, density, and stream function fields.
fields = Fields(mesh, phys_prop)
fields.build_v()
fields.build_P()
fields.build_T()
fields.build_rho()
fields.build_Psi()

# Initialize stream function coefficients: aP, aE, aW, aN, aS.
coeffs = DiscretCoeffs(mesh, phys_prop, fields)
coeffs.build_coeffs()
toc('Preprocessing')

# Solve for the stream function.
tic('Solving')
solver = SolverGaussSeidelPsi(epsilon, mesh, fields, coeffs)
solver.solve_psi()
toc('Solving')

# Postprocess by solving the interesting flow physical variables.
tic('Postprocessing')
postprocessor = SolverPhysicalQuantities(mesh, fields, coeffs, phys_prop)
postprocessor.solve_velocity()
postprocessor.solve_temperature()
postprocessor.solve_pressure()
postprocessor.solve_rho()

# Analytical solution.
analytical_solution = SolverAnalyticalPsi(mesh, fields, phys_prop)
analytical_solution.solve_psi()

plotter = Plotter(mesh, fields, output_files_name)
plotter.plot_mesh()
# plotter.plot_stream_lines()
# plotter.plot_pressure(p_cp=False, p_relative=True)
# plotter.plot_temperature(p_relative=True)
# plotter.plot_rho(p_relative=True)
# plotter.plot_numerical_vs_analytical_contour()
# plotter.plot_error_vs_iterations(solver)
toc('Postprocessing')

print('Elapsed time:', timer)
print('abs(Psi - Psi_ana):', np.sum(fields.Psi) - np.sum(fields.Psi_analytical))
