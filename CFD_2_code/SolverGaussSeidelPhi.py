# Imports
from Mesher import Mesher
from Fields import Fields
from DiscretCoeffs import DiscretCoeffs
import numpy as np

"""
Class SolverGaussSeidelPhi.
- In charge of: solving the flow equations.
Physical units in International System.
"""


class SolverGaussSeidelPhi:
    def __init__(self, a_user_input: dict, a_the_mesh:Mesher = None, a_the_fields:Fields = None,
                 a_the_coeffs:DiscretCoeffs = None):
        self.usr = a_user_input  # The user inputs.
        self.msh = a_the_mesh  # The mesh.
        self.fields = a_the_fields  # The fields.
        self.coeffs = a_the_coeffs  # The coeffs.
        self.solution_error = []
        self.solution_iterations = []
        self.solution_divergence = False
        self.timer = None
    
    def solve_phi(self, p_verbose_level=1):
        # Fill initial guess Phi with the boundary conditions already applied in Phi field.
        phi0 = self.fields.Phi.copy()
        # Fill initial guess Phi internal nodes.
        phi0_initial_guess = 1
        for row in range(1,self.msh.cv_y):
            for col in range(1,self.msh.cv_x):
                phi0[row,col] = phi0_initial_guess

        # Solve for self.fields.Phi.
        error = np.zeros((self.msh.M, self.msh.N))
        iterations_counter = 0
        all_errors_below_tolerance = False
        epsilon = self.usr.get('epsilon', 1e-6)

        while not all_errors_below_tolerance:
            # Protection.
            if len(self.solution_error) > 1:
                last_error = self.solution_error[-1]
                if last_error >= 1e6:
                    self.solution_divergence = True
                    self.log_error(iterations_counter, 999, p_last_iteration=True,
                                   p_solution_diverged=True, p_verbose_level=p_verbose_level)
                    return

            # Compute phi (excluding boundary nodes).
            self.compute_phi(phi0)

            # Compute error (excluding boundary nodes).
            self.compute_error(error, phi0)

            total_iteration_error = np.sum(error)
            # Logging.
            self.log_error(iterations_counter, total_iteration_error, p_verbose_level=p_verbose_level)

            # The sum of all errors should be less than n_nodes * epsilon.
            if total_iteration_error > error.size * epsilon:
                phi0 = self.fields.Phi.copy()
                iterations_counter = iterations_counter + 1

                # Instability check. If the solution diverges, then abort before overflow error.
                if len(self.solution_error) > 3:
                    last_error = self.solution_error[-1]
                    last_but_one_error = self.solution_error[-2]
                    last_but_two_error = self.solution_error[-3]
                    if last_error >= 2 * last_but_one_error:
                        if last_but_one_error >= 2 * last_but_two_error:
                            self.solution_divergence = True
                            self.log_error(iterations_counter, total_iteration_error, p_last_iteration=True, p_solution_diverged=True, p_verbose_level=p_verbose_level)
                            return


            else:
                all_errors_below_tolerance = True
                self.log_error(iterations_counter, total_iteration_error, p_last_iteration=True, p_verbose_level=p_verbose_level)

        # Solution converged successfully. Apply boundary condition to outlet, d(Phi)/dy = 0.
        end = len(self.fields.Phi[0,:].tolist())
        for col in range(self.get_outlet_indices(),end):
            self.fields.Phi[0,col] = self.fields.Phi[1,col]

    def log_error(self, iterations_counter, total_iteration_error, p_last_iteration=False, p_solution_diverged=False, p_verbose_level=1):
        if p_verbose_level:
            if ((not iterations_counter % 100) or (iterations_counter == 0) or p_last_iteration) and not p_solution_diverged:
                print('(Iter: ', iterations_counter, ', Error: ', total_iteration_error, ')')
                if p_last_iteration:
                    print('Solution converged successfully!')
            if p_solution_diverged:
                print('(Iter: ', iterations_counter, ', Error: ', total_iteration_error, ')')
                print('Aborting due to solution divergence.')
        self.solution_iterations.append(iterations_counter)
        self.solution_error.append(total_iteration_error)

    def compute_error(self, error, phi0):
        for row in range(1, self.msh.cv_y):
            for col in range(1, self.msh.cv_x):
                error[row, col] = abs(self.fields.Phi[row, col] - phi0[row, col])

    def compute_phi(self, phi0):
        for row in range(1, self.msh.cv_y):
            for col in range(1, self.msh.cv_x):
                # try:
                self.fields.Phi[row, col] = (self.coeffs.aE[row, col] * phi0[row, col + 1] +
                                             self.coeffs.aW[row, col] * self.fields.Phi[row, col - 1] +
                                             self.coeffs.aN[row, col] * phi0[row + 1, col] +
                                             self.coeffs.aS[row, col] * self.fields.Phi[row - 1, col]) / \
                                             self.coeffs.aP[row, col]
                # except:
                #     self.solution_divergence = True
                #     self.log_error(999, 999, p_last_iteration=True,
                #                    p_solution_diverged=True,
                #                    p_verbose_level=0)
                #     return

    def get_outlet_indices(self):
        # DUPLICATED CODE :XXXX
        # Get the outlet first index, cutremente.
        x_coords = self.msh.cells_x_coords[0, :]
        outlet_ix = 0
        for i in range(len(x_coords)):
            if x_coords[i] >= 0:
                outlet_ix = i
                break
        return outlet_ix