# Imports
from Mesher import Mesher
from Fields import Fields
from DiscretCoeffs import DiscretCoeffs
from SolverGaussSeidelPhi import SolverGaussSeidelPhi
from Plotter import Plotter
from CFDUtils import load_reference_solution, tic, toc, timer

import plotly.graph_objects as go

"""
NumericalStudyGenerator

Script which does the mesh convergence, mesh size influence, 
and numerical scheme influence studies, for the Smith-Hutton case only.

Physical units in International System.
"""


# Utility methods.
def get_user_input_for_mesh_size_influence_study(a_rho_gamma_studied, a_meshes_studied, a_study_name):
    user_input = dict()
    user_input['problem'] = 'Smith-Hutton'  # 'Smith-Hutton' | 'Diagonal flow'
    user_input['alpha_diagonal_flow'] = 45
    user_input['v0_diagonal_flow'] = 10
    user_input['numerical_scheme'] = 'UDS'  # 'CDS' | 'UDS' | 'EDS'
    user_input['rho'] = a_rho_gamma_studied
    user_input['gamma'] = 1  # Let it be 1.
    user_input['cv_y'] = a_meshes_studied
    user_input['cv_x'] = [2*i for i in a_meshes_studied]  # For Smith-Hutton problem, to keep the mesh constant, as the domain is twice as long as high.
    user_input['epsilon'] = 1e-6  # Convergence criteria.
    user_input['output_files_name'] = a_study_name
    return user_input

def perform_mesh_size_influence_study(a_usr_input: dict):
    study_result = dict()
    print('Performing mesh size influence study --- START.')


    for rho_gamma in a_usr_input['rho']:
        print('rho_gamma = ' + str(rho_gamma))
        i = 0
        for mesh in a_usr_input['cv_y']:
            i = i+1
            print('mesh = ' + str(mesh))

            # Prepare.
            current_study = {'mesh': None, 'fields': None, 'coeffs': None, 'solver': None}
            current_study_name = 'rho_gamma' + str(rho_gamma) + '_' + 'mesh' + str(i)
            current_study_user_input = a_usr_input.copy()
            current_study_user_input['rho'] = rho_gamma
            current_study_user_input['cv_y'] = mesh
            current_study_user_input['cv_x'] = 2*mesh  # For Smith-Hutton problem, to keep the mesh constant, as the domain is twice as long as high.

            # Perform.
            current_study['mesh'] = Mesher(current_study_user_input)
            current_study['mesh'].build_mesh()

            current_study['fields'] = Fields(current_study_user_input,
                                             current_study['mesh'])
            current_study['fields'].build_velocity()
            current_study['fields'].build_phi()
            current_study['fields'].build_convection_strength()
            current_study['fields'].build_diffusion_strength()
            current_study['fields'].build_peclet_number()

            current_study['coeffs'] = DiscretCoeffs(current_study_user_input,
                                                    current_study['mesh'],
                                                    current_study['fields'])
            current_study['coeffs'].build_coeffs()

            current_study['solver'] = SolverGaussSeidelPhi(current_study_user_input,
                                                           current_study['mesh'],
                                                           current_study['fields'],
                                                           current_study['coeffs'])
            tic('Solving')
            current_study['solver'].solve_phi(p_verbose_level=0)
            toc('Solving')
            current_study['solver'].timer = timer.copy()

            if (current_study['solver'].solution_divergence):
                print(current_study_name + ' - Solution has diverged.')

            # Save the current study.
            study_result[current_study_name] = current_study

    print('Performing mesh size influence study ---   END.')
    return study_result

def plot_three_meshes(a_the_meshes: list):
    # a_the_meshes = [study1['rho_gamma10.0_mesh1']['mesh'],
    #                 study1['rho_gamma10.0_mesh2']['mesh'],
    #                 study1['rho_gamma10.0_mesh3']['mesh']]

    usr = study1['rho_gamma10.0_mesh1']['mesh'].usr

    msh = a_the_meshes[0]
    plotter = Plotter(usr, msh)
    plotter.plot_mesh()
    usr = study1['rho_gamma10.0_mesh2']['mesh'].usr
    msh = a_the_meshes[1]
    plotter = Plotter(usr, msh)
    plotter.plot_mesh()
    usr = study1['rho_gamma10.0_mesh3']['mesh'].usr
    msh = a_the_meshes[2]
    plotter = Plotter(usr, msh)
    plotter.plot_mesh()

def plot_three_quivers(a_the_meshes: list, a_the_fields: list):
    # a_the_meshes = [study1['rho_gamma10.0_mesh2']['mesh'],
    #                 study1['rho_gamma1000.0_mesh2']['mesh'],
    #                 study1['rho_gamma1000000.0_mesh2']['mesh']]

    # a_the_fields = [study1['rho_gamma10.0_mesh2']['fields'],
    #                 study1['rho_gamma1000.0_mesh2']['fields'],
    #                 study1['rho_gamma1000000.0_mesh2']['fields']]

    usr = study1['rho_gamma10.0_mesh2']['mesh'].usr

    for i in range(3):
        msh = a_the_meshes[i]
        field = a_the_fields[i]
        plotter = Plotter(usr, msh, field)
        plotter.plot_quiver()

def plot_3mesh_one_gamma_phi_vs_x(a_the_usr: list, a_the_meshes: list, a_the_fields: list):
    # a_the_usr = [study1['rho_gamma10.0_mesh1']['mesh'].usr,
    #              study1['rho_gamma10.0_mesh2']['mesh'].usr,
    #              study1['rho_gamma10.0_mesh3']['mesh'].usr]
    # a_the_meshes = [study1['rho_gamma10.0_mesh1']['mesh'],
    #                 study1['rho_gamma10.0_mesh2']['mesh'],
    #                 study1['rho_gamma10.0_mesh3']['mesh']]
    # a_the_fields = [study1['rho_gamma10.0_mesh1']['fields'],
    #                 study1['rho_gamma10.0_mesh2']['fields'],
    #                 study1['rho_gamma10.0_mesh3']['fields']]


    reference_solution = load_reference_solution()
    plotter = Plotter(a_the_usr, a_the_meshes, a_the_fields, reference_solution)
    plotter.plot_phi_vs_x_mesh_size_influence_study()

def get_user_input_for_mesh_convergence_study(a_rho_gamma_studied, a_meshes_studied, a_study_name):
    user_input = dict()
    user_input['problem'] = 'Smith-Hutton'  # 'Smith-Hutton' | 'Diagonal flow'
    user_input['alpha_diagonal_flow'] = 45
    user_input['v0_diagonal_flow'] = 10
    user_input['numerical_scheme'] = 'UDS'  # 'CDS' | 'UDS' | 'EDS'
    user_input['rho'] = a_rho_gamma_studied
    user_input['gamma'] = 1  # Let it be 1.
    user_input['cv_y'] = a_meshes_studied
    user_input['cv_x'] = [2*i for i in a_meshes_studied]  # For Smith-Hutton problem, to keep the mesh constant, as the domain is twice as long as high.
    user_input['epsilon'] = 1e-6  # Convergence criteria.
    user_input['output_files_name'] = a_study_name
    return user_input

def perform_mesh_convergence_study(a_usr_input: dict):
    study_result = dict()
    print('Performing mesh convergence study --- START.')


    for rho_gamma in a_usr_input['rho']:
        print('rho_gamma = ' + str(rho_gamma))
        i = 0
        for mesh in a_usr_input['cv_y']:
            i = i+1
            print('mesh = ' + str(i))

            # Prepare.
            current_study = {'mesh': None, 'fields': None, 'coeffs': None, 'solver': None}
            current_study_name = 'rho_gamma' + str(rho_gamma) + '_' + 'mesh' + str(i)
            current_study_user_input = a_usr_input.copy()
            current_study_user_input['rho'] = rho_gamma
            current_study_user_input['cv_y'] = mesh
            current_study_user_input['cv_x'] = 2*mesh  # For Smith-Hutton problem, to keep the mesh constant, as the domain is twice as long as high.
            current_study_user_input['number_of_nodes'] = mesh * (2 * mesh)

            # Perform.
            current_study['mesh'] = Mesher(current_study_user_input)
            current_study['mesh'].build_mesh()

            current_study['fields'] = Fields(current_study_user_input,
                                             current_study['mesh'])
            current_study['fields'].build_velocity()
            current_study['fields'].build_phi()
            current_study['fields'].build_convection_strength()
            current_study['fields'].build_diffusion_strength()
            current_study['fields'].build_peclet_number()

            current_study['coeffs'] = DiscretCoeffs(current_study_user_input,
                                                    current_study['mesh'],
                                                    current_study['fields'])
            current_study['coeffs'].build_coeffs()

            current_study['solver'] = SolverGaussSeidelPhi(current_study_user_input,
                                                           current_study['mesh'],
                                                           current_study['fields'],
                                                           current_study['coeffs'])
            tic('Solving')
            current_study['solver'].solve_phi(p_verbose_level=0)
            toc('Solving')
            current_study['solver'].timer = timer.copy()

            if (current_study['solver'].solution_divergence):
                print(current_study_name + ' - Solution has diverged.')

            # Save the current study.
            study_result[current_study_name] = current_study

    print('Performing mesh convergence study ---   END.')
    return study_result

def get_mesh_convergence_error(the_case, reference_solution):
    total_reference_solution_outlet_phi = sum(reference_solution)
    total_numerical_solution_outlet_phi = 0

    first_outlet_node_index = 0
    x_coords = the_case['mesh'].cells_x_coords[0, :]
    for i in range(len(x_coords)):
        if x_coords[i] >= 0:
            first_outlet_node_index = i
            break
    numerical_solution_outlet_phi = the_case['fields'].Phi[0,:].tolist()[first_outlet_node_index::]

    total_numerical_solution_outlet_phi = sum(numerical_solution_outlet_phi)

    return abs(total_reference_solution_outlet_phi - total_numerical_solution_outlet_phi)

def plot_mesh_convergence(a_the_cases: list, a_reference_solution, a_rho_gamma):
    # a_the_cases = [study2['rho_gamma10.0_mesh1'],
    #                 study2['rho_gamma10.0_mesh2'],
    #                 study2['rho_gamma10.0_mesh3'],
    #                 study2['rho_gamma10.0_mesh4'],
    #                 study2['rho_gamma10.0_mesh5'],
    #                 study2['rho_gamma10.0_mesh6'],
    #                 study2['rho_gamma10.0_mesh7'],
    #                 study2['rho_gamma10.0_mesh8'],
    #                 study2['rho_gamma10.0_mesh9']]

    # Get the data.
    n = len(a_the_cases)
    x = [0 for i in range(n)]
    y = [0 for i in range(n)]
    for i in range(n):
        number_of_nodes = a_the_cases[i]['mesh'].usr['number_of_nodes']
        error = get_mesh_convergence_error(a_the_cases[i], a_reference_solution)
        x[i] = number_of_nodes
        y[i] =  error

    # Plot the data.
    fig = go.Figure(data=go.Scatter(x=x, y=y, line=dict(color='firebrick'),
                                    name='Error vs Number of nodes'))
    separator = '-'
    space_m = '\:'
    # Text display hack in case of Smith-Hutton problem.
    problem = 'Smith \: Hutton'
    title = 'Mesh' + space_m +  'convergence'+ space_m + 'study' + separator + '\\rho / \\Gamma = ' + str(a_rho_gamma)
    title = '$' + title + '$'
    fig.update_layout( title={
                            'text': title,
                            'x':0.05,
                            'xanchor': 'left',
                            'yanchor': 'top'},
                      xaxis_title="Number of nodes (total)",
                      yaxis_title="Error",
                      font={'size':15})

    fig.show()

def get_user_input_for_numerical_scheme_study(a_rho_gamma_studied, a_convective_scheme, a_study_name):
    user_input = dict()
    user_input['problem'] = 'Smith-Hutton'  # 'Smith-Hutton' | 'Diagonal flow'
    user_input['alpha_diagonal_flow'] = 45
    user_input['v0_diagonal_flow'] = 10
    user_input['numerical_scheme'] = a_convective_scheme  # 'CDS' | 'UDS' | 'EDS'
    user_input['rho'] = a_rho_gamma_studied
    user_input['gamma'] = 1  # Let it be 1.
    user_input['cv_y'] = 21
    user_input['cv_x'] = 2 * user_input['cv_y'] # For Smith-Hutton problem, to keep the mesh constant, as the domain is twice as long as high.
    user_input['epsilon'] = 1e-6  # Convergence criteria.
    user_input['output_files_name'] = a_study_name
    return user_input

def perform_numerical_scheme_study(a_usr_input):
    study_result = dict()
    print('Performing numerical scheme study --- START.')

    for scheme in a_usr_input['numerical_scheme']:
        print('Numerical scheme = ' + scheme)
        for rho_gamma in a_usr_input['rho']:
            print('rho_gamma = ' + str(rho_gamma))
            # Prepare.
            current_study = {'mesh': None, 'fields': None, 'coeffs': None, 'solver': None}
            current_study_name = 'rho_gamma_' + str(rho_gamma) + '_scheme_' + scheme
            current_study_user_input = a_usr_input.copy()
            current_study_user_input['numerical_scheme'] = scheme
            current_study_user_input['rho'] = rho_gamma

            # Perform.
            current_study['mesh'] = Mesher(current_study_user_input)
            current_study['mesh'].build_mesh()

            current_study['fields'] = Fields(current_study_user_input,
                                             current_study['mesh'])
            current_study['fields'].build_velocity()
            current_study['fields'].build_phi()
            current_study['fields'].build_convection_strength()
            current_study['fields'].build_diffusion_strength()
            current_study['fields'].build_peclet_number()

            current_study['coeffs'] = DiscretCoeffs(current_study_user_input,
                                                    current_study['mesh'],
                                                    current_study['fields'])
            current_study['coeffs'].build_coeffs()

            current_study['solver'] = SolverGaussSeidelPhi(current_study_user_input,
                                                           current_study['mesh'],
                                                           current_study['fields'],
                                                           current_study['coeffs'])
            tic('Solving')
            current_study['solver'].solve_phi(p_verbose_level=0)
            toc('Solving')
            current_study['solver'].timer = timer.copy()

            if (current_study['solver'].solution_divergence):
                print(current_study_name + ' - Solution has diverged.')

            # Save the current study.
            study_result[current_study_name] = current_study

    print('Performing numerical scheme study --- END.')
    return study_result

def plot_numerical_scheme(a_schemes, a_reference_solution, a_rho_gamma):
    the_cases = []
    for scheme in a_schemes:
        the_key = 'rho_gamma_' + str(a_rho_gamma) + '_scheme_' + scheme
        the_cases.append(study3[the_key])

    usr = [elem.get('mesh').usr for elem in the_cases]
    msh = [elem.get('mesh') for elem in the_cases]
    fields = [elem.get('fields') for elem in the_cases]

    plotter = Plotter(usr, msh, fields, a_reference_solution)
    plotter.plot_phi_vs_x_scheme_study()

if __name__ == '__main__':
    tal =    __name__ == '__main__'
    no_tal = __name__ != '__main__'
    if no_tal:
        # ----------------------- Mesh size influence study.
        # For each rho/gamma, and for each each mesh.

        rho_gamma_studied = [1e1, 1e3, 1e6]
        meshes_studied = [10, 15, 20]  # Represents the cv_y. cv_x=2*cv_y, to achieve an uniform mesh.
        study_name = 'Mesh_size_influence_study'
        user_input = get_user_input_for_mesh_size_influence_study(rho_gamma_studied, meshes_studied, study_name)

        study1 = perform_mesh_size_influence_study(user_input)

        # Plot the meshes (3 figs).
        the_meshes = [study1['rho_gamma10.0_mesh1']['mesh'],
                      study1['rho_gamma10.0_mesh2']['mesh'],
                      study1['rho_gamma10.0_mesh3']['mesh']]
        plot_three_meshes(the_meshes)

        # Plot the quiver plots, for one mesh only, and for each rho_gamma (3 figs).
        # the_meshes = [study1['rho_gamma10.0_mesh2']['mesh'],
        #                 study1['rho_gamma1000.0_mesh2']['mesh'],
        #                 study1['rho_gamma1000000.0_mesh2']['mesh']]

        # the_fields = [study1['rho_gamma10.0_mesh2']['fields'],
        #                 study1['rho_gamma1000.0_mesh2']['fields'],
        #                 study1['rho_gamma1000000.0_mesh2']['fields']]
        # plot_three_quivers(the_meshes, the_fields)

        # Plot the phi @ outlet, for all meshes,and for each rho_gamma (3 figs).
        the_usr = [study1['rho_gamma10.0_mesh1']['mesh'].usr,
                   study1['rho_gamma10.0_mesh2']['mesh'].usr,
                   study1['rho_gamma10.0_mesh3']['mesh'].usr]
        the_meshes = [study1['rho_gamma10.0_mesh1']['mesh'],
                      study1['rho_gamma10.0_mesh2']['mesh'],
                      study1['rho_gamma10.0_mesh3']['mesh']]
        the_fields = [study1['rho_gamma10.0_mesh1']['fields'],
                      study1['rho_gamma10.0_mesh2']['fields'],
                      study1['rho_gamma10.0_mesh3']['fields']]
        plot_3mesh_one_gamma_phi_vs_x(the_usr, the_meshes, the_fields)

        the_usr = [study1['rho_gamma1000.0_mesh1']['mesh'].usr,
                   study1['rho_gamma1000.0_mesh2']['mesh'].usr,
                   study1['rho_gamma1000.0_mesh3']['mesh'].usr]
        the_meshes = [study1['rho_gamma1000.0_mesh1']['mesh'],
                      study1['rho_gamma1000.0_mesh2']['mesh'],
                      study1['rho_gamma1000.0_mesh3']['mesh']]
        the_fields = [study1['rho_gamma1000.0_mesh1']['fields'],
                      study1['rho_gamma1000.0_mesh2']['fields'],
                      study1['rho_gamma1000.0_mesh3']['fields']]
        plot_3mesh_one_gamma_phi_vs_x(the_usr, the_meshes, the_fields)

        the_usr = [study1['rho_gamma1000000.0_mesh1']['mesh'].usr,
                   study1['rho_gamma1000000.0_mesh2']['mesh'].usr,
                   study1['rho_gamma1000000.0_mesh3']['mesh'].usr]
        the_meshes = [study1['rho_gamma1000000.0_mesh1']['mesh'],
                      study1['rho_gamma1000000.0_mesh2']['mesh'],
                      study1['rho_gamma1000000.0_mesh3']['mesh']]
        the_fields = [study1['rho_gamma1000000.0_mesh1']['fields'],
                      study1['rho_gamma1000000.0_mesh2']['fields'],
                      study1['rho_gamma1000000.0_mesh3']['fields']]
        plot_3mesh_one_gamma_phi_vs_x(the_usr, the_meshes, the_fields)

        # Result of mesh size influence study.
        selected_mesh = 20
        # del rho_gamma_studied, the_usr, study1,the_fields, meshes_studied, the_meshes, study_name, user_input, selected_mesh
        # ----------------------- END

    if tal:
        # ----------------------- Mesh convergence study.
        rho_gamma_studied = [1e1, 1e3, 1e6] # Represents the cv_y. cv_x=2*cv_y, to achieve an uniform mesh.
        meshes_studied = [10, 15, 20, 25, 30, 35, 40, 45, 50]
        # meshes_studied = [5, 7, 9, 11, 13, 15, 17, 19, 20]
        reference_solution = load_reference_solution()
        study_name = 'Mesh_convergence_study'
        user_input = get_user_input_for_mesh_convergence_study(rho_gamma_studied, meshes_studied, study_name)

        study2 = perform_mesh_convergence_study(user_input)

        the_cases = [study2['rho_gamma10.0_mesh1'],
                     study2['rho_gamma10.0_mesh2'],
                     study2['rho_gamma10.0_mesh3'],
                     study2['rho_gamma10.0_mesh4'],
                     study2['rho_gamma10.0_mesh5'],
                     study2['rho_gamma10.0_mesh6'],
                     study2['rho_gamma10.0_mesh7'],
                     study2['rho_gamma10.0_mesh8'],
                     study2['rho_gamma10.0_mesh9']]
        plot_mesh_convergence(the_cases, reference_solution['rho_gamma_1e1'], 1e1)

        the_cases = [study2['rho_gamma1000.0_mesh1'],
                     study2['rho_gamma1000.0_mesh2'],
                     study2['rho_gamma1000.0_mesh3'],
                     study2['rho_gamma1000.0_mesh4'],
                     study2['rho_gamma1000.0_mesh5'],
                     study2['rho_gamma1000.0_mesh6'],
                     study2['rho_gamma1000.0_mesh7'],
                     study2['rho_gamma1000.0_mesh8'],
                     study2['rho_gamma1000.0_mesh9']]
        plot_mesh_convergence(the_cases, reference_solution['rho_gamma_1e3'], 1e3)

        the_cases = [study2['rho_gamma1000000.0_mesh1'],
                     study2['rho_gamma1000000.0_mesh2'],
                     study2['rho_gamma1000000.0_mesh3'],
                     study2['rho_gamma1000000.0_mesh4'],
                     study2['rho_gamma1000000.0_mesh5'],
                     study2['rho_gamma1000000.0_mesh6'],
                     study2['rho_gamma1000000.0_mesh7'],
                     study2['rho_gamma1000000.0_mesh8'],
                     study2['rho_gamma1000000.0_mesh9']]
        plot_mesh_convergence(the_cases, reference_solution['rho_gamma_1e6'], 1e6)

        # Result of mesh convergence study.
        selected_mesh_convergence = 21
        # del rho_gamma_studied, meshes_studied, study_name, user_input, study2, the_cases
        # ----------------------- END

    if no_tal:
        # ----------------------- Numerical scheme influence study.
        convective_scheme = ['UDS', 'CDS', 'EDS']
        rho_gamma_studied = [1e1, 1e3, 1e6]  # Represents the cv_y. cv_x=2*cv_y, to achieve an uniform mesh.
        reference_solution = load_reference_solution()
        study_name = 'Numerical_scheme_study'
        user_input = get_user_input_for_numerical_scheme_study(rho_gamma_studied, convective_scheme, study_name)

        study3 = perform_numerical_scheme_study(user_input)

        # Plot only the numerical schemes which converged.
        schemes = ['UDS', 'CDS']
        plot_numerical_scheme(schemes, reference_solution, 1e1)

        schemes = ['UDS', 'EDS']
        plot_numerical_scheme(schemes, reference_solution, 1e3)

        schemes = ['UDS']
        plot_numerical_scheme(schemes, reference_solution, 1e6)

        # Result of the numerical scheme study.
        selected_scheme = 'UDS'
        # ----------------------- END

    # ----------------------- Numerical results for a selected configuration, and for each gamma.
    pass
    # ----------------------- END
