import time

files_output_path = 'C:\\workcopy_computational_engineering\\docs\\'

timer = {'Preprocessing': 0,
         'Solving': 0,
         'Postprocessing': 0}
start_preprocessing = None
end_preprocessing = None
start_solving = None
end_solving = None
start_postprocessing = None
end_postprocessing = None

def tic(p_operation=None):
    """
    Start the timer.
    :param p_operation: 'Preprocessing'/'Solving'/'Postprocessing'
    """
    if p_operation is None:
        print('CFDUtils.tic: operation not specified.')
        return

    if p_operation == 'Preprocessing':
        print('----------- START preprocessing -----------')
        global start_preprocessing
        start_preprocessing = time.time()
    elif p_operation == 'Solving':
        # print('----------- START solving -----------')
        global start_solving
        start_solving = time.time()
    elif p_operation == 'Postprocessing':
        print('----------- START postprocessing -----------')
        global start_postprocessing
        start_postprocessing = time.time()
    else:
        print('CFDUtils.tic: wrong operation.')
        return

def toc(p_operation=None):
    """
    Stop the timer.
    :param p_operation: 'Preprocessing'/'Solving'/'Postprocessing'
    """
    if p_operation is None:
        print('CFDUtils.toc: operation not specified.')
        return

    global timer
    if p_operation == 'Preprocessing':
        print('----------- END preprocessing -----------')
        global end_preprocessing
        end_preprocessing = time.time()
        timer['Preprocessing'] = end_preprocessing - start_preprocessing
    elif p_operation == 'Solving':
        # print('----------- END solving -----------')
        global end_solving
        end_solving = time.time()
        timer['Solving'] = end_solving - start_solving
    elif p_operation == 'Postprocessing':
        print('----------- END postprocessing -----------')
        global end_postprocessing
        end_postprocessing = time.time()
        timer['Postprocessing'] = end_postprocessing - start_postprocessing
    else:
        print('CFDUtils.tic: wrong operation.')
        return

def load_reference_solution():
    reference_solution = dict()

    reference_solution['x'] =  [i/10 for i in range(11)]
    reference_solution['rho_gamma_1e1'] = [1.989, 1.402, 1.146, 0.946, 0.775, 0.621, 0.480, 0.349, 0.227, 0.111, 0]
    reference_solution['rho_gamma_1e3'] = [2, 1.9990, 1.9997, 1.9850, 1.8410, 0.9510, 0.1540, 0.0010, 0, 0, 0]
    reference_solution['rho_gamma_1e6'] = [2, 2, 2, 1.999, 1.964, 1, 0.036, 0.001, 0, 0, 0]

    return reference_solution