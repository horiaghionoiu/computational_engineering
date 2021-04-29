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
        print('----------- START solving -----------')
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
        print('----------- END solving -----------')
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