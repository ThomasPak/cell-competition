import os
from sys import argv

from collections import Counter

def getSimulationDirectories(test_dir):
    """
    Assumes that there are only simulation directories under test_dir, no other
    directories.
    """
    return next(os.walk(test_dir))[1]

def getSimulationDirectoryKeys(sim_dir):
    """
    Assumes that sim_dir is of the format key1_val1_key2_val2_...keyN_valN.
    """
    tokens = sim_dir.split('_')
    assert len(tokens) % 2 == 0

    return [tokens[2 * i] for i in range(len(tokens) // 2)]

def extractSimulationDirectoryData(sim_dir, sim_data):
    """
    Assumes that sim_dir is of the format key1_val1_key2_val2_...keyN_valN.
    """
    tokens = sim_dir.split('_')
    assert len(tokens) % 2 == 0

    for i in range(len(tokens) // 2):
        sim_data[tokens[2 * i]] = tokens[2 * i + 1]

def extractVizancestorsData(test_dir, sim_dir, sim_data, results_dir="results_from_time_0"):
    """
    Assumes that test_dir/sim_dir/results_dir/results.vizancestors exists

    Updates sim_data with
    - initial_timestep
    - final_timestep
    - initial_cell_count
    - final_cell_count
    - last_alive_timestep
    - initial_ancestors
    - For every ancestor x
      - initial_cell_count_x
      - final_cell_count_x
      - last_alive_timestep_x
    """
    ancestors_file = test_dir + "/" + sim_dir + "/" + results_dir + "/results.vizancestors"
    assert os.path.isfile(ancestors_file)

    with open(ancestors_file) as fh:

        # Get timestep and cell count from first line
        line = next(fh)
        timestep, ancestors = line.split('\t')
        cell_count = len(ancestors.split())
        cell_count_per_ancestor = Counter(ancestors.split())

        # Get initial ancestors
        initial_ancestors = ' '.join(sorted(cell_count_per_ancestor.keys()))

        # Get initial timestep and cell count
        initial_timestep = timestep
        initial_cell_count = cell_count
        initial_cell_count_per_ancestor = cell_count_per_ancestor

        # Initialise last_alive_timestep
        last_alive_timestep = timestep
        last_alive_timestep_per_ancestor = {ancestor : timestep for ancestor in initial_ancestors.split()}

        # Loop over remaining lines
        for line in fh:

            # Get timestep and cell count
            timestep, ancestors = line.split('\t')
            cell_count = len(ancestors.split())
            cell_count_per_ancestor = Counter(ancestors.split())

            # If there are cells, set last_alive_timestep
            if cell_count > 0:
                last_alive_timestep = timestep

            # Counter object implicitly only includes keys with nonzero count
            for ancestor in cell_count_per_ancestor.keys():
                last_alive_timestep_per_ancestor[ancestor] = timestep

        # Set last_timestep
        final_timestep = timestep
        final_cell_count = cell_count
        final_cell_count_per_ancestor = cell_count_per_ancestor

    # Store data
    sim_data['initial_ancestors'] = initial_ancestors
    sim_data['initial_timestep'] = initial_timestep
    sim_data['final_timestep'] = final_timestep
    sim_data['initial_cell_count'] = str(initial_cell_count)
    sim_data['final_cell_count'] = str(final_cell_count)
    sim_data['last_alive_timestep'] = last_alive_timestep

    for ancestor in sim_data['initial_ancestors'].split():
        if ancestor in initial_cell_count_per_ancestor:
            sim_data['initial_cell_count_{}'.format(ancestor)] = str(initial_cell_count_per_ancestor[ancestor])
        else:
            raise Exception('This should not be reached')

        if ancestor in final_cell_count_per_ancestor:
            sim_data['final_cell_count_{}'.format(ancestor)] = str(final_cell_count_per_ancestor[ancestor])
        else:
            sim_data['final_cell_count_{}'.format(ancestor)] = str(0)

        if ancestor in last_alive_timestep_per_ancestor:
            sim_data['last_alive_timestep_{}'.format(ancestor)] = str(last_alive_timestep_per_ancestor[ancestor])
        else:
            raise Exception('This should not be reached')

def extractAncestorApoptoticsData(test_dir, sim_dir, sim_data, results_dir="results_from_time_0"):
    """
    Assumes that test_dir/sim_dir/results_dir/T2AncestorsIsApoptotics.dat exists

    Assumes that sim_data contains initial_ancestors, extracted from
    vizancestors data

    Updates sim_data with
    - num_apoptosis
    - num_extrusion
    - For every ancestor x
      - num_apoptosis_x
      - num_extrusion_x
    """
    is_apoptotic_file = test_dir + "/" + sim_dir + "/" + results_dir + "/T2AncestorsIsApoptotics.dat"
    assert os.path.isfile(is_apoptotic_file)

    num_apoptosis = 0
    num_extrusion = 0

    num_apoptosis_per_ancestor = {}
    num_extrusion_per_ancestor = {}
    with open(is_apoptotic_file) as fh:

        for line in fh:
            for word in line.split('\t'):

                if word == 'apoptosis':
                    num_apoptosis += 1
                elif word == 'extrusion':
                    num_extrusion += 1

            # Loop over ancestor, death_type pairs
            for ancestor, death_type in zip(line.split('\t')[2::2], line.split('\t')[3::2]):

                if death_type == 'apoptosis':
                    if ancestor in num_apoptosis_per_ancestor:
                        num_apoptosis_per_ancestor[ancestor] += 1
                    else:
                        num_apoptosis_per_ancestor[ancestor] = 1
                elif death_type == 'extrusion':
                    if ancestor in num_extrusion_per_ancestor:
                        num_extrusion_per_ancestor[ancestor] += 1
                    else:
                        num_extrusion_per_ancestor[ancestor] = 1
                else:
                    raise Exception('Death type "{}" not recognised'.format(death_type))

    sim_data['num_apoptosis'] = str(num_apoptosis)
    sim_data['num_extrusion'] = str(num_extrusion)

    for ancestor in sim_data['initial_ancestors'].split():
        if ancestor in num_apoptosis_per_ancestor:
            sim_data['num_apoptosis_{}'.format(ancestor)] = str(num_apoptosis_per_ancestor[ancestor])
        else:
            sim_data['num_apoptosis_{}'.format(ancestor)] = str(0)

        if ancestor in num_extrusion_per_ancestor:
            sim_data['num_extrusion_{}'.format(ancestor)] = str(num_extrusion_per_ancestor[ancestor])
        else:
            sim_data['num_extrusion_{}'.format(ancestor)] = str(0)

def getParallelJoblogKeys(test_dir):
    """
    Assumes that test_dir/test_dir.joblog exists
    Assumes that the simulation number is the Seq field minus one.

    Returns keys.
    """
    joblog_file = test_dir + "/" + os.path.basename(test_dir) + ".joblog"
    assert os.path.isfile(joblog_file)

    with open(joblog_file) as fh:

        # Read keys from header
        header = next(fh)
        return header.split()

def readParallelJoblog(test_dir):
    """
    Assumes that test_dir/test_dir.joblog exists
    Assumes that the simulation number is the Seq field minus one.

    Returns mapping from simulation number to a dictionary containing the
    simulation's joblog data.
    """
    joblog_file = test_dir + "/" + os.path.basename(test_dir) + ".joblog"
    assert os.path.isfile(joblog_file)

    joblog_dicts = []

    with open(joblog_file) as fh:

        # Read keys from header
        header = next(fh)
        keys = header.split()

        for line in fh:
            # Extract values
            vals = line.split('\t')

            # Create and insert dictionary
            joblog_dicts.append(dict(zip(keys, vals)))

    return {int(d['Seq']) - 1 : d for d in joblog_dicts}

def extractData(data, sim_dir, sim_data):

    # Get simulation number fom sim_dir
    tokens = sim_dir.split('_')
    d = {tokens[i * 2] : tokens[i * 2 + 1] for i in range(len(tokens) // 2)}
    sim_number = int(d['sim'])

    # Update sim_data with appropriate data dict
    sim_data.update(data[sim_number])

def getSacctOutputKeys(test_dir):
    """
    Assumes that test_dir/test_dir.sacct.txt exists
    Assumes that the simulation number is the part of the JobID after the
    underscore

    Returns keys.
    """
    sacct_file = test_dir + "/" + os.path.basename(test_dir) + ".sacct.txt"
    assert os.path.isfile(sacct_file)

    with open(sacct_file) as fh:

        # Read keys from header
        header = next(fh)
        return header.split()

def readSacctOutput(test_dir):
    """
    Assumes that test_dir/test_dir.sacct.txt exists
    Assumes that the simulation number is the part of the JobID after the
    underscore

    Returns mapping from simulation number to a dictionary containing the
    simulation's sacct output data.
    """
    sacct_file = test_dir + "/" + os.path.basename(test_dir) + ".sacct.txt"
    assert os.path.isfile(sacct_file)

    sacct_dicts = []

    with open(sacct_file) as fh:

        # Read keys from header
        header = next(fh)
        keys = header.split()

        # Skip separating line
        dummy = next(fh)

        for line in fh:
            # Extract values
            vals = line.split()

            # Create and insert dictionary
            sacct_dicts.append(dict(zip(keys, vals)))

    return {int(d['JobID'].split('_')[1]) : d for d in sacct_dicts}

def getParamKeys(test_dir):
    """
    Assumes that test_dir/test_dir.parameters.csv exists
    Assumes that the simulation number is given by simulation-id=N

    Returns keys.
    """
    param_file = test_dir + "/" + os.path.basename(test_dir) + ".parameters.csv"
    assert os.path.isfile(param_file)

    with open(param_file) as fh:

        # Read keys from header
        header = next(fh)
        return header.strip().split(',')

def readParamOutput(test_dir):
    """
    Assumes that test_dir/test_dir.parameters.csv exists
    Assumes that the simulation number is given by simulation-id=N

    Returns mapping from simulation number to a dictionary containing the
    simulation's param output data.
    """
    param_file = test_dir + "/" + os.path.basename(test_dir) + ".parameters.csv"
    assert os.path.isfile(param_file)

    param_dicts = []

    with open(param_file) as fh:

        # Read keys from header
        header = next(fh)
        keys = getParamKeys(test_dir)

        for line in fh:
            # Extract values
            vals = [item.split('=')[1] for item in line.strip().split(',')]

            # Create and insert dictionary
            param_dicts.append(dict(zip(keys, vals)))

    return {int(d['simulation-id']) : d for d in param_dicts}

def printData(sim_datas, keys):

    # print header
    print(','.join(keys))

    for sim_data in sim_datas:
        vals = [sim_data[key].strip() for key in keys]
        print(','.join(vals))

if __name__ == '__main__':

    # Help message
    if len(argv) < 2:
        print("Usage: {} TEST_DIR [KEY_LIST] [SORT_BY_KEY]\n"
                "TEST_DIR       is the test directory\n"
                "KEY_LIST       is a comma-separated list of keys to output\n"
                "               By default, all keys are outputted\n"
                "SORT_BY_KEY    is the key by which to sort the simulations\n"
                "               The field to sort by must be convertible to a float\n"
                "               By default, the simulations are sorted by 'sim'".format(argv[0]))
        exit(1)

    # Set test_dir
    test_dir = argv[1].strip().rstrip('/')
    assert os.path.isdir(test_dir)

    # If there is a GNU parallel joblog file
    has_joblog_file =  os.path.isfile(test_dir + '/' + os.path.basename(test_dir) + '.joblog')
    if has_joblog_file:
        joblog_data = readParallelJoblog(test_dir)

    # If there is a sacct file
    has_sacct_file =  os.path.isfile(test_dir + '/' + os.path.basename(test_dir) + '.sacct.txt')
    if has_sacct_file:
        sacct_data = readSacctOutput(test_dir)

    # If there is a parameters file
    has_param_file =  os.path.isfile(test_dir + '/' + os.path.basename(test_dir) + '.parameters.csv')
    if has_param_file:
        param_data = readParamOutput(test_dir)

    #  Loop over simulation directories
    sim_datas = []
    sim_dirs = getSimulationDirectories(test_dir)
    for sim_dir in sim_dirs:

        sim_data = {}

        # Extract data from simulation directory name
        extractSimulationDirectoryData(sim_dir, sim_data)

        # Extract data from results.vizancestors
        extractVizancestorsData(test_dir, sim_dir, sim_data)

        # Extract data from T2AncestorsIsApoptotics.dat
        extractAncestorApoptoticsData(test_dir, sim_dir, sim_data)

        # Extract data from parallel joblog data
        if has_joblog_file:
            extractData(joblog_data, sim_dir, sim_data)

        # Extract data from parallel sacct data
        if has_sacct_file:
            extractData(sacct_data, sim_dir, sim_data)

        # Extract data from parameters data
        if has_param_file:
            extractData(param_data, sim_dir, sim_data)

        sim_datas.append(sim_data)

    # Get keys from simulation directory name
    sim_dir_keys = getSimulationDirectoryKeys(sim_dirs[0])

    # Get keys from parallel joblog file
    if has_joblog_file:
        joblog_keys = getParallelJoblogKeys(test_dir)
    else:
        joblog_keys = []

    # Get keys from sacct file
    if has_sacct_file:
        sacct_keys = getSacctOutputKeys(test_dir)
    else:
        sacct_keys = []

    # Get keys from param file
    if has_param_file:
        param_keys = getParamKeys(test_dir)
    else:
        param_keys = []

    # Determine keys and sort by key
    keys = sim_dir_keys + \
            ['initial_timestep', 'final_timestep',
            'initial_cell_count', 'final_cell_count',
            'last_alive_timestep',
            'num_apoptosis', 'num_extrusion'] + \
            joblog_keys + sacct_keys + param_keys

    # Determine per ancestor keys if there is more than one initial ancestor in
    # the simulation 0 (we assume that the initial ancestors of simulation 0 is
    # representative for the whole simulation suite).
    initial_ancestors = sim_datas[0]['initial_ancestors'].split()
    if len(initial_ancestors) > 1:
        for initial_ancestor in initial_ancestors:
            keys += ['initial_cell_count_{}'.format(initial_ancestor)]
            keys += ['final_cell_count_{}'.format(initial_ancestor)]
            keys += ['last_alive_timestep_{}'.format(initial_ancestor)]
            keys += ['num_apoptosis_{}'.format(initial_ancestor)]
            keys += ['num_extrusion_{}'.format(initial_ancestor)]

    if len(argv) >= 3:
        keys = argv[2].split(',')

    sort_by_key = 'sim'

    if len(argv) >= 4:
        sort_by_key = argv[3].strip()

    # Sort sim_datas
    sim_datas.sort(key = lambda d: int(d[sort_by_key]))

    # Print data
    printData(sim_datas, keys)
