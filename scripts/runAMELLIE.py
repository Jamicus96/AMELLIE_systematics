import os
import numpy as np
import argparse
import subprocess
import json
import time

# repo address is not getting the full path, only the path from where I run the script

def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run AMELLIE simulation and subsequent analysis code for list of sim info')

    parser.add_argument('--macro_repo', '-mar', type=str, dest='macro_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/macros/', help='Folder to save Region-selected root files in.')
    parser.add_argument('--sim_repo', '-sir', type=str, dest='sim_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/sims/', help='Folder to save intial root files from AMELLIE simulations in.')
    parser.add_argument('--splithist_repo', '-shir', type=str, dest='splithist_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/split_hists/', help='Folder to save uncombined root files with tracking information in.')
    parser.add_argument('--tothist_repo', '-thir', type=str, dest='tothist_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/tot_hists/', help='Folder to save recombined root files with tracking information in.')
    parser.add_argument('--stats_repo', '-str', type=str, dest='stats_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/stats/', help='Folder to save stats txt files in.')
    parser.add_argument('--json_repo', '-jsr', type=str, dest='json_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/Sim/Slopes/final_stats/',
                        help='Folder to json file with final stats in.')

    parser.add_argument('--nevts_total', '-N', type=int, dest='nevts_total',
                        default=1000, help='Number of events to simulate for each setting, total')
    parser.add_argument('--nevts_persim', '-n', type=int, dest='nevts_persim',
                        default=500, help='Max number of events to simulate per macro (simulations will be split up to this amount).')
    parser.add_argument('--max_jobs', '-m', type=int, dest='max_jobs',
                        default=70, help='Max number of tasks in an array running at any one time.')
    parser.add_argument('--region_lims', '-r', type=list, dest='region_lims',
                        default=[-0.95, 0.85, -0.95, 20,  -4, 20],
                        help='List of region limits: [direct_x_max, reflected_x_min, direct_y_centre, direct_dy, reflected_y_centre, reflected_dy]')

    parser.add_argument('--list', '-l', type=str, dest='list_file',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/Sim/Slopes/info_lists/list.txt', help='Text file with list of sim stats. Format per line:\n\
                            geo_file.geo, LED, fibre, reemis, abs')
    parser.add_argument('---step', '-s', type=str, dest='step',
                    default='allSlope', choices=['sim', 'hist', 'sim-hist', 'FOM', 'slope', 'hist-FOM',
                                            'hist-Slopes', 'allFOM', 'allSlope'],
                    help='which step of the process is it in?')
    parser.add_argument('---verbose', '-v', type=bool, dest='verbose',
                    default=True, help='print and save extra info')

    args = parser.parse_args()
    return args


### Miscelaneous functions ###

def getRepoAddress():
    '''Returns the full address of the git repo containing with script'''
    repo_address = __file__[:-len('scripts/runAMMELIE.py')]
    if repo_address == '':
        firt_char = None
    else:
        firt_char = repo_address[0]
    if firt_char != '/':
        repo_address = os.getcwd() + '/' + repo_address
    return repo_address

def checkRepo(repo_address, verbose=False):
    '''Check format of repo address is right to be used here. Also if repo does not exist, create it.'''

    # Format address
    new_address = repo_address
    if new_address == '.':
        new_address = ''
    elif new_address != '':
        if (new_address[int(len(new_address) - 1)] != '/'):
            new_address += '/'
        
        # If directory does not exists, make it
        if not os.path.isdir(new_address):
            os.mkdir(new_address)
            if verbose:
                print('Created new directory: ', new_address)

    return new_address

def filename_format(info, withoutAbs=False):
    '''returns string with simulation info to use in filenames'''
    # geo_file[:-4] + '_' + LED + '_' + fibre + 'reemis' + reemis + '_abs' + abs
    name = 'AMELLIE_' + info[0][:-4] + '_' + info[1] + '_' + info[2] + '_reemis' + info[3]
    if withoutAbs:
        return name
    else:
        return name + '_abs' + info[4]

def job_str_map(jobName_str, info, isArray):
    '''Create code string to put at the start of the job file name, so that it can
    be recognised in the job list (the job list only displays the first 10 characters).'''

    map = {
        'job_name': {
            'sims_': 'M',
            'splitHist_': 'H',
            'tot_hists_': 'T',
            'slopes_': 'S'
        },
        'geo_file': {
            'snoplusnative_te.geo': 'a'
        },
        'LED': {
            'LED403': 'A'
        },
        'fibre': {
            'FA108': 'a'
        }
    }

    # 'geo_file=', info[0], ', wavelength=', info[1], ', fibre=', info[2], ', reemis=', info[3], ', abs=', info[4]
    string = map['job_name'][jobName_str] + map['geo_file'][info[0]] + map['LED'][info[1]] + map['fibre'][info[2]] + info[3][2:]
    # If need one more character, can remove geo file character I guess. Anyway, only 9 characters needed max so far.
    # the '0.' from the reemission fraction is removed since it will always be between 0 and 1
    if isArray:
        string += info[4]  # 1.05 = 4 characters

    return string

def makeJobArrayScript(jobName_str, example_jobScript, overall_folder, commandList_address, info, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += job_str_map(jobName_str, info, True) + '.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_' + jobName_str + filename_format(info) + '.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif '/Address/CommandList.txt' in line:
            new_line = line.replace('/Address/CommandList.txt', commandList_address, 1)
        else:
            new_line = line

        new_jobScript.append(new_line)

    # Create job file
    with open(new_job_address, "w") as f:
        new_jobScript = "".join(new_jobScript)
        f.write(new_jobScript)

    return new_job_address

def makeJobSingleScript(jobName_str, example_jobScript, overall_folder, commands, info, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += job_str_map(jobName_str, info, False) + '.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  'log_' + jobName_str + filename_format(info) + '.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif 'your commands' in line:
            new_line = line.replace('your commands', commands, 1)
        else:
            new_line = line

        new_jobScript.append(new_line)

    # Create job file
    with open(new_job_address, "w") as f:
        new_jobScript = "".join(new_jobScript)
        f.write(new_jobScript)

    return new_job_address

def getNevtsPerMacro(nevts_total, nevts_persim):
    '''Create array with number of events to simulate per macro'''

    n_macros = nevts_total // nevts_persim
    remainder = nevts_total % nevts_persim
    if remainder == 0:
        n_evts = np.an_array = np.full(n_macros, nevts_persim)
    else:
        n_evts = np.an_array = np.full(n_macros + 1, nevts_persim)
        n_evts[n_macros] = remainder

    return n_evts

def checkJobsDone(jobName_substr, input_info, wait_time, isArray):
    '''Wait until submitted jobs of certain forma are finished. Wait time in seconds.'''

    # Turns out the name of the job is only the 10 first characters of the job file name

    running = True
    while running:
        running = False
        output = subprocess.Popen('qstat -u $USER', stdout=subprocess.PIPE, shell=True).communicate()[0]
        lines = output.decode("utf-8").split('\n')
        for line in lines:
            if running:
                break
            else:
                for info in input_info:
                    if job_str_map(jobName_substr, info, isArray) in line:
                        running = True
                        break
        time.sleep(wait_time)

    return True


### Simulation functions ###

def makeMacros(example_macro, save_macro_folder, save_sims_folder, abs_nominal, info, nevts, idx):
    '''Make and save macro according to provided parameters'''

    # Modify by factor abs_scaling_fact
    abs_scaling_fact = float(info[4])
    abs_scaling = '[' + str(abs_nominal[0] * abs_scaling_fact)
    abs_scaling += ', ' + str(abs_nominal[1] * abs_scaling_fact)
    abs_scaling += ', ' + str(abs_nominal[2] * abs_scaling_fact)
    abs_scaling += ', ' + str(abs_nominal[3] * abs_scaling_fact) + ']'

    # AMELLIE_geoFile_LEDnum_fibre_reemis_abs.root
    new_macro_address = save_macro_folder + 'macro_' + filename_format(info) + '_' + str(idx) + '.mac'
    output_address = save_sims_folder + 'simOut_' + filename_format(info) + '_' + str(idx) + '.root'

    new_macro = []
    for line in example_macro:
        # Replace placeholders in macro
        if 'geo/YOUR_FAV_GEO_FILE.geo' in line:
            new_line = line.replace('YOUR_FAV_GEO_FILE.geo', info[0], 1)
        elif 'YOUR_FAV_LED' in line:
            new_line = line.replace('YOUR_FAV_LED', info[1], 1)
        elif 'YOUR_FAV_FIBRE' in line:
            new_line = line.replace('YOUR_FAV_FIBRE', info[2], 1)
        elif '/PATH/TO/YOUR/OUTPUT/FILE.root' in line:
            new_line = line.replace('/PATH/TO/YOUR/OUTPUT/FILE.root', output_address, 1)
        elif '/rat/run/start 2000' in line:
            new_line = line.replace('/rat/run/start 2000', '/rat/run/start ' + str(nevts), 1)
        elif 'ABSLENGTH_SCALING [1.9, 2.0, 2.0, 0.0]' in line:
            new_line = line.replace('[1.9, 2.0, 2.0, 0.0]', abs_scaling, 1)
        else:
            new_line = line

        new_macro.append(new_line)
    
    # Create new macro file
    with open(new_macro_address, "w") as f:
        new_macro = "".join(new_macro)
        f.write(new_macro)

    return new_macro_address

def runSims(args, input_info):
    '''Runs AMELLIE simulations based in input information'''
    print('Running runSims().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    example_macro_address = repo_address + 'macros/runSimulation.mac'
    with open(example_macro_address, "r") as f:
        example_macro = f.readlines()

    example_jobScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(example_jobScript_address, "r") as f:
        example_jobScript = f.readlines()

    # Make sure folders are of the correct format to  use later
    save_macro_folder = checkRepo(args.macro_repo, args.verbose)
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_macro_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    # Nominal ABSLENGTH_SCALING (see rat/data/OPTICS_TeDiol_0p5_bismsb_dda.ratdb)
    abs_nominal = [0.95, 1.0, 1.0, 2.3]  #(fyi for ABSLENGTH_SCALING: 0 = LAB, 1 = PPO, 2 = Te-Diol, 3 = bisMSB)

    # How to split up sims into manageable macros
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)
    
    ### MAKE MACROS AND JOB SCRIPTS TO RUN THE SIMULATIONS ###
    print('Creating macros and job scripts...')
    job_addresses = []
    #random_seed = '-1829418327'  # option to make sure all sims have the same random seed, to isolate changes
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', LED=', info[1], ', fibre=', info[2], ', reemis=', info[3], ', abs=', info[4])
        # Make list of commands for job array to call
        commandList_address = jobScript_repo + 'sim_commandList_' + filename_format(info) + '.txt'
        commandList_file = open(commandList_address, 'w')
        for i in range(len(n_evts)):
            # Create all the macros
            macro_address = makeMacros(example_macro, save_macro_folder, save_sims_folder, abs_nominal, info, n_evts[i], i)
            log_file_address = save_macro_folder + 'log_files/ratLog_' + filename_format(info) + '.log'
            macro_command = 'rat ' + macro_address + ' -l ' + log_file_address # + ' -s ' + random_seed
            commandList_file.write(macro_command + '\n')
        commandList_file.close()

        # Create the job script to run all these macros in an array
        new_job_address = makeJobArrayScript('sims_', example_jobScript, save_macro_folder, commandList_address, info, args.verbose)
        job_addresses.append(new_job_address)


    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for job_address in job_addresses:
        command = 'qsub -t 1-' + str(len(n_evts)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address 
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True


### Histogram functions ###

def getHists(args, input_info):
    '''Get histograms from split up simulation outputs, and recombine them'''
    print('Running getHists().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    jobArrayScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(jobArrayScript_address, "r") as f:
        example_jobArrayScript = f.readlines()

    jobSingleScript_address = repo_address + 'job_scripts/jobSingle.job'
    with open(jobSingleScript_address, "r") as f:
        example_jobSingleScript = f.readlines()

    # Make sure folders are of the correct format to  use later
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)
    save_splithists_folder = checkRepo(args.splithist_repo, args.verbose)
    save_tothists_folder = checkRepo(args.tothist_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_sims_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    # How to split up sims into manageable macros
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)

    ### Check that all simulations have finished running ###
    checkJobsDone('sims_', input_info, 10, True)

    ### MAKE JOB SCRIPTS TO CREATE HISTOGRAMS ###
    print('Creating split hist job scripts...')
    hist_command_base = repo_address + 'scripts/GetHists.exe '
    job_addresses = []
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', LED=', info[1], ', fibre=', info[2], ', reemis=', info[3], ', abs=', info[4])
        # Make list of commands for job array to call
        commandList_address = jobScript_repo + 'hist_commandList_' + filename_format(info) + '.txt'
        commandList_file = open(commandList_address, 'w')
        wavelength = info[1][3:]
        for i in range(len(n_evts)):
            # Create all the histogram making commands
            hist_command = hist_command_base + save_sims_folder + 'simOut_' + filename_format(info) + '_' + str(i) + '.root '\
                                             + save_splithists_folder + 'splitHist_' + filename_format(info) + '_' + str(i) + '.root '\
                                             + info[2] + ' ' + wavelength + ' ' + str(int(args.verbose))
            commandList_file.write(hist_command + '\n')
        commandList_file.close()

        # Create the job script to run all these macros in an array
        new_job_address = makeJobArrayScript('splitHist_', example_jobArrayScript, save_sims_folder, commandList_address, info, args.verbose)
        job_addresses.append(new_job_address)

    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for job_address in job_addresses:
        command = 'qsub -t 1-' + str(len(n_evts)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    # Wait until these job arrays are done
    checkJobsDone('splitHist_', input_info, 10, True)

    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

    ### MAKE JOB SCRIPTS TO COMBINE HISTOGRAMS ###
    print('Creating combine hists job scripts...')
    hist_command_base = 'hadd ' + save_tothists_folder + 'tot_hists_'
    job_addresses = []
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', LED=', info[1], ', fibre=', info[2], ', reemis=', info[3], ', abs=', info[4])
        # Make list of command for job to call
        command = hist_command_base + filename_format(info) + '.root ' + save_splithists_folder + 'splitHist_' + filename_format(info) + '_*.root'

        # Create the job script to run all these macros in an array
        new_job_address = makeJobSingleScript('tot_hists_', example_jobSingleScript, save_splithists_folder, command, info, args.verbose)
        job_addresses.append(new_job_address)

    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for job_address in job_addresses:
        command = 'qsub ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True


### Analysis functions ###

def MakeSlopeCommand(input_info, line1, line2, region_lims, repo_address, example_jobScript, save_stats_folder, save_tothists_folder, args):
    '''Create string of command to run slope code'''

    # Create command
    slope_command_base = repo_address + 'scripts/GetFOMabsSlope.exe '
    info_str = filename_format(input_info[line1, :], True)  # Same as usual but without absorption
    output_stats_file = save_stats_folder + 'slopeStats_' + info_str + '.txt'
    output_root_file = save_stats_folder + 'Regions_' + info_str + '.root'
    slope_command = slope_command_base + args.list_file + ' ' + str(line1) + ' ' + str(line2) + ' ' + save_tothists_folder + ' '\
                    + output_stats_file + ' ' + str(int(args.verbose)) + ' ' + output_root_file + ' ' + region_lims

    # Create the job script to run all these macros in an array
    new_job_address = makeJobSingleScript('slopes_', example_jobScript, save_tothists_folder, slope_command, input_info[0], args.verbose)
    return new_job_address, output_stats_file

def getSlopes(args, input_info):
    '''Compute slope from FOM computed from different absorption scalings. If verbose flag is
    true, also records FOM for each different abs file. So, if only an FOM is wanted, just run this
    and ignore the slope result.'''
    print('Running getSlopes().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    jobScript_address = repo_address + 'job_scripts/jobSingle.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()

    # Make sure folders are of the correct format to  use later
    save_tothists_folder = checkRepo(args.tothist_repo, args.verbose)
    save_stats_folder = checkRepo(args.stats_repo, args.verbose)
    json_stats_folder = checkRepo(args.json_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_tothists_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    # Package region limits in string format to use in commands
    # args.region_lims = [direct_x_max, reflected_x_min, direct_y_centre, direct_dy, reflected_y_centre, reflected_dy]
    direct_x_max = args.region_lims[0]
    direct_x_min = -1.0
    direct_y_max = args.region_lims[2] + 0.5 * args.region_lims[3]
    direct_y_min = args.region_lims[2] - 0.5 * args.region_lims[3]
    reflected_x_max = 1.0
    reflected_x_min = args.region_lims[1]
    reflected_y_max = args.region_lims[4] + 0.5 * args.region_lims[5]
    reflected_y_min = args.region_lims[0] - 0.5 * args.region_lims[5]
    region_lims = str(direct_x_max) + ' ' + str(direct_x_min) + ' ' + str(direct_y_max) + ' ' + str(direct_y_min) + ' '\
                + str(reflected_x_max) + ' ' + str(reflected_x_min) + ' ' + str(reflected_y_max) + ' ' + str(reflected_y_min)

    # Wait until previous jobs are done
    checkJobsDone('tot_hists_', input_info, 10, False)

    ### MAKE JOB SCRIPTS TO RUN ANALYSIS ###
    print('Creating analysis job scripts...')
    
    line1 = 0
    line2 = 0
    i = 0
    lines = []
    job_addresses = []
    output_stats_files = []
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', LED=', info[1], ', fibre=', info[2], ', reemis=', info[3], ', abs=', info[4])

        # Check if all info except absorption is still the same
        is_same = input_info[line1, 0] == info[0] and input_info[line1, 1] == info[1] and input_info[line1, 2] == info[2] and input_info[line1, 3] == info[3]
        if not is_same:
            new_job_address, new_output_stats_file = MakeSlopeCommand(input_info, line1, line2, region_lims, repo_address, example_jobScript, save_stats_folder, save_tothists_folder, args)
            job_addresses.append(new_job_address)
            output_stats_files.append(new_output_stats_file)
            lines.append(line1)
            line1 = line2 + 1
        line2 = i
        i += 1
    new_job_address, new_output_stats_file = MakeSlopeCommand(input_info, line1, line2, region_lims, repo_address, example_jobScript, save_stats_folder, save_tothists_folder, args)
    job_addresses.append(new_job_address)
    output_stats_files.append(new_output_stats_file)
    lines.append(line1)

    ### RUN JOB SCRIPTS ###
    print('Submitting job(s)...')
    for job_address in job_addresses:
        command = 'qsub ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    # Wait until these job arrays are done
    checkJobsDone('slopes_', input_info, 10, False)

    i = 0
    for output_stats_file in output_stats_files:
        # Read in results to make json table for easier use
        with open(output_stats_file, "r") as f:
            stats = f.readlines()
        table = {}

        # Get region limits
        outer_lims = stats[0].split(' ')
        # From: [direct_x_max, direct_x_min, direct_y_max, direct_y_min, reflected_x_max, reflected_x_min, reflected_y_max, reflected_y_min]
        # To: [direct_x_max, reflected_x_min, direct_y_centre, direct_dy, reflected_y_centre, reflected_dy]
        table['region_lims'] = {}
        table['region_lims']['direct_x_max'] = float(outer_lims[0])
        table['region_lims']['reflected_x_min'] = float(outer_lims[5])
        table['region_lims']['direct_y_centre'] = 0.5 * (float(outer_lims[2]) + float(outer_lims[3]))
        table['region_lims']['direct_dy'] = float(outer_lims[2]) - float(outer_lims[3])
        table['region_lims']['reflected_y_centre'] = 0.5 * (float(outer_lims[4]) + float(outer_lims[5]))
        table['region_lims']['reflected_dy'] = float(outer_lims[4]) - float(outer_lims[5])

        # Get slope, and FOM of inividual absorptions if want extra info
        table['slope'] = float(stats[1])
        if args.verbose:
            table['FOMs'] = {}
            for j in range(2, len(stats)):
                info = stats[j].split(' ')  # [abs, FOM]
                table['FOMs'][info[0]] = float(info[1])

        # Save table to json file
        save_file = json_stats_folder + 'FinalStats_' + filename_format(input_info[lines[i], :], True) + '.json'
        with open(save_file, 'w') as f:
            json.dump(table, f)
        i += 1

    return True


### Combined functions ###

def runSims_getHists(args, input_info):
    '''Run AMELLIE simulations, then get their histograms and recombine'''
    result = runSims(args, input_info)
    return result and getHists(args, input_info)

def getHists_Slopes(args, input_info):
    '''Get histograms from simulation outputs, combine and then analyse them.'''
    result = getHists(args, input_info)
    return result and getSlopes(args, input_info)

def runAllSlopes(args, input_info):
    '''Run the whole AMELLIE code, from simulation, through histograms, to analysis'''
    result = runSims(args, input_info)
    result = result and getHists(args, input_info)
    return result and getSlopes(args, input_info)


### MAIN ###

def main():
    # read in argument
    args = argparser()

    # read in file info
    textFile = open(args.list_file, "r")
    lines = [line.split(', ') for line in textFile]
    textFile.close()
    input_info = np.asarray(lines)

    # geo_files = input_info[:, 0]
    # wavelengths = input_info[:, 1]
    # fibres = input_info[:, 2]
    # reemissions = input_info[:, 3]
    # abs_factors = input_info[:, 4]

    work_modes = {
        'sim': runSims,
        'hist': getHists,
        'sim-hist': runSims_getHists,
        'slope': getSlopes,
        'hist-Slopes': getHists_Slopes,
        'allSlope': runAllSlopes
    }

    result = work_modes[args.step](args, input_info)


if __name__ == '__main__':
    main()