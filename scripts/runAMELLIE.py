import os
import numpy as np
import argparse
import subprocess
import json
import time

# If need to rerun jobs from an array, remember that files are labelled with integer one lower than
# the array number. For example, `qsub -t 3-3 jobScript.job` will give `outputFile_2.root`.

def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Run AMELLIE simulation and subsequent analysis code for list of sim info')

    parser.add_argument('--macro_repo', '-mar', type=str, dest='macro_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/2p2gL/MC/macros/', help='Folder to save Region-selected root files in.')
    parser.add_argument('--sim_repo', '-sir', type=str, dest='sim_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/2p2gL/MC/sims/', help='Folder to save intial root files from AMELLIE simulations in.')
    parser.add_argument('--splithist_repo', '-shir', type=str, dest='splithist_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/2p2gL/MC/split_hists/', help='Folder to save uncombined root files with tracking information in.')
    parser.add_argument('--tothist_repo', '-thir', type=str, dest='tothist_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/2p2gL/MC/tot_hists/', help='Folder to save recombined root files with tracking information in.')
    parser.add_argument('--stats_repo', '-str', type=str, dest='stats_repo',
                        default='/mnt/lustre/scratch/epp/jp643/AMELLIE/2p2gL/MC/stats/', help='Folder to save stats txt files in.')
    parser.add_argument('--json_repo', '-jsr', type=str, dest='json_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/2p2gL/results/',
                        help='Folder to json file with final stats in.')

    parser.add_argument('--nevts_total', '-N', type=int, dest='nevts_total',
                        default=1000, help='Number of events to simulate for each setting, total')
    parser.add_argument('--nevts_persim', '-n', type=int, dest='nevts_persim',
                        default=500, help='Max number of events to simulate per macro (simulations will be split up to this amount).')
    parser.add_argument('--max_jobs', '-m', type=int, dest='max_jobs',
                        default=70, help='Max number of tasks in an array running at any one time.')
    parser.add_argument('--abs_nominal', '-ab', type=str, dest='abs_nominal',
                        default='[ 1.0, 1.1, 1.0, ]',
                        help='Original absorption scaling from optics file being used. Ex: ABSLENGTH_SCALING: [ 1.5, 0.55, ],.\n\
                            (fyi for ABSLENGTH_SCALING: 0 = LAB, 1 = PPO, 2 = Te-Diol, 3 = bisMSB).')
    parser.add_argument('--rs_nominal', '-rs', type=str, dest='rs_nominal',
                        default='[1.176,]',
                        help='Original reemission and scattering scaling from optics file being used. Ex: RSLENGTH_SCALING: [1.176,],.\n\
                            (fyi for RSLENGTH_SCALING: 0 = scale full scintilaltor mix I think).')

    parser.add_argument('--list', '-l', type=str, dest='list_file',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/2p2gL/list_info/list.txt', help='Text file with list of sim stats. Format per line:\n\
                            geo_file.geo, LED, fibre, reemis, abs')
    parser.add_argument('--region_lims', '-r', type=str, dest='region_lims',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/2p2gL/list_info/lim_list.txt', help='Text file with list of region limits to apply. Format per line:\n\
                            direct_x_max, reflected_x_min, direct_y_centre, direct_dy, reflected_y_centre, reflected_dy')
    parser.add_argument('---step', '-s', type=str, dest='step',
                    default='all', choices=['sim', 'hist', 'combi', 'stats', 'sim-hist', 'hist-combi',
                                            'combi-stats', 'hist-combi-stats',  'all'],
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

def filename_format(info, Analysis=False):
    '''returns string with simulation info to use in filenames'''
    # geo_file[:-4] + '_' + inner_av_material + '_' + LED + '_' + fibre + 'reemis' + reemis + '_abs' + abs
    name = 'AMELLIE_' + info[0][:-4] + '_' + info[1] + '_' + info[2] + '_' + info[3] + '_reemis' + info[4]
    if Analysis:
        return name
    else:
        return name + '_abs' + info[5]

def job_str_map(jobName_str, info):
    '''Create code string to put at the start of the job file name, so that it can
    be recognised in the job list (the job list only displays the first 10 characters).'''

    map = {
        'job_name': {
            'sims_': 'M',
            'splitHist_': 'H',
            'tot_hists_': 'T',
            'stats_': 'S'
        },
        'geo_file': {
            'snoplusnative.geo': 'a',
            'snoplusnative_te.geo': 'b'
        },
        'inner_av_material': {
            '': 'A',
            'labppo_1p1_berkeley_scintillator': 'B',
            'labppo_2p2_scintillator': 'C'
        },
        'LED': {
            'LED403': 'a'
        },
        'fibre': {
            'FA108': 'A',
            'FA093': 'B'
        }
    }

    # 'geo_file=', info[0], ', inner_av_material' = info[1], ', wavelength=', info[2], ', fibre=', info[3], ', reemis=', info[4], ', abs=', info[5]
    string = map['job_name'][jobName_str] + map['geo_file'][info[0]] + map['inner_av_material'][info[1]] + map['LED'][info[2]] + map['fibre'][info[3]] + info[4][2:]
    # If need one more character, can remove geo file character I guess. Anyway, only 9 characters needed max so far.
    # the '0.' from the reemission fraction is removed since it will always be between 0 and 1
    string += info[5]  # 1.05 = 4 characters

    return string

def makeJobArrayScript(jobName_str, example_jobScript, overall_folder, commandList_address, info, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'job_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += job_str_map(jobName_str, info) + '.job'

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
    new_job_address += job_str_map(jobName_str, info) + '.job'

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

def checkJobsDone(jobName_substr, input_info, wait_time, verbose):
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
                    map_str = job_str_map(jobName_substr, info)
                    if len(map_str) > 10:
                        map_str = map_str[:9]
                    if map_str in line:
                        running = True
                        if verbose:
                            print('Waiting for jobs to finish...')
                        break
        time.sleep(wait_time)

    return True


### Simulation functions ###

def makeMacros(example_macro, save_macro_folder, save_sims_folder, abs_nominal, rs_nominal, info, nevts, idx):
    '''Make and save macro according to provided parameters'''

    # Modify absorption by factor abs_scaling_fact
    abs_scaling_fact = float(info[5])
    abs_nominal = abs_nominal.strip('[').strip(']').strip(' ').split(',')
    abs_scaling = '['
    for absorb in abs_nominal:
        if absorb != '':
            abs_scaling += str(float(absorb) * abs_scaling_fact) + ', '
    abs_scaling += ']'

    # Modify reemission and scattering by factor rs_scaling_fact
    rs_scaling_fact = float(info[4])
    rs_nominal = rs_nominal.strip('[').strip(']').strip(' ').split(',')
    rs_scaling = '['
    for reemission in rs_nominal:
        if reemission != '':
            rs_scaling += str(float(reemission) * rs_scaling_fact) + ', '
    rs_scaling += ']'

    # AMELLIE_geoFile_LEDnum_fibre_reemis_abs.root
    new_macro_address = save_macro_folder + 'macro_' + filename_format(info) + '_' + str(idx) + '.mac'
    output_address = save_sims_folder + 'simOut_' + filename_format(info) + '_' + str(idx) + '.root'

    new_macro = []
    for line in example_macro:
        new_line = line
        # Replace placeholders in macro
        if 'geo/YOUR_FAV_GEO_FILE.geo' in line:
            new_line = new_line.replace('YOUR_FAV_GEO_FILE.geo', info[0], 1)
        if 'YOUR_INNER_AV_MATERIAL' in line:
            new_line = new_line.replace('YOUR_INNER_AV_MATERIAL', info[1], 1)
        if 'YOUR_FAV_LED' in line:
            new_line = new_line.replace('YOUR_FAV_LED', info[2], 1)
        if 'YOUR_FAV_FIBRE' in line:
            new_line = new_line.replace('YOUR_FAV_FIBRE', info[3], 1)
        if '/PATH/TO/YOUR/OUTPUT/FILE.root' in line:
            new_line = new_line.replace('/PATH/TO/YOUR/OUTPUT/FILE.root', output_address, 1)
        if '/rat/run/start 2000' in line:
            new_line = new_line.replace('/rat/run/start 2000', '/rat/run/start ' + str(nevts), 1)
        if 'ABSLENGTH_SCALING [LAB, PPO, Te-Diol, bisMSB,]' in line:
            new_line = new_line.replace('[LAB, PPO, Te-Diol, bisMSB,]', abs_scaling, 1)
        if 'RSLENGTH_SCALING [SCINT-MIX,]' in line:
            new_line = new_line.replace('[SCINT-MIX,]', rs_scaling, 1)
            

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
    #abs_nominal = [0.95, 1.0, 1.0, 2.3]  #(fyi for ABSLENGTH_SCALING: 0 = LAB, 1 = PPO, 2 = Te-Diol, 3 = bisMSB)

    # How to split up sims into manageable macros
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)
    
    ### MAKE MACROS AND JOB SCRIPTS TO RUN THE SIMULATIONS ###
    print('Creating macros and job scripts...')
    job_addresses = []
    #random_seed = '-1829418327'  # option to make sure all sims have the same random seed, to isolate changes
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', inner_av_material=', info[1], ', LED=', info[2], ', fibre=', info[3], ', reemis=', info[4], ', abs=', info[5])
        # Make list of commands for job array to call
        commandList_address = jobScript_repo + 'sim_commandList_' + filename_format(info) + '.txt'
        commandList_file = open(commandList_address, 'w')
        for i in range(len(n_evts)):
            # Create all the macros
            macro_address = makeMacros(example_macro, save_macro_folder, save_sims_folder, args.abs_nominal, args.rs_nominal, info, n_evts[i], i)
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

    # Make sure folders are of the correct format to  use later
    save_sims_folder = checkRepo(args.sim_repo, args.verbose)
    save_splithists_folder = checkRepo(args.splithist_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_sims_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    # How to split up sims into manageable macros
    n_evts = getNevtsPerMacro(args.nevts_total, args.nevts_persim)

    ### Check that all simulations have finished running ###
    checkJobsDone('sims_', input_info, 10, args.verbose)

    ### MAKE JOB SCRIPTS TO CREATE HISTOGRAMS ###
    print('Creating split hist job scripts...')
    hist_command_base = repo_address + 'scripts/GetHists.exe '
    job_addresses = []
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', inner_av_material=', info[1], ', LED=', info[2], ', fibre=', info[3], ', reemis=', info[4], ', abs=', info[5])
        # Make list of commands for job array to call
        commandList_address = jobScript_repo + 'hist_commandList_' + filename_format(info) + '.txt'
        commandList_file = open(commandList_address, 'w')
        wavelength = info[2][3:]
        for i in range(len(n_evts)):
            # Create all the histogram making commands
            hist_command = hist_command_base + save_sims_folder + 'simOut_' + filename_format(info) + '_' + str(i) + '.root '\
                                             + save_splithists_folder + 'splitHist_' + filename_format(info) + '_' + str(i) + '.root '\
                                             + info[3] + ' ' + wavelength + ' ' + str(int(args.verbose))
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

    return True

def combiHists(args, input_info):
    '''Combine split histograms together'''
    print('Running combiHists().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    # Make sure folders are of the correct format to  use later
    save_splithists_folder = checkRepo(args.splithist_repo, args.verbose)
    save_tothists_folder = checkRepo(args.tothist_repo, args.verbose)

    jobSingleScript_address = repo_address + 'job_scripts/jobSingle.job'
    with open(jobSingleScript_address, "r") as f:
        example_jobSingleScript = f.readlines()

    # Making job scripts
    hist_command_base = 'hadd ' + save_tothists_folder + 'tot_hists_'
    job_addresses = []
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', inner_av_material=', info[1], ', LED=', info[2], ', fibre=', info[3], ', reemis=', info[4], ', abs=', info[5])
        # Make list of command for job to call
        command = hist_command_base + filename_format(info) + '.root ' + save_splithists_folder + 'splitHist_' + filename_format(info) + '_*.root'

        # Create the job script to run all these macros in an array
        new_job_address = makeJobSingleScript('tot_hists_', example_jobSingleScript, save_splithists_folder, command, info, args.verbose)
        job_addresses.append(new_job_address)

    # Wait until these job arrays are done
    checkJobsDone('splitHist_', input_info, 10, args.verbose)
    
    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for job_address in job_addresses:
        command = 'qsub ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

### Analysis functions ###

def getOuterLims(region_lims_list, verbose):
    '''Create string of outer region limits with a space between each, to be used in commands.'''

    outer_lims_list = []
    for region_lims in region_lims_list:
        if verbose:
            print('direct_x_max=', region_lims[0], ', reflected_x_min=', region_lims[1], ', direct_y_centre=', region_lims[2], ', direct_dy=', region_lims[3], ', reflected_y_centre=', region_lims[4], 'reflected_dy=', region_lims[5])
        # Package region limits in string format to use in commands
        # region_lims = [direct_x_max, reflected_x_min, direct_y_centre, direct_dy, reflected_y_centre, reflected_dy]
        direct_x_max = str(region_lims[0])
        direct_x_min = str(-1.0)
        direct_y_max = str(region_lims[2] + 0.5 * region_lims[3])
        direct_y_min = str(region_lims[2] - 0.5 * region_lims[3])
        reflected_x_max = str(1.0)
        reflected_x_min = str(region_lims[1])
        reflected_y_max = str(region_lims[4] + 0.5 * region_lims[5])
        reflected_y_min = str(region_lims[0] - 0.5 * region_lims[5])
        outer_lims_list.append(direct_x_max + ' ' + direct_x_min + ' ' + direct_y_max + ' ' + direct_y_min + ' ' + reflected_x_max + ' ' + reflected_x_min + ' ' + reflected_y_max + ' ' + reflected_y_min)
    
    return outer_lims_list

def makeStatsJobScript(hist_files, info, outer_lims_list, repo_address, save_stats_folder, save_tothists_folder, verbose):
    '''Create job script and associated command list file to run Analysis on array of same simulations
    (except absorption factor), with different region limits.'''

    # Folder for job scripts and command lists they use
    jobScript_repo = save_tothists_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, verbose)

    # For each set of absorptions, apply all the region limits
    commandList_address = jobScript_repo + 'stats_commandList_' + filename_format(info, True) + '.txt'
    commandList_file = open(commandList_address, 'w')
    output_stats_files = []
    for i in range(len(outer_lims_list)):
        outer_lims = outer_lims_list[i]
        if verbose:
            print('outer_lims=', outer_lims)

        output_stats_file = save_stats_folder + 'stats_' + filename_format(info, True) + '_' + outer_lims.replace(' ', '_') + '.txt'
        slope_command = repo_address + 'scripts/GetStats.exe ' + output_stats_file + ' ' + outer_lims + ' ' + str(int(verbose)) + ' ' + hist_files
        commandList_file.write(slope_command)
        output_stats_files.append(output_stats_file)
        if i != len(outer_lims_list) - 1:
            commandList_file.write('\n')
    commandList_file.close()

    jobArrayScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(jobArrayScript_address, "r") as f:
        example_jobArrayScript = f.readlines()

    job_address = makeJobArrayScript('stats_', example_jobArrayScript, save_tothists_folder, commandList_address, info, verbose)

    return job_address, output_stats_files

def writeJsonStatsFile(output_info, region_lims_list, json_stats_folder):
    '''Write all stats to json files'''

    for output_stats_file_list, info, abs_factors in output_info:
        # Write sim info
        table = {}
        table['geo_file'] = info[0]
        table['inner_av_material'] = info[1]
        table['LED'] = info[2]
        table['fibre'] = info[3]
        table['reemis'] = info[4]

        # Write results for each region applied
        table['regions'] = []
        for i in range(len(output_stats_file_list)):
            output_stats_file = output_stats_file_list[i]
            region_lims = region_lims_list[i]
            table['regions'].append({})

            # Read in results to make json table for easier use
            with open(output_stats_file, "r") as f:
                stats = f.readlines()

            table['regions'][i]['region_lims'] = {}
            table['regions'][i]['region_lims']['direct_x_max'] = region_lims[0]
            table['regions'][i]['region_lims']['reflected_x_min'] = region_lims[1]
            table['regions'][i]['region_lims']['direct_y_centre'] = region_lims[2]
            table['regions'][i]['region_lims']['direct_dy'] = region_lims[3]
            table['regions'][i]['region_lims']['reflected_y_centre'] = region_lims[4]
            table['regions'][i]['region_lims']['reflected_dy'] = region_lims[5]

            # Get stats for each individual absorption factor
            table['regions'][i]['stats'] = {}
            for j in range(len(stats)):
                absorb = abs_factors[j]
                stats_info = stats[j].split(' ')  # [tot_hits, direct_hits, reflected_hits]
                table['regions'][i]['stats'][absorb] = {}
                table['regions'][i]['stats'][absorb]['tot_hits'] = int(stats_info[0])
                table['regions'][i]['stats'][absorb]['direct_hits'] = int(stats_info[1])
                table['regions'][i]['stats'][absorb]['reflected_hits'] = int(stats_info[2])

        # Save table to json file
        save_file = json_stats_folder + 'FinalStats_' + filename_format(info, True) + '.json'
        with open(save_file, 'w') as f:
            json.dump(table, f)

def getStats(args, input_info):
    '''Compute stats (number of hits in direct and reflected regions, and total hits) from different absorption scalings.
    If verbose flag is true, also create 2-D hists with regions overlain.'''
    print('Running getStats().')

    # Read in example macro and job script + info
    repo_address = getRepoAddress()

    # Make sure folders are of the correct format to  use later
    save_tothists_folder = checkRepo(args.tothist_repo, args.verbose)
    save_stats_folder = checkRepo(args.stats_repo, args.verbose)
    json_stats_folder = checkRepo(args.json_repo, args.verbose)

    # Wait until previous jobs are done
    checkJobsDone('tot_hists_', input_info, 10, args.verbose)

    # Apply region limits from list
    if args.verbose:
        print('Reading in region limits...')
    textFile = open(args.region_lims, "r")
    lines = [line.split(', ') for line in textFile]
    textFile.close()
    region_lims_list = np.asarray(lines).astype(np.float)

    # Get outer region limit strings
    outer_lims_list = getOuterLims(region_lims_list, args.verbose)


    ### MAKE JOB SCRIPTS TO RUN ANALYSIS ###
    print('Creating analysis job scripts...')
    jobScript_addresses = []
    output_info = []
    abs_factors = []
    info_1 = input_info[0]
    hist_files = ''
    for i in range(len(input_info)):
        info = input_info[i]
        if args.verbose:
            print('geo_file=', info[0], ', inner_av_material=', info[1], ', LED=', info[2], ', fibre=', info[3], ', reemis=', info[4], ', abs=', info[5])
        # Check whether to add arguments to current command, or start new one
        is_same = info_1[0] == info[0] and info_1[1] == info[1] and info_1[2] == info[2] and info_1[3] and info_1[4] == info[4]
        tot_hist_file_address = save_tothists_folder + 'tot_hists_' + filename_format(info) + '.root'

        # Group up files which are the same except for absorption scaling
        if is_same:
            hist_files += ' ' + tot_hist_file_address
            abs_factors.append(info[5])
        else:
            jobScript_address, output_stats_files = makeStatsJobScript(hist_files, info, outer_lims_list, repo_address, save_stats_folder, save_tothists_folder, args.verbose)
            jobScript_addresses.append(jobScript_address)
            output_info.append((output_stats_files, info_1, abs_factors))
            abs_factors = []
            hist_files = tot_hist_file_address
            info_1 = info
            abs_factors.append(info[5])
        if i == len(input_info) - 1:
            jobScript_address, output_stats_files = makeStatsJobScript(hist_files, info, outer_lims_list, repo_address, save_stats_folder, save_tothists_folder, args.verbose)
            jobScript_addresses.append(jobScript_address)
            output_info.append((output_stats_files, info_1, abs_factors))


    ### RUN JOB SCRIPTS ###
    print('Submitting job(s)...')
    for job_address in jobScript_addresses:
        command = 'qsub -t 1-' + str(len(outer_lims_list)) + ' -tc ' + str(args.max_jobs) + ' ' + job_address
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    # Wait until these job arrays are done
    checkJobsDone('stats_', input_info, 10, args.verbose)

    ### READ STATS AND WRITE TO JSON FILE ###
    print('Writing stats to json file...')
    writeJsonStatsFile(output_info, region_lims_list, json_stats_folder)

    return True


### Combined functions ###

def runSims_getHists(args, input_info):
    '''Run AMELLIE simulations, then get their histograms'''
    result = runSims(args, input_info)
    return result and getHists(args, input_info)

def getHists_combine(args, input_info):
    '''Get simulations histograms and recombine them'''
    result = getHists(args, input_info)
    return result and combiHists(args, input_info)

def combine_getStats(args, input_info):
    '''Recombine histograms and analyse them'''
    result = combiHists(args, input_info)
    return result and getStats(args, input_info)

def getHists_Combine_getStats(args, input_info):
    '''Get histograms from simulation outputs, combine and then analyse them.'''
    result = getHists(args, input_info)
    result = result and combiHists(args, input_info)
    return result and getStats(args, input_info)

def runAll(args, input_info):
    '''Run the whole AMELLIE code, from simulation, through histograms, combining, to analysis'''
    result = runSims(args, input_info)
    result = result and getHists(args, input_info)
    result = result and combiHists(args, input_info)
    return result and getStats(args, input_info)


### MAIN ###

def main():
    # read in argument
    args = argparser()

    # read in file info
    textFile = open(args.list_file, "r")
    lines = [line.replace('\n', '').split(', ') for line in textFile]
    textFile.close()
    input_info = np.asarray(lines)

    # geo_files = input_info[:, 0]
    # inner_av_material = input_info[:, 1]
    # wavelengths = input_info[:, 2]
    # fibres = input_info[:, 3]
    # reemissions = input_info[:, 4]
    # abs_factors = input_info[:, 5]

    work_modes = {
        'sim': runSims,
        'hist': getHists,
        'combi': combiHists,
        'stats': getStats,

        'sim-hist': runSims_getHists,
        'hist-combi': getHists_combine,
        'combi-stats': combine_getStats,
        'hist-combi-stats': getHists_Combine_getStats,
        'all': runAll
    }

    result = work_modes[args.step](args, input_info)


if __name__ == '__main__':
    main()