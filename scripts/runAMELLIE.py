import os
import numpy as np
import argparse
import subprocess
import json
import time


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

    parser.add_argument('--nevts_total', '-N', type=int, dest='nevts_total',
                        default=1000, help='Number of events to simulate for each setting, total')
    parser.add_argument('--nevts_persim', '-n', type=int, dest='nevts_persim',
                        default=500, help='Max number of events to simulate per macro (simulations will be split up to this amount).')
    parser.add_argument('--region_lims', '-r', type=list, dest='region_lims',
                        default=[-0.95, 0.85, -0.95, 20,  -4, 20],
                        help='List of region limits: [direct_x_max, reflected_x_min, direct_y_centre, direct_dy, reflected_y_centre, reflected_dy]')

    parser.add_argument('--list', '-l', type=str, dest='list_file',
                        default='list.txt', help='Text file with list of sim stats. Format per line:\n\
                            geo_file.geo, wavelength, fibre, reemis, abs')
    parser.add_argument('---step', '-s', type=str, dest='step',
                    default='all', choices=['sim', 'hist', 'sim-hist', 'FOM', 'slope', 'hist-FOM',
                                            'hist-Slopes', 'allFOM', 'allSlope'],
                    help='which step of the process is it in?')
    parser.add_argument('---verbose', '-v', type=bool, dest='verbose',
                    default=False, help='print and save extra info')

    args = parser.parse_args()
    return args


### Miscelaneous functions ###

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

def filename_format(info):
    '''returns string with simulation info to use in filenames'''
    # geo_file[:-4] + '_' + wavelength + '_' + fibre + 'reemis' + reemis + '_abs' + abs
    return 'AMELLIE_' + info[0][:-4] + '_' + info[1] + '_' + info[2] + 'reemis' + info[3] + '_abs' + info[4]

def makeJobArrayScript(jobName_str, example_jobScript, overall_folder, commandList_address, info, verbose):
    '''Create job script to run array of rat macros'''

    new_job_address = overall_folder + 'jobs_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += jobName_str + '.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  jobName_str + '.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif 'Command here' in line:
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

    new_job_address = overall_folder + 'jobs_scripts/'
    new_job_address = checkRepo(new_job_address, verbose)
    new_job_address += jobName_str + '.job'

    output_logFile_address = overall_folder + 'log_files/'
    output_logFile_address = checkRepo(output_logFile_address, verbose)
    output_logFile_address +=  jobName_str + '.txt'

    new_jobScript = []
    for line in example_jobScript:
        # Replace placeholders in macro
        if 'output_log.txt' in line:
            new_line = line.replace('output_log.txt', output_logFile_address, 1)
        elif 'Command here' in line:
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
    n_evts = np.an_array = np.full(n_macros + 1, nevts_persim)
    n_evts[n_macros] = remainder

    return n_evts

def checkJobsDone(jobName_substr, input_info, wait_time):
    '''Wait until submitted jobs of certain forma are finished. Wait time in seconds.'''

    running = True
    while running:
        running = False
        output = subprocess.Popen('qstat -u $USER', stdout=subprocess.PIPE).communicate()[0]
        lines = output.decode("utf-8").split('\n')
        for line in lines:
            if running:
                break
            else:
                for info in input_info:
                    if jobName_substr + filename_format(info) in line:
                        running = True
                        break
        time.wait(wait_time)

    return True


### Simulation functions ###

def makeMacros(example_macro, save_macro_folder, save_sims_folder, abs_nominal, info, nevts, idx):
    '''Make and save macro according to provided parameters'''

    # Modify by factor abs_scaling_fact
    abs_scaling_fact = float(abs)
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
            new_line = line.replace('/rat/run/start 2000', '/rat/run/start ' + nevts, 1)
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

    # Read in example macro and job script + info
    repo_address = __file__[:-len('scripts/runAMMELIE.py')]

    macro_address = repo_address + 'macros/runSimulation.mac'
    with open(macro_address, "r") as f:
        example_macro = f.readlines()

    jobScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(jobScript_address, "r") as f:
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
            print('geo_file=', info[0], ', wavelength=', info[1], ', fibre=', info[2], ', reemis=', info[3], ', abs=', info[4])
        # Make list of commands for job array to call
        commandList_address = jobScript_repo + 'sim_commandList_' + filename_format(info) + '_' + str(i) + '.root'
        commandList_file = open(commandList_address, 'w')
        for i in range(n_evts):
            # Create all the macros
            new_macro_address = makeMacros(example_macro, save_macro_folder, save_sims_folder, abs_nominal, info, n_evts[i], i)
            #macro_addresses.append(new_macro_address)
            macro_command = 'rat ' + new_macro_address #+ ' -s ' + random_seed
            commandList_file.write(macro_command + '\n')
        commandList_file.close()

        # Create the job script to run all these macros in an array
        new_job_address = makeJobArrayScript('sims_' + filename_format(info), example_jobScript, save_macro_folder, commandList_address, info, args.verbose)
        job_addresses.append(new_job_address)
        
    if len(n_evts) != len(job_addresses):
        print('ERROR: Number of job scripts does not match number needed to run macros')
        print('len(n_evts) = ', len(n_evts), ', len(job_addresses) = ', len(job_addresses))
        exit()


    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for i in range(len(n_evts)):
        command = 'qsub -t 1-' + str(n_evts[i]) + ' ' + job_addresses[i] 
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True


### Histogram functions ###

def getHists(args, input_info):
    '''Get histograms from split up simulation outputs, and recombine them'''

    # Read in example macro and job script + info
    repo_address = __file__[:-len('scripts/runAMMELIE.py')]

    jobScript_address = repo_address + 'job_scripts/jobArray.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()

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
    checkJobsDone('sims_jobScript_', filename_format(info), 10)

    ### MAKE JOB SCRIPTS TO CREATE HISTOGRAMS ###
    print('Creating macros and job scripts...')
    hist_command_base = repo_address + 'scripts/./GetHists.exe '
    job_addresses = []
    for info in input_info:
        if args.verbose:
            print('geo_file=', info[0], ', wavelength=', info[1], ', fibre=', info[2], ', reemis=', info[3], ', abs=', info[4])
        # Make list of commands for job array to call
        commandList_address = jobScript_repo + 'hist_commandList_' + filename_format(info) + '_' + str(i) + '.root'
        commandList_file = open(commandList_address, 'w')
        for i in range(n_evts):
            # Create all the histogram making commands
            hist_command = hist_command_base + save_sims_folder + 'simOut_' + filename_format(info) + '_' + str(i) + '.root '\
                                             + save_splithists_folder + 'splitHist_' + filename_format(info) + '_' + str(i) + '.root '\
                                             + info[2] + ' ' + info[1] + ' ' + str(int(args.verbose))
            commandList_file.write(hist_command + '\n')
        commandList_file.close()

        # Create the job script to run all these macros in an array
        new_job_address = makeJobArrayScript('hists_' + filename_format(info), example_jobScript, save_tothists_folder, commandList_address, info, args.verbose)
        job_addresses.append(new_job_address)

    if len(n_evts) != len(job_addresses):
        print('ERROR: Number of job scripts does not match number needed to run macros')
        print('len(n_evts) = ', len(n_evts), ', len(job_addresses) = ', len(job_addresses))
        exit()


    ### RUN JOB SCRIPTS ###
    print('Submitting jobs...')
    for i in range(len(n_evts)):
        command = 'qsub -t 1-' + str(n_evts[i]) + ' ' + job_addresses[i] 
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    # Wait until these job arrays are done
    checkJobsDone('hists_jobScript_', filename_format(input_info), 10)

    ### COMBINE SPLIT UP HISTS ###
    print('Combining hists...')
    combi_command_base = 'hadd tot_hists_'
    current_wd = os.getcwd()
    os.chdir(save_splithists_folder)
    for info in input_info:
        command = combi_command_base + filename_format(info) + '.root hists_' + filename_format(info) + '_*.root'
        if args.verbose:
            print('Running command: ', command)
        subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    return True


### Analysis functions ###

def getSlopes(args, input_info):
    '''Compute slope from FOM computed from different absorption scalings. If verbose flag is
    true, also records FOM for each different abs file. So, if only an FOM is wanted, just run this
    and ignore the slope result.'''

    # Read in example macro and job script + info
    repo_address = __file__[:-len('scripts/runAMMELIE.py')]

    jobScript_address = repo_address + 'job_scripts/jobSingle.job'
    with open(jobScript_address, "r") as f:
        example_jobScript = f.readlines()

    # Make sure folders are of the correct format to  use later
    save_tothists_folder = checkRepo(args.tothist_repo, args.verbose)
    save_stats_folder = checkRepo(args.stats_repo, args.verbose)

    # Folder for job scripts and command lists they use
    jobScript_repo = save_tothists_folder + 'job_scripts/'
    jobScript_repo = checkRepo(jobScript_repo, args.verbose)

    ### Check that all simulations have finished running ###
    # checkJobsDone('sims_jobScript_', filename_format(info), 10)

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

    ### MAKE JOB SCRIPTS TO RUN ANALYSIS ###
    print('Creating macros and job scripts...')

    # Check all input info that should be the same is in fact the same (everything except absorption for the moment)
    # geo_files = input_info[:, 0], wavelengths = input_info[:, 1], fibres = input_info[:, 2], reemissions = input_info[:, 3], abs_factors = input_info[:, 4]
    result = np.all(input_info[:, 0] == input_info[0, 0]) and np.all(input_info[:, 1] == input_info[0, 1]) and np.all(input_info[:, 2] == input_info[0, 2]) and np.all(input_info[:, 3] == input_info[0, 3])

    slope_command_base = repo_address + 'scripts/./GetFOMabsSlope.exe '
    

    # Create command
    info_str = input_info[0, 0] + '_' + input_info[0, 1] + '_' + input_info[0, 2] + '_' + input_info[0, 3]  # Same as usual but without absorption
    output_stats_file = save_stats_folder + 'slopeStats_' + info_str + '.txt'
    output_root_file = save_stats_folder + 'Regions_' + info_str + '.root'
    slope_command = slope_command_base + args.list_file + ' ' + save_tothists_folder + ' '\
                    + output_stats_file + ' ' + str(int(args.verbose)) + ' ' + output_root_file + ' ' + region_lims

    # Create the job script to run all these macros in an array
    new_job_address = makeJobSingleScript('slopes_' + info_str, example_jobScript, save_tothists_folder, slope_command, input_info[0], args.verbose)

    ### RUN JOB SCRIPTS ###
    print('Submitting job...')
    command = 'qsub ' + new_job_address
    if args.verbose:
        print('Running command: ', command)
    subprocess.call(command, stdout=subprocess.PIPE, shell=True) # use subprocess to make code wait until it has finished

    # Wait until these job arrays are done
    checkJobsDone('slopes_', filename_format(input_info), 10)

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
        table['FOMs']
        for i in range(2, len(stats)):
            info = stats[i].split(' ')  # [abs, FOM]
            table[info[0]] = float(info[1])

    # Save table to json file
    save_file = save_stats_folder + 'FinalStats_' + info_str + '.json'
    with open(save_file, 'w') as f:
        json.dump(table, f)

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