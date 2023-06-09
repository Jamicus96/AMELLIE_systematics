import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import argparse
import json
import os


def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Plot shit')

    parser.add_argument('--stats_repo', '-r', type=str, dest='stats_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/Sim/Slopes/final_stats/',
                        help='Folder json stats files are saved in.')
    parser.add_argument('--json_file', '-f', type=str, dest='json_file',
                        default='FinalStats_AMELLIE_snoplusnative_labppo_1p1_berkeley_scintillator_LED403_FA108_reemis1.0.json',
                        #default='FinalStats_AMELLIE_snoplusnative_te__LED403_FA108_reemis0.7.json',
                        help='json stats file.')
    parser.add_argument('--plot_repo', '-p', type=str, dest='plot_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/Sim/Slopes/final_stats/plots/',
                        help='json stats file.')

    args = parser.parse_args()
    return args


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

def getSlope(x, y):
    '''Get slope of y vs x data where y(1)=1, using weighted least squares fitting.'''
    # S_xy = ((x - 1.0) * (y - 1.0)) / (y_err**2)
    # S_xx = ((x - 1.0) * (y - 1.0)) / (y_err**2)
    S_xy = np.sum((x - 1.0) * (y - 1.0))
    S_xx = np.sum((x - 1.0)**2)
    return S_xy / S_xx

def plot(data, save_address):
    '''Plot the DATA!!'''

    # Upack data
    title = data['geo_file'] + ', ' + data['inner_av_material'] + ', ' + data['LED'] + ', ' + data['fibre']

    region_lims = data['regions'][0]['region_lims']
    stats = data['regions'][0]['stats']

    abs_scalings = []
    tot_hits = []
    direct_hits = []
    reflected_hits = []
    for abs in stats:
        recorded_info = stats[abs]
        abs_scalings.append(float(abs))
        tot_hits.append(float(recorded_info['tot_hits']))
        direct_hits.append(float(recorded_info['direct_hits']))
        reflected_hits.append(float(recorded_info['reflected_hits']))

        if float(abs) == 1.0:
            nominal_tot_hits = float(recorded_info['tot_hits'])
            nominal_direct_hits = float(recorded_info['direct_hits'])
            nominal_reflected_hits = float(recorded_info['reflected_hits'])

    abs_scalings = np.asarray(abs_scalings)
    tot_hits = np.asarray(tot_hits) / nominal_tot_hits
    direct_hits = np.asarray(direct_hits) / nominal_direct_hits
    reflected_hits = np.asarray(reflected_hits) / nominal_reflected_hits
    FOM =  reflected_hits / direct_hits

    slope = getSlope(abs_scalings, FOM)
    print('Slope = ', slope)

    # Plot data
    plt.figure(figsize=(10, 7), dpi=100)

    plt.scatter(abs_scalings, tot_hits, marker='+', label='tot_hits')
    plt.scatter(abs_scalings, direct_hits, marker='+', label='direct_hits')
    plt.scatter(abs_scalings, reflected_hits, marker='+', label='reflected_hits')
    plt.scatter(abs_scalings, FOM, marker='x', label='FOM (reflected_hits/direct_hits)')

    plt.title(title + ' => Slope = %.2f' % slope)
    plt.xlabel('Absorption coefficient scaling')
    plt.ylabel('Number of hits normalised at abs=1\n(nominal: tot_hits=' + str(nominal_tot_hits) + ', direct_hits='\
                + str(nominal_direct_hits) + ', reflected_hits=' + str(nominal_reflected_hits))
    #plt.ylim([0.95, 1.06])
    plt.legend(loc='best')
    plt.xscale('log')
    plt.yscale('log')

    plt.savefig(save_address, dpi=1000)
    #plt.show()


### MAIN ###

def main():
    # read in arguments
    args = argparser()
    stats_repo = checkRepo(args.stats_repo, True)
    plot_repo = checkRepo(args.plot_repo, True)

    # Load data
    f = open(stats_repo + args.json_file)
    data = json.load(f)

    # plot shit
    plot_saveFile_name = plot_repo + args.json_file.replace('.json', '') + '.png'
    plot(data, plot_saveFile_name)

if __name__ == '__main__':
    main()