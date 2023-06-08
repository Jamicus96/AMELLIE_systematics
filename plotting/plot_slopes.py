import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import argparse


def argparser():
    '''Define arguments'''

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Plot shit')

    parser.add_argument('--slopes_file', '-r', type=str, dest='slopes_file',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/2p2gL/text_slopes_out.txt',
                        help='File with slopes stats saved in.')
    parser.add_argument('--plot_repo', '-p', type=str, dest='plot_repo',
                        default='/mnt/lustre/projects/epp/general/neutrino/jp643/rat_dependent/AMELLIE/2p2gL/',
                        help='repo to save output plots in.')

    args = parser.parse_args()
    return args

def create_lims(min, max, N):
    step = (max - min) / float(N - 1)
    # step = (max - min) / float(N)

    lims = np.zeros(N)
    for n in range(N):
        lims[n] = min + n * step

    return lims

def plot_ax(ax, vals, x, y, x_name, y_name, minimum, maximum):

    # Get axis ticks
    x_ticks = []
    for x_val in x:
        x_ticks.append('%.2f' % x_val)
    y_ticks = []
    for y_val in y:
        y_ticks.append('%.2f' % y_val)

    ax.set_xlabel(x_name)
    ax.set_xticks(list(range(len(x)))) # values
    ax.set_xticklabels(x_ticks, rotation=45) # labels

    ax.set_ylabel(y_name)
    ax.set_yticks(list(range(len(y)))) # values
    ax.set_yticklabels(y_ticks) # labels

    # Plot
    im = ax.imshow(vals, vmin = minimum, vmax = maximum)
    return im


def plotting(lims, max_idx, slopes, args):

    # Set up figure
    fig = plt.figure()
    fig, ax = plt.subplots(1, 3, figsize=(17, 6))
    title_1 = 'FOM vs Absorption Slope, with Different Cut Limits\n'
    title_3 = '\ndirect x max = %.2f, direct t centre = %.2f, direct t width = %.2f, reflected x min = %.2f, reflected t centre = %.2f, reflected t width = %.2f'\
             % (lims[0, max_idx[0]], lims[1, max_idx[1]], lims[2, max_idx[2]], lims[3, max_idx[3]], lims[4, max_idx[4]], lims[5, max_idx[5]])
    fig.suptitle(title_1 + r'centered on (x = cos$\theta$):' + title_3)

    # Get data
    data_30 = slopes[:, max_idx[1], max_idx[2], :, max_idx[4], max_idx[5]]
    data_21 = slopes[max_idx[0], :, :, max_idx[3], max_idx[4], max_idx[5]]
    data_54 = slopes[max_idx[0], max_idx[1], max_idx[2], max_idx[3], :, :]

    # Find Maximum and Minimum
    maximum = np.amax(slopes)
    minimum = np.amin(np.array([data_30, data_21, data_54]))

    # Iterate over plots
    im1 = plot_ax(ax[0], data_30, lims[3, :], lims[0, :], 'reflected x min', 'direct x max', minimum, maximum)
    im2 = plot_ax(ax[1], data_21, lims[2, :], lims[1, :], 'direct t width', 'direct t centre', minimum, maximum)
    im3 = plot_ax(ax[2], data_54, lims[5, :], lims[4, :], 'reflected t width', 'reflected t centre', minimum, maximum)

    # Make legend bar
    fig.subplots_adjust(right=0.9, left=0.05)
    cbar_ax = fig.add_axes([0.91, 0.23, 0.03, 0.5])
    colbar = fig.colorbar(im3, cax=cbar_ax)
    colbar.ax.set_ylabel('       slope', rotation=0)

    # Save it
    plt.savefig(args.plot_repo + 'test_plot.png')


def main():
    args = argparser()

    # Create reagion lims to iterate over (6 values):
    dir_x_max_min = -0.95;  dir_x_max_max = -0.8
    ref_x_min_min = 0.8;    ref_x_min_max = 0.95
    dir_t_cen_min = -7.0;  dir_t_cen_max = 15.0
    ref_t_cen_min = -7.0;  ref_t_cen_max = 15.0
    dir_t_wid_min = 3.0;    dir_t_wid_max = 20.0
    ref_t_wid_min = 3.0;    ref_t_wid_max = 20.0
    lims_min_max = [(dir_x_max_min, dir_x_max_max), (dir_t_cen_min, dir_t_cen_max), (dir_t_wid_min, dir_t_wid_max),
                    (ref_x_min_min, ref_x_min_max), (ref_t_cen_min, ref_t_cen_max), (ref_t_wid_min, ref_t_wid_max)]
    
    N = 10
    lims = np.full((6, N), -999.0)
    for i, lim_lims in enumerate(lims_min_max):
        min, max = lim_lims
        lims[i] = create_lims(min, max, N)

    print(lims)

    # Get info
    slopes = np.full((N, N, N, N, N, N), -99.0)
    indices = np.full(6, -1)
    with open(args.slopes_file) as file:
        for line in file:
            line_splt = np.array(line.split(' '), dtype=float)  # dir_x_max_lims dir_t_cen_lims dir_t_wid_lims ref_x_min_lims ref_t_cen_lims ref_t_wid_lims slope
            lims_temp = line_splt[:-1]
            slope = line_splt[-1]

            for i, lim in enumerate(lims_temp):
                indices[i] = np.argmin(abs(lims[i] - lim))
            
            slopes[indices[0], indices[1], indices[2], indices[3], indices[4], indices[5]] = slope

    lim_names = ['dir_x_max_lims', 'dir_t_cen_lims', 'dir_t_wid_lims', 'ref_x_min_lims', 'ref_t_cen_lims', 'ref_t_wid_lims']

    max_slope = np.max(slopes)
    max_idx = np.argwhere(slopes == max_slope)[0]

    print('Max slope = {}, for:'.format(max_slope))
    for i, lim in enumerate(lim_names):
        print(lim + ' = ' + str(lims[i, max_idx[i]]))

    RIP = np.argwhere(slopes == -99.0)
    print('RIP:', len(RIP))

    plotting(lims, max_idx, slopes, args)


if __name__ == '__main__':
    main()