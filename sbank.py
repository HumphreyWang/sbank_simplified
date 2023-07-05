"""
This is the main part of stochastic bank generation, we use Parser to parse input parameters, you can type
```shell
python3 sbank.py --help
```
in command line for further help docs.

NOTE: It'd better set an environment variable: (Do this before import numpy!)
>>>import os
>>>os.environ["OMP_NUM_THREADS"] = "1"  # a better way than os.system("export OMP_NUM_THREADS=1")
In pycharm you can set OMP_NUM_THREADS=1 in [run/debug configurations]
(Run - Edit configurations... or in the upper right toolbar)
"""
from collections import deque
from optparse import OptionParser
from time import strftime
import os

import numpy as np

from bank import Bank
from templates import coord_frames, proposals

import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection


def parse_command_line():
    parser = OptionParser()
    # coord_frames parameter options
    parser.add_option("--coord-frame", choices=list(coord_frames.keys()), metavar='|'.join(coord_frames.keys()),
                      help="Required. Specify the coord_frames to use for template generation.")
    parser.add_option("--x1-min", type="float", metavar="FLOAT",
                      help="Required. Set minimum x of the first axis.")
    parser.add_option("--x1-max", type="float", metavar="FLOAT",
                      help="Required. Set maximum x of the first axis.")
    parser.add_option("--x2-min", type="float", metavar="FLOAT",
                      help="Set minimum x of the second axis. If not specified, the x limits provided on the first axis will be assumed for the second axis.")
    parser.add_option("--x2-max", type="float", metavar="FLOAT",
                      help="Set maximum x of the second axis. If not specified, the x limits provided on the first axis will be assumed for the second axis.")
    # initial condition options
    parser.add_option("--seed", type="int", metavar="INT", default=42,
                      help="Set the seed for the random number generator used for parameter(x1, x2) generation.")
    parser.add_option("--bank-seed", metavar="FILE", action="append", default=[],
                      help="Add templates from FILE to the initial bank. Can be specified multiple times. Only the additional templates will be outputted.")
    # distance calculation options
    parser.add_option("--distance-max", type="float", default=0.1,
                      help="Set maximum distance of the bank. Note that since this is a stochastic process, the requested maximal distance may not be strictly guaranteed but should be fulfilled on a statistical basis. Default: 0.1.")
    parser.add_option("--convergence-threshold", type="int", metavar="N", default=1000,
                      help="Set the criterion for convergence of the stochastic bank. The code terminates when there are N rejected proposals for each accepted proposal, averaged over the last ten acceptances. Default 1000.")
    parser.add_option("--max-new-templates", type="int", metavar="N", default=float('inf'),
                      help="Use this option to force the code to exit after accepting a specified number N of new templates. Note that the code may exit with fewer than N templates if the convergence criterion is met first.")
    parser.add_option("--neighborhood-size", type="float", metavar="N", default=0.25,
                      help="Specify the window size to define \"nearby\" templates used to compute the distance against each proposed template. The neighborhood is chosen symmetric about the proposed template; \"nearby\" is defined using the option --neighborhood-param. The default value of 0.25 is *not a guarantee of performance*. Choosing the neighborhood too small will lead to larger banks (but also higher bank coverage).")
    parser.add_option("--neighborhood-param", choices=["x1", "x2", "norm"], default="x1",
                      help="Choose how the neighborhood is sorted for match calculations.")
    # output options
    parser.add_option("--output-filename", default=None,
                      help="Required. Name for output template bank. May not clash with seed bank.")
    parser.add_option("--generate-full-plots", action="store_true", default=False,
                      help="Generate plots during calculation. DO NOT use it for bank size >~ 40, be sure that you know how much storage space it will take!")
    parser.add_option("--verbose", action="store_true", default=False,
                      help="Be verbose and write diagnostic information out to file.")

    opts_, args_ = parser.parse_args()

    # check for required arguments
    for opt in ("coord_frame", "x1_min", "x1_max", "output_filename"):
        if getattr(opts_, opt) is None:
            parser.error("--%s is required" % opt.replace("_", "-"))

    if opts_.x2_min is None:
        if opts_.coord_frame == 'Polar':
            opts_.x2_min = 0.
        else:
            opts_.x2_min = opts_.x1_min
    if opts_.x2_max is None:
        if opts_.coord_frame == 'Polar':
            opts_.x2_max = 2 * np.pi
        else:
            opts_.x2_max = opts_.x1_max

    if opts_.coord_frame == 'Polar':
        if (not 0 <= opts_.x1_min <= opts_.x1_max) or (not 0 <= opts_.x2_min <= opts_.x2_max <= 2 * np.pi):
            parser.error("For polar coordinates, (x1, x2) denotes for (r, theta), where r>=0, theta in [0, 2*pi].")

    for seed in opts_.bank_seed:
        if seed == opts_.output_filename:
            raise ValueError("Bank seed %s would be overwritten by output file. Choose a different output name." % seed)

    return opts_, args_


print(strftime('%Y/%m/%d %H:%M:%S'))
opts, args = parse_command_line()
if opts.generate_full_plots:
    fig_dir = f'fig_{opts.output_filename[:-4]}'
    os.makedirs(fig_dir, exist_ok=False)
else:
    fig_dir = ''

# choose coord_frame
tmplt_class = coord_frames[opts.coord_frame]

# initialize the bank
bank = Bank(opts.neighborhood_size, opts.neighborhood_param, if_plot=opts.generate_full_plots)

# add templates to bank
for seed_file in opts.bank_seed:
    arr = np.load(seed_file).T
    bank.add_from_array(arr, tmplt_class)
    if opts.verbose:
        print("Added %d seed templates from %s to initial bank." % (len(arr), seed_file))

if opts.verbose:
    print("Initialized the template bank to seed with %d precomputed templates." % len(bank))

np.random.mtrand.seed(opts.seed)

constraints = {'x1': (opts.x1_min, opts.x1_max),
               'x2': (opts.x2_min, opts.x2_max)}
proposal = proposals[opts.coord_frame](tmplt_class, **constraints)

# For robust convergence, ensure that an average of k_max/len(ks) of
# the last len(ks) proposals have been rejected by SBank.
ks = deque(10 * [1], maxlen=10)
k = 0  # k is n_prop per iteration
n_prop = 0  # count total number of proposed templates

# main working loop
for tmplt in proposal:

    # check if stopping criterion has been reached
    if not (((k + float(sum(ks))) / len(ks) < opts.convergence_threshold) and
            (len(bank) < opts.max_new_templates)):
        break
    # accounting for number of proposals
    k += 1  # since last acceptance
    n_prop += 1  # total throughout lifetime of process
    prefix = f'{fig_dir}/{n_prop:0>5d}'

    # check if proposal is already covered by existing templates
    distance, matcher = bank.covers(tmplt, opts.distance_max, prefix)
    if distance > opts.distance_max:
        bank.insort(tmplt, prefix)
        ks.append(k)
        if opts.verbose:
            print("\nbank size: %d\t\tproposed: %d\trejection rate: %.6f / (%.6f)" %
                  (len(bank), n_prop, 1-float(len(ks))/float(sum(ks)), 1-1./opts.convergence_threshold))
            print("accepted:\t\t", tmplt)
            if matcher is not None:
                print("min distance (%.4f):\t" % distance, matcher)
        k = 0

scatter_points = np.array([i.params for i in iter(bank)]).T
np.save(opts.output_filename, scatter_points)
print("\ntotal number of proposed templates: %d" % n_prop)
# print("total number of match calculations: %d" % bank._nmatch)
print("final bank size: %d" % len(bank))
print(strftime('%Y/%m/%d %H:%M:%S'))


fig, ax = plt.subplots(figsize=(6, 6))
w, h = 2*opts.distance_max*np.sqrt(tmplt_class.vals)
ec = EllipseCollection(w, h, tmplt_class.ang, units='xy', color='darkblue', alpha=0.2,
                       offsets=np.column_stack(scatter_points), offset_transform=ax.transData)
ax.add_collection(ec)
ax.scatter(*scatter_points, s=1, c='darkblue', marker='.')
# from matplotlib.patches import Arc
# ax.add_patch(Arc((0, 0), 2, 2, 0, 0, 90, color='navy', fill=False))
ax.set_xlabel('$x_1$', fontsize=24)
ax.set_ylabel('$x_2$', fontsize=24)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)  # TODO: improve these hard-coded settings
ax.tick_params(labelsize=16)
ax.grid(which='both', zorder=1, alpha=0.5)
ax.set_aspect('equal')
plt.tight_layout()
plt.show()
plt.savefig(f'{opts.output_filename[:-4]}.png')
