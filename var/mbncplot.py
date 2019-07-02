import sys
import argparse

# DEBUG
# import pdb

import netCDF4 as nc
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(\
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Plot variables in MBDyn output using Matplotlib.',\
        epilog='Multiple --var arguments can be provided.\n' + \
                'Components should be indicated by a single integer and are 0-based.\n' + \
                'For matrices, indexes 0-8 can be used, indicating the components in.\n' + \
                'row-major ordering\n' + \
                'Please note that it is not mandatory to place the --comps option directly \n' + \
                'after the related variables, but variables and components will be regarded \n' + \
                'as ordered lists. In other works, these two inputs are equivalent: \n\n' + \
                'python mbncplot.py --var node.struct.1.X --comps 1 --var node.struct.2.X --comps 2\n' + \
                'python mbncplot.py --var node.struct.1.X --var node.struct.2.X --comps 1 --comps 2\n'\
                )

parser.add_argument('ncfile', metavar='ncfile', help='MBDyn NetCDF output file')
parser.add_argument('--time', '-t', action='store_true',\
        help='use simulation time on X axis')
parser.add_argument('--var', '-v', nargs=1, \
        help='variable[s] to be plotted', action='append', required=True)
parser.add_argument('--comps', '-c', type=int, \
        help='component[s] of the variable to be plotted', \
        nargs=1, action='append')

args = parser.parse_args()

# row-major indexes
rows = {0: 0, 1: 0, 2: 0, 3: 1, 4: 1, 5: 1, 6: 2, 7: 2, 8: 2}
cols = {0: 0, 1: 1, 2: 2, 3: 0, 4: 1, 5: 2, 6: 0, 7: 1, 8: 2}

try:
    nd = nc.Dataset(args.ncfile, 'r')
except FileNotFoundError:
    print('Error: NetCDF file not found.')
    sys.exit(1)

# fig = plt.figure()
if args.time:
    X = nd.variables['time']
    plt.xlabel('Time [s]')
else:
    X = nd.variables['run.step']
    plt.xlabel('Step')

try:
    for idx in range(len(args.var)):
        ncVar = nd.variables[args.var[idx][0]]
        if len(ncVar.shape) == 1:
            plt.plot(X, ncVar, label = args.var[idx][0])
        elif len(ncVar.shape) == 2:
            if args.comps[idx]:
                plt.plot(X, ncVar[:, args.comps[idx][0]], \
                        label = args.var[idx][0] + '[:,' + str(args.comps[idx][0]) + ']'\
                        + '  [' + ncVar.units +']')
            else:
                for jdx in range(ncVar.shape[1]):
                    plt.plot(X, ncVar[:, jdx], \
                            label = args.var[idx][0] + '[:, ' + str(jdx) + ']'\
                            + '  [' + ncVar.units +']')
        elif len(ncVar.shape) == 3:
            if args.comps[idx][0]:
                plt.plot(X, ncVar[:, rows[args.comps[idx][0]], cols[args.comps[idx][0]]],\
                        label = args.var[idx][0] + '[:, ' + \
                        str[rows[args.comps[idx][0]]] + ' , ' + 
                        str[cols[args.comps[idx][0]]] + ']'\
                        + '  [' + ncVar.units +']')
            else:
                for jdx in range(9):
                    plt.plot(X, ncVar[:, rows[jdx], cols[jdx]],\
                        label = args.var[idx][0] + '[:, ' + \
                        str[rows[jdx]] + ' , ' + 
                        str[cols[jdx]] + ']'\
                        + '  [' + ncVar.units +']')
except KeyError:
    print('Error: variable ' + args.var[idx][0] + ' not found')
    sys.exit(2)

plt.legend()
plt.show()
