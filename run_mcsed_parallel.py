""" Script for running MCSED in parallel

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np
import os
import sys
from astropy.table import Table, vstack
from multiprocessing import cpu_count, Manager, Process
#from run_mcsed_fit import main as run_mcsed_ind
#from run_mcsed_fit import parse_args
from ssp import read_ssp
from distutils.dir_util import mkpath
import run_mcsed_fit
run_mcsed_ind = run_mcsed_fit.main
parse_args = run_mcsed_fit.parse_args

def worker(f, i, chunk, ssp_info, out_q, err_q, kwargs):
    ''' Simple design to catch exceptions from the given call '''
    try:
        result = f(argv=chunk, ssp_info=ssp_info)
    except Exception as e:
        err_q.put(e)
        return

    # output the result and task ID to output queue
    out_q.put((i, result))


def parallel_map(func, argv, args, ncpu, ssp_info, clean=True, **kwargs):
    '''
    Make multiple calls to run_mcsed_fit's main function for either
    test or real data to parallelize the computing effort.  Collect the info
    at the end.

    Inputs
    ------
    func : callable function
        This is meant to be run_mcsed_fit's main function
    argv : list
        Arguments (command line or otherwise) list.
        python run_mcsed_fit.py -h
    args : class
        Built arguments from argv
    ncpu : int
        Number of parallelized cpus
    ssp_info : list
        SSP data for spectra, ages, metallicities, etc.
    clean : bool
        Remove temporary files
    '''
    if isinstance(ncpu, (int, np.integer)) and ncpu == 1:
        return [func(0, argv, **kwargs)]

    manager = Manager()
    out_q = manager.Queue()
    err_q = manager.Queue()
    jobs = []

    if args.test:
        ncpu = min(args.nobjects, ncpu)
        x = np.arange(args.nobjects)
        v = [len(i) for i in np.array_split(x, ncpu)]
        counts = np.full(len(v), 1)
        counts[1:] += np.cumsum(v)[:-1]
        chunks = [argv + ['--nobjects', '%i' % vi, '--count', '%i' % cnt, '--already_parallel']
                  for vi, cnt in zip(v, counts)]
    else:
        mkpath('temp')
        ind = argv.index('-f')
        data = Table.read(argv[ind+1], format='ascii')
        # split up the input file into NCPU files and save as temporary files
        ncpu = min(len(data), ncpu)
        datachunks = np.array_split(data, ncpu)
        for i, chunk in enumerate(datachunks):
            T = Table(chunk)
            T.write('temp/temp_%i.dat' % i, format='ascii', overwrite=True)
        # create a separate argument list for each temporary file
        chunks = [argv + ['-f', 'temp/temp_%i.dat' % i, '--already_parallel']
                  for i, chunk in enumerate(datachunks)]
## WPBWPB delete:
#        print(chunks)
#        return

    for i, chunk in enumerate(chunks):
        p = Process(target=worker, args=(func, i, chunk, ssp_info, out_q,
                                         err_q, kwargs))
        jobs.append(p)
        p.start()

    # gather the results
    for proc in jobs:
        proc.join()

    if not err_q.empty():
        # kill all on any exception from any one slave
        raise err_q.get()

    # Processes finish in arbitrary order. Process IDs double
    # as index in the resultant array.
    results = [None] * len(jobs)
    while not out_q.empty():
        idx, result = out_q.get()
        results[idx] = result

    # Remove the temporary (divided) input files
    if (clean) & (not args.test):
        for i, chunk in enumerate(chunks):
            if os.path.exists('temp/temp_%i.dat' % i):
                os.remove('temp/temp_%i.dat' % i)
        try:
            os.rmdir('temp/')
        except OSError:
            pass
    return results


def main_parallel(argv=None):

## WPBWPB delete
#    print('this is argv from parallel.main_parallel:')
#    print(argv)
##    print(type(argv))
#    return

    # read command line arguments, if not already calling
    # from within run_mcsed_fit.py
    if argv == None:
        argv = sys.argv
        argv.remove('run_mcsed_parallel.py')
        argv = argv + ['--parallel'] 

    args = parse_args(argv=argv)


## WPBWPB delete
#    print('this is argv:')
#    print(argv)
#    print(type(argv))
#    args = parse_args(argv=argv)
#    print(vars(args).keys())
#    print(args)
#    print(type(argv))

    ssp_info = None  # read_ssp(args)
    NCPU = cpu_count()
    ncpu = np.max([1, NCPU - args.reserved_cores])
## WPBWPB delete
#    print('this is args.reserved_cores:  ' + str(args.reserved_cores))
#    print('this is ncpu: ' + str(ncpu))
#    print('this is argv from parallel.main_parallel:')
#    print(argv)
#    print(type(argv))
#    return
## WPBWPB delete
#    return
    results = parallel_map(run_mcsed_ind, argv, args, ncpu, ssp_info)
    table = vstack([result[0] for result in results])
    if args.output_dict['parameters']:
        table.write('output/%s' % args.output_filename,
                    format='ascii.fixed_width_two_line',
                    formats=results[0][1], overwrite=True)
    if args.output_dict['settings']:
        filename = open('output/%s.args' % args.output_filename, 'w')
        filename.write( str( vars(args) ) )
        filename.close()


if __name__ == '__main__':
    main_parallel()
