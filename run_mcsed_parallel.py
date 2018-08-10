""" Script for running MCSED in parallel

.. moduleauthor:: Greg Zeimann <gregz@astro.as.utexas.edu>

"""

import numpy as np
import os
import sys
from astropy.table import Table, vstack
from multiprocessing import cpu_count, Manager, Process
from run_mcsed_fit import main as run_mcsed_ind
from run_mcsed_fit import parse_args
from ssp import read_ssp
from distutils.dir_util import mkpath


def worker(f, i, chunk, ssp_info, out_q, err_q, kwargs):
    ''' Simple design to catch exceptions from the given call '''
    try:
        result = f(argv=chunk, ssp_info=ssp_info)
    except Exception as e:
        err_q.put(e)
        return

    # output the result and task ID to output queue
    out_q.put((i, result))


def parallel_map(func, argv, args, ncpu, ssp_info, clean=False, **kwargs):
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
        Built argsuments from argv
    ncpu : int
        Number of parallelized cpus
    ssp_info : list
        SSP data for spectra, ages, metallicities, ect.
    '''
    if isinstance(ncpu, (int, np.integer)) and ncpu == 1:
        return [func(0, argv, **kwargs)]

    manager = Manager()
    out_q = manager.Queue()
    err_q = manager.Queue()
    jobs = []

    if '-t' in argv:
        x = np.arange(args.nobjects)
        v = [len(i) for i in np.array_split(x, ncpu)]
        counts = np.cumsum(v, dtype=int) - v[0]
        chunks = [argv + ['--nobjects', '%i' % vi, '--count', '%i' % cnt]
                  for vi, cnt in zip(v, counts)]
    else:
        mkpath('temp')
        ind = argv.index('-f')
        data = Table.read(argv[ind+1], format='ascii')
        datachunks = np.array_split(data, ncpu)
        for i, chunk in enumerate(datachunks):
            T = Table(chunk)
            T.write('temp/temp_%i.dat' % i, format='ascii', overwrite=True)
        chunks = [argv + ['-f', 'temp/temp_%i.dat' % i]
                  for i, chunk in enumerate(datachunks)]

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

    if clean:
        for i, chunk in enumerate(chunks):
            if os.path.exists('temp/temp_%i.dat' % i):
                os.remove('temp/temp_%i.dat' % i)

    return results


def main_parallel(argv=None):
    argv = sys.argv
    argv.remove('run_mcsed_parallel.py')
    argv = argv + ['--parallel', 'True']
    args = parse_args(argv=argv)
    ssp_info = None  # read_ssp(args)
    NCPU = cpu_count()
    ncpu = np.max([1, NCPU - 2])
    results = parallel_map(run_mcsed_ind, argv, args, ncpu, ssp_info)
    table = vstack([result[0] for result in results])
    table.write('output/%s' % args.output_filename,
                format='ascii.fixed_width_two_line',
                formats=results[0][1], overwrite=True)

if __name__ == '__main__':
    main_parallel()
