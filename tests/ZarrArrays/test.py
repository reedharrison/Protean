import os as os
import numpy as np
import zarr as zarr
import multiprocessing as mp
import itertools as it

n_threads = 4
datafile = 'array.zarr'

shape = (1000, 1000)
chunks = (100, 100)
synchronizer = zarr.ProcessSynchronizer('array.sync')

a = zarr.open(datafile, mode='a', shape=shape, chunks=chunks, dtype=float, synchronizer=synchronizer)
a[0, :] = np.arange(1000)
a[:, 0] = np.arange(1000)


def f(args):
	(i, j), x = args
	return (i, j), x+x

pool = mp.Pool(processes=n_threads)
indices = it.product(f)



r = zarr.open(datafile, mode='r', shape=shape, chunks=chunks, dtype=float, synchronizer=synchronizer)





from joblib import Parallel, delayed
import numpy as np

if __name__ is '__main__':
	result = Parallel(n_jobs=2)(delayed(np.sum)([i, j]) for i in range(1, 6) for j in range(11, 16))