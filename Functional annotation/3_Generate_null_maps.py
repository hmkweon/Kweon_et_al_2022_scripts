
# module load Python/3.9.5-GCCcore-10.3.0-bare
# PYTHONPATH=/sw/arch/Centos8/EB_production/2021/software/Python/3.9.5-GCCcore-10.3.0-bare/easybuild/python:/home/hkweon/.local/lib/python3.9/site-packages:/home/hkweon/tools/brainsmash:/home/hkweon/conda/envs/brainsmash/lib/python3.9/site-packages:/home/hkweon/conda/envs/brainsmash/lib:/home/hkweon/conda/lib/python3.8/site-packages

# import os
import numpy
import matplotlib.pyplot as plt
from brainsmash.workbench.geo import volume
from brainsmash.mapgen.eval import sampled_fit
from brainsmash.mapgen.sampled import Sampled
from joblib import Parallel, delayed
import multiprocessing

output_dir = "./"
coord_file = "./COORD_TO_PERM.txt"

# Compute distance: need to be done once
filenames = volume(coord_file, output_dir)

# These are three of the key parameters affecting the variogram fit
kwargs = {'ns': 1500,
          'knn': 800,
          'pv': 25,
          }

def run_func(i):
    input = "./" + i + "_VAL_TO_PERM.txt"
    gen = Sampled(x=input, D="./distmat.npy", index="./index.npy", **kwargs, resample=True)
    surrogate_maps = gen(n=10000)
    file = "./NULL_" + i + "_MAPS.txt"
    numpy.savetxt(file, surrogate_maps)

LIST = ["SES", "PGS", "COND_PGS"]
Parallel(n_jobs=3)(delayed(run_func)(i) for i in LIST)

