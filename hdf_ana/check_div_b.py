#imports
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
import sys
import hdf_class as hc

data_num=sys.argv[1]
DATA=hc.HDF5_DATA(data_num=data_num)
fig,ax=plt.subplots(1,1)
DATA.check_div_b(ax)
print 'saving'
plt.savefig("./test.png",dpi=200)
