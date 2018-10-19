import yt
import numpy as np
import hdf_class as hc
import sys

data_num=sys.argv[1]
DATA=hc.HDF5_DATA(data_num=data_num)
den=DATA.den[...,np.newaxis]
vx=DATA.vx[...,np.newaxis]
vy=DATA.vy[...,np.newaxis]
vz=DATA.vz[...,np.newaxis]
data = dict(density = (den, "g/cm**3"),
            vx      = (vx , "cm/s"   ),
            vy      = (vy , "cm/s"   ),
            vz      = (vz , "cm/s"   ),
           )
bbox = np.array([[-300, 300], [-300, 300], [-1, 1]])
ds = yt.load_uniform_grid(data, den.shape, length_unit="AU", bbox=bbox, nprocs=4)
slc = yt.SlicePlot(ds, "z", "density")
slc.annotate_streamlines('vy', 'vx')
#slc.annotate_quiver('vy', 'vx',16)
#slc.annotate_line_integral_convolution('vy', 'vx', lim=(0.4,0.55))
#slc.set_cmap("density", "Blues")
slc.set_width((400,'AU'))
slc.annotate_grids(cmap=None)
slc.save()
