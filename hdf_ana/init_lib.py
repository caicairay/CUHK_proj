import numpy as np
import pickle
import hdf_class as hc

data_num=0
DATA=hc.HDF5_DATA(data_num=data_num)

def init_ene():
    Ek,Eb=DATA.calc_ene()
    Ek_avr=DATA.cir_avr(np.log10(Ek))
    Eb_avr=DATA.cir_avr(np.log10(Eb))
    Ekl,Ekr=np.split(Ek_avr,2) 
    Ebl,Ebr=np.split(Eb_avr,2) 
    return Ekr,Ebr
def init_v():
    vr,vt=DATA.calc_v()
    vr_avr=DATA.cir_avr(np.log10(vr))
    vt_avr=DATA.cir_avr(np.log10(vt))
    vrl,vrr=np.split(vr_avr,2) 
    vtl,vtr=np.split(vt_avr,2) 
    return vrr,vtr

vr,vt=init_v()
f=open("../init_data/init_vr_log.txt","w+")
pickle.dump(vr,f)
f.close()
f=open("../init_data/init_vt_log.txt","w+")
pickle.dump(vt,f)
f.close()

ek,eb=init_ene()
f=open("../init_data/init_ek_log.txt","w+")
pickle.dump(ek,f)
f.close()
f=open("../init_data/init_eb_log.txt","w+")
pickle.dump(eb,f)
f.close()

