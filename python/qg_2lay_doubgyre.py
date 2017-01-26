import numpy as np
import matplotlib.pyplot as plt
import my_toolbox as my
from netCDF4 import Dataset
import sys
import calc_mitgcm as mitgcm

savefig = False
path_fig = '../pics_nac/'
nnf=0

run = 'run_008'
path_data = '/home/nbruggemann/work/proj_nac_hp/qg_wind_gyre/'+run+'/'
#fname = 'pp_isi_timeave.cdf'
fname = 'pybintocdf.cdf'

it = 500

execfile(path_data+'parameter.m')

A = np.pi*tau0/1000./(beta*Hk_1)

# load grid
f = Dataset(path_data+fname, 'r')
lvars = ['xt', 'yt', 'zt', 'xu', 'yu', 'zu', 't']
for var in lvars:
  exec('%s = f.variables[var][:]'%(var))

# load variables 
lvars = []
lvars += ['psi', 'pv']
for var in lvars:
  if it is None:
    exec('%s = f.variables[var][:]'%(var))
  else:
    exec('%s = f.variables[var][it,...]'%(var))
  #try:
  #  exec('%s = %s.filled(0)'%(var,var))
  #except:
  #  pass
f.close()

# Coriolis parameter
fcort = (f0+beta*yt)[:,np.newaxis]
fcoru = (f0+beta*yu)[:,np.newaxis]
fcent = f0+beta*Ly/2.0

psi1 = psi[0]
psi2 = psi[1]

# 
q = np.ma.zeros((nz,ny,nx))
F1 = f0**2/(gred_2*Hk_1)
F2 = f0**2/(gred_2*Hk_2) 
q[0,:,:] = fcort + F1*(psi2-psi1) - f0**2/(gred_1*Hk_1)*psi1
q[1,:,:] = fcort + F2*(psi1-psi2)
psib = psi1+F1/F2*psi2
qb = fcort + F2*psib

qc = f0+beta*yt[-2]/2.0

# ================================================================================ 
# Here starts plotting
# ================================================================================ 
plt.close('all')
xpt = np.concatenate((xu,[xu[-1]+dx]))
ypt = np.concatenate((yu,[yu[-1]+dy]))
xpu = np.concatenate(([xt[0]-dx], xt))
ypu = np.concatenate(([yt[0]-dx], yt))

 
hca, hcb = my.arrange_axes(4,2, plot_cb=True, sasp=1.0, fig_size_fac=2., \
                            sharex=True, sharey=True, xlabel="", ylabel="")
ii=-1

clim = 'auto' 

for iz in [0,1]:
  ii+=1; ax=hca[ii]; cax=hcb[ii]
  my.shade(xpt/1e3, ypt/1e3, psi[iz,:,:], contfs='auto', clim=clim, ax=ax, cax=cax)
  ax.set_title('psi iz = %d'%iz)
  
  ii+=1; ax=hca[ii]; cax=hcb[ii]
  my.shade(xpt/1e3, ypt/1e3, pv[iz,:,:], contfs='auto', clim=clim, ax=ax, cax=cax)
  ax.set_title('pv iz = %d'%iz)

  ii+=1; ax=hca[ii]; cax=hcb[ii]
  my.shade(xpt/1e3, ypt/1e3, q[iz,:,:], contfs='auto', clim=clim, ax=ax, cax=cax)
  ax.set_title('q iz = %d'%iz)

  ii+=1; ax=hca[ii]; cax=hcb[ii]
 
ii=3; ax=hca[ii]; cax=hcb[ii]
my.shade(xpt/1e3, ypt/1e3, psib, contfs='auto', clim=clim, ax=ax, cax=cax)
ax.set_title('psib')

ii=7; ax=hca[ii]; cax=hcb[ii]
my.shade(xpt/1e3, ypt/1e3, qb, contfs='auto', clim=clim, ax=ax, cax=cax)
ax.set_title('qb')


nnf+=1
if savefig:
  print 'save figure: %s_%02d.pdf' % (__file__.split('/')[-1][:-3], nnf)
  plt.savefig("%s_%02d.pdf" % (path_fig+__file__.split('/')[-1][:-3], nnf))

plt.show()
