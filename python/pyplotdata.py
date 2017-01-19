#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys, os, struct
# -------------------------------------------------------------------------------- 
#       read_data
# -------------------------------------------------------------------------------- 
def read_data(fname, dims, numskip=0, binformat=">f"):
  if binformat[1] == "f":
    nbytes = 4
  elif binformat[1] == "d":
    nbytes = 8
  
  nnums = dims.prod()

  # read data file
  fobj = open(fname,'r')
  fobj.seek( nbytes*numskip )
  data = fobj.read( nbytes*nnums )
  fobj.close()

  # unpack data
  data = np.array( struct.unpack( (binformat[0]+"%i"+binformat[1])%(nnums), data ) )
  # reshape data
  data = np.array(data).reshape( dims )

  return(data)

# -------------------------------------------------------------------------------- 
#       read_meta
# -------------------------------------------------------------------------------- 
def read_meta(fname):
  # change from .data to .meta
  fname = fname[:-4]+"meta" 

  # check if file exists
  if not os.path.isfile(fname):
    sys.exit(":: Error file %s does not exist! ::" % (fname))

  # read metafile and execute metafile
  fobj = open(fname,'r')
  metafile = fobj.read()
  metafile = metafile.replace("{","[")
  metafile = metafile.replace("}","]")
  metafile = "if True:\n"+metafile
  fobj.close()
  exec(metafile)

  dims = np.array(dims)

  return(dims)

# -------------------------------------------------------------------------------- 
#       main part
# -------------------------------------------------------------------------------- 
fname = "test.data"

"""
nx = 10
ny = 12
nz = 2
dims = np.array([nz, ny, nx])
"""

dims = read_meta(fname)
data = read_data(fname, dims, binformat=">f")
dataplt = data[0,:,:]

plt.close("all")
ii = -1
hca = [0]

ii = ii+1
hca[ii] = plt.axes()
hp = hca[ii].pcolormesh(dataplt)
plt.colorbar(hp)

plt.show()
