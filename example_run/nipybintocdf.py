#!/usr/bin/env python
# ================================================================================ 
#       input:
# ================================================================================ 
def user_input(par, contlist):
  """
  for help see "help" below
  """
  #par.path_data    = "/scratch/uni/ifmto/u241161/smm/th_smm_47/run_000/"
  #par.contlist_names = ['pv', 'psi', 'u', 'v', 'hpr', 'wek']
  par.contlist_names = ['pv', 'psi', 'u', 'v', 'hpr', 'pvr', 'pvp', 'Gpvadv', 'Gpvdif', 'Gpvfor']
  par.level_by_level = False
  par.myqg = True
  par.mask_land = False

  """
  contlist.append( binfile  ( 'Tbal', nrecords=5) )
  contlist[-1].set_names    ( ['TOTTTEND', 'ADVr_TH ', 'ADVx_TH ', 'ADVy_TH ', 'DFrI_TH '] )
  contlist[-1].set_longnames( 'temp tend', 'vertical adv temp flux', 'zonal adv temp flux', 'meridional adv temp flux', 'implicite diffusive flux')
  contlist[-1].set_grids    ( [['xt', 'yt', 'zt'], ['xt', 'yt', 'zu'], ['xu', 'yt', 'zt'], ['xt', 'yu', 'zt'], ['xt', 'yt', 'zu']] )
  contlist[-1].set_units    ( [ 'degree C / sec', 'degree C m^2 / sec', 'degree C m^2 / sec', 'degree C m^2 / sec', 'degree C m^2 / sec'])

  contlist.append( binfile  ( 'T' ) )
  contlist[-1].set_names    ( ['temp'] )
  contlist[-1].set_units    ( ['degree C'] )
  contlist[-1].set_grids    ( [['xt', 'yt', 'zt']] )
  contlist[-1].set_longnames( ['temperature'] )

  contlist.append( binfile  ( 'U' ) )
  contlist[-1].set_units    ( ['m / sec'] )
  contlist[-1].set_grids    ( [['xu', 'yt', 'zt']] )
  contlist[-1].set_longnames( ['zonal velocity'] )

  contlist.append( binfile  ( 'V' ) )
  contlist[-1].set_units    ( ['m / sec'] )
  contlist[-1].set_grids    ( [['xt', 'yu', 'zt']] )
  contlist[-1].set_longnames( ['meridional velocity'] )

  contlist.append( binfile  ( 'W' ) )
  contlist[-1].set_units    ( ['m / sec'] )
  contlist[-1].set_grids    ( [['xt', 'yt', 'zu']] )
  contlist[-1].set_longnames( ['vertical velocity'] )

  contlist.append( binfile  ( 'Eta' ) )
  contlist[-1].set_units    ( ['m'] )
  contlist[-1].set_grids    ( [['xt', 'yt']] )
  contlist[-1].set_longnames( ['SSH'] )
  """
  return(par, contlist)
# ================================================================================ 
#       end input
# ================================================================================ 

# ===============================================================================
#       help
# ===============================================================================
  """
    Define here all parameters to convert the MITgcm binary data into netcdf 4 data.
  Parameters and their defaults:
    par.path_data     = "pwd"                           # path to binary data 
    par.path_outdata  = "path_data"                     # path where netcdf file is saved
    par.contlist_names       = ["T"]                    # list of file prefixes                       
    par.netcdfname    = "pybintocdf.cdf"                # name of netcdf file
    par.deltaT        = "auto"                          # time step [s] to infer on model time: tstep*deltaT
    par.tsteps        = "auto"                          # list of tsteps can be chosen explicitly 
    par.indt          = "all"                           # integer that says that every indt'th time step is chosen ("all"=1)
    par.binformat     = ">f"                            # > big and < little endian; f for single and d for double precission
    par.mask_land     = True                            # if True: values of data that are exactly zero are masked
    par.lonlat        = False                           # if True: name of x and y variables are degree lon and lat otherwise m
    par.fm            = "not to be set by user"         # handle of netcdf file; do not set this parameter
    par.dims          = "not to be set by user"         # dimension of file (nt, nz, ny, nx)

  - Minimum example of parameters that have to be set:
      par.path_data   = "/scratch/uni/ifmto/u241161/smm/th_smm_47/run_000"
      par.netcdfname  = "pybintopdf.cdf"
      par.contlist_names     = ['T', 'Eta']
    
  - More extensive example:
      par.path_data   = "/scratch/uni/ifmto/u241161/smm/th_smm_47/run_000"
      par.netcdfname  = "pybintopdf.cdf"
      par.contlist_names     = ['T', 'Eta']
      par.tsteps      = np.arange(50, 100, 10)*1000     # start from tstep 50000 go in steps of 10000 up to tstep 100000

      contlist.append( binfile  ( 'T' ) )
      contlist[-1].set_names    ( ['temp'] )
      contlist[-1].set_units    ( ['degree C'] )
      contlist[-1].set_grids    ( [['xt', 'yt', 'zt']] )

      contlist.append( binfile  ( 'Eta' ) )
      contlist[-1].set_units    ( ['m'] )
      contlist[-1].set_grids    ( [['xt', 'yt']] )
      contlist[-1].set_longnames( ['SSH'] )
      contlist[-1].set_longnames( ['temperature'] )

  - With multiple records in each file:
      par.contlist_names = ['Tbal']

      contlist.append( binfile  ( 'Tbal', nrecords=5) )
      contlist[-1].set_names    ( ['TOTTTEND', 'ADVr_TH ', 'ADVx_TH ', 'ADVy_TH ', 'DFrI_TH '] )
      contlist[-1].set_longnames( 'temp tend', 'vertical adv temp flux', 'zonal adv temp flux', 'meridional adv temp flux', 'implicite diffusive flux')
      contlist[-1].set_grids    ( [['xt', 'yt', 'zt'], ['xt', 'yt', 'zu'], ['xu', 'yt', 'zt'], ['xt', 'yu', 'zt'], ['xt', 'yt', 'zu']] )
      contlist[-1].set_units    ( [ 'degree C / sec', 'degree C m^2 / sec', 'degree C m^2 / sec', 'degree C m^2 / sec', 'degree C m^2 / sec'])
  """


# -------------------------------------------------------------------------------- 
#       myout
# -------------------------------------------------------------------------------- 
def myout(string):
  print(string)
  return

# -------------------------------------------------------------------------------- 
#       class binfile
# -------------------------------------------------------------------------------- 
class binfile(object):
  def __init__(self, fprfx, path_data, names="auto", nrecords="auto", dims="auto", units="", grids="auto", longnames=""):

    self.set_fprfx(fprfx)

    # find one meta file
    flist = glob.glob(path_data + fprfx + ".*.meta")  
    if len(flist) == 0 and (names=="auto" or nrecords=="auto" or dims=="auto"):
      sys.exit(":: Error: Could not find any meta file but it is needed since names, nrecords and dims have not been specified! ::")
    else:
      fname = flist[0]
      meta_names, meta_dims, meta_nrecords = read_meta(fname)

    if names == "auto":
      self.set_names(meta_names)

    if dims == "auto":
      self.dims = meta_dims

    if nrecords == "auto":
      self.nrecords = meta_nrecords[0]

    if longnames == "":
      longnames  = [""]*self.nrecords 
    self.set_longnames(longnames)
    
    if units == "" or len(units)==1:
      units = [units]*self.nrecords 
    self.set_units(units)

    if grids == "auto":
      grids = self.set_auto_grids(self.names) 
    self.set_grids(grids)

  def set_auto_grids(self, names):
    grids = []
    for var in names: 
      if self.dims[0]==1:
        grids.append(['xt', 'yt'])
      elif var[0] == "T":
        grids.append(['xt', 'yt', 'zt'])
      elif var[0] == "U":
        grids.append(['xu', 'yt', 'zt'])
      elif var[0] == "V":
        grids.append(['xt', 'yu', 'zt'])
      elif var[0] == "W":
        grids.append(['xt', 'yt', 'zu'])
      else:
        grids.append(['xt', 'yt', 'zt'])
    return(grids)

  def set_fprfx(self, fprfx):
    self.fprfx = fprfx
    return

  def set_names(self, names):
    for nn, name in enumerate(names):
      # change t or T to temp
      if name == 't' or name == 'T':
        myout("Change variable name \"%s\" to \"temp\" to satisfy Ferret." % (name))
        names[nn] = 'temp'
      # replace spaces
      names[nn] = name.replace(" ", "")
    self.names = names
    return

  def set_units(self, units):
    if len(units)==1:
      units = units*self.nrecords
    self.units = units
    return

  def set_grids(self, grids):
    if len(grids)==1:
      grids = grids*self.nrecords
    for nn, grid in enumerate(grids):
      ngrid = ['t']
      grid.sort(reverse=True)
      for gg in grid:
        ngrid.append(gg)
      grids[nn] = tuple(ngrid)
    self.grids = grids
    return

  def set_longnames(self, longnames):
    if len(longnames)==1:
      longnames = longnames*self.nrecords
    self.longnames = longnames
    return

# -------------------------------------------------------------------------------- 
#       pseudo parameter class
# -------------------------------------------------------------------------------- 
class parameter(object):
  """ 
  pseudo parameter class 
  """
  def __init__(self):
    self.path_data     = "pwd"                           # path to binary data 
    self.path_outdata  = "path_data"                     # path where netcdf file is saved
    self.varlist       = []                              # list of file prefixes                       
    self.netcdfname    = "pybintocdf.cdf"                # name of netcdf file
    self.deltaT        = "auto"                          # time step [s] to infer on model time: tstep*deltaT
    self.tsteps        = "auto"                          # list of tsteps can be chosen explicitly 
    self.indt          = "all"                           # integer that says that every indt'th time step is chosen ("all"=1)
    self.level_by_level= True                            # if True z-levels are read sepearatly (better for large files)
    self.binformat     = ">f"                            # > big and < little endian; f for single and d for double precission
    self.mask_land     = True                            # if True: values of data that are exactly zero are masked
    self.lonlat        = False                           # if True: name of x and y variables are degree lon and lat otherwise m
    self.fm            = "not to be set by user"         # handle of netcdf file; do not set this parameter
    self.dims          = [0,0,0]                         # dimension of file (nz, ny, nx)
    self.myqg          = False                           # decides if data comes from MITgcm or myqg


# -------------------------------------------------------------------------------- 
#       eval_user_input
# -------------------------------------------------------------------------------- 
def eval_user_input(par, contlist):
  # check if path names are valid
  # if path_data is not given
  if par.path_data == "pwd":
    par.path_data = subprocess.check_output(["pwd"])[:-1]

  # if path_outdata is not given make it to be path_data
  if par.path_outdata == "path_data": 
    par.path_outdata = par.path_data

  # check path names
  if par.path_data[-1] != "/":
    par.path_data = par.path_data + "/"
  if par.path_outdata[-1] != "/":
    par.path_outdata = par.path_outdata + "/"

  # initilize binfile-objects in contlist if not done so far
  if len(par.contlist_names) == 0 and len(contlist) == 0:
    sys.exit(":: Error: Either set par.contlist_names or defince contlist! ::")
  elif len(contlist) == 0 and len(par.contlist_names)!=0:
    par.contlist = []
    for cont in par.contlist_names:
      # initialize objects of binfile-class
      contlist.append( binfile  ( cont, par.path_data ) )

  # initialize time steps if not done so far
  if par.tsteps == "auto": 
    pref = contlist[0].fprfx
    flist = glob.glob(par.path_data + pref + ".*.meta")  
    if len(flist) == 0:
      sys.exit(":: Error: Could not find any file that matches %s! ::" % (par.path_data + pref + ".*.meta"))
    par.tsteps = np.zeros((len(flist)), dtype=int)
    for l, fname in enumerate(flist):
      par.tsteps[l] = int(fname[-15:-5])
    par.tsteps = np.sort(par.tsteps)

  # to allow par.tsteps to be lists or np.arrays in the input; in the following par.tsteps are np.arrays
  if type(par.tsteps) == list:
    par.tsteps = np.array(par.tsteps)

  # take only indt values of tsteps
  if isinstance(par.indt, int):
    par.tsteps = par.tsteps[::par.indt]
  elif isinstance(par.indt, list):
    par.indt = [ it for it in par.indt if it < par.tsteps.size ]
    par.tsteps = par.tsteps[indt]
  elif par.indt == "all":
    pass
  else:
    sys.exit(":: Error: Wrong value in indt! ::")
  
  # find deltaT
  if par.deltaT == "auto":
    if os.path.isfile(par.path_data + "STDOUT.0000"):
      with open(par.path_data + "STDOUT.0000") as fobj:
        lnum = -10
        for ll, line in enumerate(fobj):
          if "deltaTClock" in line:
            lnum = ll
          if ll == lnum + 1:
            par.deltaT = float(line[21:]) 
            break
    else:
      par.deltaT = 10.
      myout("Attention could not find STDOUT.0000 and no deltaT was specified so dummy deltaT=10. is used!")

  # get grid information
  varlist, xdims, nrecords = read_meta(par.path_data+"XC.meta")
  varlist, ydims, nrecords = read_meta(par.path_data+"YC.meta")
  varlist, zdims, nrecords = read_meta(par.path_data+"RC.meta")
  par.dims = np.array([0]*3)
  par.dims[0] = zdims[0]
  par.dims[1] = ydims[1]
  par.dims[2] = xdims[2]

  # strings for printing some basic information
  if par.lonlat:
    gridname = "Polar"
  else:
    gridname = "Cartesian"

  if par.tsteps.size > 1:
    tstepstr = "{:d}: {:d} : {:d}".format(par.tsteps[0], par.tsteps[1]-par.tsteps[0], par.tsteps[-1]) 
  elif par.tsteps.size == 1:
    tstepstr = "{:d}".format(par.tsteps[0])

 # print some basic information
  myout(" path_data:      {:s}".format(par.path_data))
  myout(" netcdf file:    {:s}".format(par.path_outdata+par.netcdfname))
  myout(" binfiles:       {:s}".format(str(par.contlist_names)))
  myout(" tsteps:         {:s}".format(tstepstr))
  myout(" num tsteps:     {:d}".format(par.tsteps.size))
  myout(" model dim:      {:s}".format(str(par.dims)))
  myout(" grid:           {:s}".format(gridname))
  myout(" level_by_level: {:s}".format(str(par.level_by_level)))
  myout(" deltaT:         {:f} sec".format(par.deltaT))
  myout(" mask land:      {:s}".format(str(par.mask_land)))

  return(par, contlist)

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

  # get name of records
  if "fldList" in locals():
    fldList = fldList[0]
    varlist = []
    nlen = 8
    for nn in range(len(fldList)/nlen):
      varlist.append( (fldList[nn*nlen:(nn+1)*nlen]).replace(" ","") )
  else:
    varlist = fname.split("/")[-1]
    varlist = [varlist.split(".")[0]]

  # get dimensions
  dimList = np.array(dimList).reshape(nDims[0],-1)
  nx = dimList[0,0]
  ny = dimList[1,0]
  nz = dimList[2,0] if nDims[0]==3 else 1 

  if fname.split("/")[-1] == "RF.meta":
    nz = nz - 1

  dims = np.array([nz, ny, nx])

  """
  if nDims[0] == 3:
    dims = np.array([nz, ny, nx])
  elif nDims[0] == 2:
    dims = np.array([ny, nx])
  else:
    print ":: Error: nDims = %d! I cannot deal with that!" % (nDims[0]) 
    sys.exit(1)
  """

  return(varlist, dims, nrecords)
 
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


def read_data_old(fname, numrec=0, numlev=0, dims='meta', binformat=">f"):

  if dims=='meta':
    varlist, dims, nrecords = read_meta(fname)
   
  if binformat[1] == "f":
    nbytes = 4
  elif binformat[1] == "d":
    nbytes = 8

  # if vertical profile or horizontal section
  if dims[1] == 1 and dims[2] == 1:
    nnums = dims[0]
    resht = (dims[0])
  else:
    nnums = dims[1]*dims[2]
    resht = (dims[1], dims[2])

  # read data file
  fobj = open(fname,'r')
  fobj.seek( nbytes*nnums*(numlev+numrec*dims[0]) )
  data = fobj.read( nbytes*nnums )
  fobj.close()

  # unpack data
  data = np.array( struct.unpack( (binformat[0]+"%i"+binformat[1])%(nnums), data ) )
  # reshape data
  data = np.array(data).reshape( resht )

  return(data)

# -------------------------------------------------------------------------------- 
#       get_grid
# -------------------------------------------------------------------------------- 
def get_grid(path_data, binformat=">f"):
  
  if par.myqg:
    xt=read_data(path_data+'XC.data', par.dims[2], binformat=binformat)
    xu=read_data(path_data+'XG.data', par.dims[2], binformat=binformat)
    yt=read_data(path_data+'YC.data', par.dims[1], binformat=binformat)
    yu=read_data(path_data+'YG.data', par.dims[1], binformat=binformat)
  else:
    xt=read_data(path_data+'XC.data', par.dims[1:], binformat=binformat)
    xt=xt[0,:]
    xu=read_data(path_data+'XG.data', par.dims[1:], binformat=binformat)
    xu=xu[0,:]

    yt=read_data(path_data+'YC.data', par.dims[1:], binformat=binformat)
    yt=yt[:,0]
    yu=read_data(path_data+'YG.data', par.dims[1:], binformat=binformat)
    yu=yu[:,0]

  zt=read_data(path_data+'RC.data', par.dims[0], binformat=binformat)
  zt=zt[:]
  zu=read_data(path_data+'RF.data', par.dims[0], binformat=binformat)
  zu=zu[:]
  
  return(xt, xu, yt, yu, zt, zu)


# -------------------------------------------------------------------------------- 
#       measure_time
# -------------------------------------------------------------------------------- 
def measure_time(ct, STEP = "unknown"):
  ct.append(datetime.datetime.now())
  if len(ct) == 1:
    print "Start timer:"
  else:
    tdiff = (ct[-1]-ct[-2]).seconds
    print "Step % 10s took % 10.0f sec" % (STEP, tdiff)
  return ct

def write_grid(par, contlist):
 # close netcdf file (necessary only in interactive mode if errors occur)
  try:
    par.fm.close()
  except:
    pass

  # open netcdf file (overwrite potentially existing old file)
  par.fm = Dataset(par.path_outdata + par.netcdfname,'w', format='NETCDF4')
  
  # units for x and y
  if par.lonlat:
    xdim_unit = 'lon'
    ydim_unit = 'lat'
  else:
    xdim_unit = 'm'
    ydim_unit = 'm'

  # load grid data
  xt, xu, yt, yu, zt, zu = get_grid(par.path_data, binformat=par.binformat)
  t = par.tsteps*par.deltaT

  nx = xt.size
  ny = yt.size
  nz = zt.size
  nt = t.size

  # save grid data to netcdf
  par.fm.createDimension('xt', nx)
  xtdim       = par.fm.createVariable('xt','f4',('xt',))
  xtdim.units = xdim_unit
  xtdim[:]    = xt

  par.fm.createDimension('yt', ny)
  ytdim       = par.fm.createVariable('yt','f4',('yt',))
  ytdim.units = ydim_unit
  ytdim[:]    = yt

  par.fm.createDimension('zt', nz)
  ztdim       = par.fm.createVariable('zt','f4',('zt',))
  ztdim.units = 'm'
  ztdim[:]    = zt

  par.fm.createDimension('xu', nx)
  xudim       = par.fm.createVariable('xu','f4',('xu',))
  xudim.units = xdim_unit
  xudim[:]    = xu

  par.fm.createDimension('yu', ny)
  yudim       = par.fm.createVariable('yu','f4',('yu',))
  yudim.units = ydim_unit
  yudim[:]    = yu

  par.fm.createDimension('zu', nz)
  zudim       = par.fm.createVariable('zu','f4',('zu',))
  zudim.units = 'm'
  zudim[:]    = zu

  par.fm.createDimension('t', nt)
  tdim        = par.fm.createVariable('t','f4',('t',))
  tdim.units  = 'sec'
  tdim[:]     = t

  return(par, contlist)

# -------------------------------------------------------------------------------- 
#       main_loop
# -------------------------------------------------------------------------------- 
def main_loop(par, contlist):
# --- loop over all file container ---
  for nn, cobj in enumerate(contlist):
    # find out dimension of files
    fname     = "%s/%s.%010d.data" % (par.path_data, cobj.fprfx, par.tsteps[0])
    cobj.names, dims, nrecords = read_meta(fname)

    myout("--------------------------------------------------------------------------------")
    myout(" Processing file container no. {:2d}/{:2d}:        {:s}.xxxxxxxxxx.data ".format(nn+1, len(contlist), cobj.fprfx))
    myout(" Dimension ( nx x ny x nz ):                {:4d} x {:4d} x {:4d}".format(dims[2], dims[1], dims[0]))
    #myout(" First loaded file is                        {0:s}.{1:010d}.data ".format(cobj.fprfx, par.tsteps[0]))
    #myout(" Last loaded file is                         {0:s}.{1:010d}.data ".format(cobj.fprfx, par.tsteps[-1]))

# --- loop over all variables in file container ---
    ncv = [0]*len(cobj.names)
    for nv, var in enumerate(cobj.names):
      myout("----")
      #myout(" Processing variable no. {:2d}/{:2d}:              {:s} -- {:20s} ".format(nv+1, len(cobj.names), var, cobj.longnames[nv]))
      myout(" Processing variable no. {:2d}/{:2d}:              {:s} ".format(nv+1, len(cobj.names), var))
      myout(" Unit:                                       {:20s} ".format(cobj.units[nv]))
      myout(" Grid:                                       {:s}".format(str(cobj.grids[nv])))

      # create variable and add meta data
      cobj.names[nv] = cobj.names[nv].replace(" ","")
      ncv[nv]           = par.fm.createVariable(cobj.names[nv], 'f4', cobj.grids[nv], fill_value=-9e33)
      ncv[nv].units     = cobj.units[nv]
      ncv[nv].long_name = cobj.longnames[nv]
      
      #cti = measure_time([])
      #cto = measure_time([]) 

# --- loop over all time steps ---
      for l, tstep in enumerate(par.tsteps):
        sys.stdout.write( "\r Prossesed time steps {:3d} / {:3d}".format(l+1, par.tsteps.size) )
        sys.stdout.flush()
        fname     = "%s/%s.%010d.data" % (par.path_data, cobj.fprfx, tstep)

        if dims[0] == 1:
          numskip  = dims.prod()*nv
          var_data = read_data(fname, dims, numskip=numskip, binformat=par.binformat)
          if par.mask_land:
            var_data  = np.ma.array(var_data, mask=var_data==0.0)
          ncv[nv][l,:,:]  = var_data[0,:,:]

        elif not par.level_by_level:
          numskip  = dims.prod()*nv
          var_data = read_data(fname, dims, numskip=numskip, binformat=par.binformat)
          if par.mask_land:
            var_data  = np.ma.array(var_data, mask=var_data==0.0)
          ncv[nv][l,:,:,:]  = var_data[:,:,:]

        else:
# --- loop over all vertical levels ---
          for k in range(dims[0]):
            #var_data  = read_data(fname, numrec=nv, numlev=k, binformat=par.binformat)
            numskip  = par.dims[1:].prod()*(par.dims[0]*nv+k)
            var_data = read_data(fname, par.dims[1:], numskip=numskip, binformat=par.binformat)
            if par.mask_land:
              var_data  = np.ma.array(var_data, mask=var_data==0.0)
            ncv[nv][l,k,:,:]  = var_data[:,:]
      myout("")

  par.fm.close()
  return(par, contlist)

# --------------------------------------------------------------------------------
#     Execute the programm
# --------------------------------------------------------------------------------
myout("--------------------------------------------------------------------------------")
myout("{:^80s}".format("Execute "+__file__))
myout("--------------------------------------------------------------------------------")

import numpy as np
from netCDF4 import Dataset
import struct
import datetime
import sys, os, subprocess, glob

starttime = datetime.datetime.now()  

contlist = []
par   = parameter()

par, contlist  = user_input(par, contlist)
par, contlist  = eval_user_input(par, contlist)
par, contlist  = write_grid(par, contlist)
par, contlist  = main_loop(par, contlist)

endtime = datetime.datetime.now()
myout("--------------------------------------------------------------------------------")
myout("du -sh:      " + subprocess.check_output(["du","-sh", par.path_outdata+par.netcdfname])[:-1])
myout("started at:      {:s}".format( starttime.isoformat() ))
myout("finished at:     {:s}".format( endtime.isoformat() ))
myout("conversion took: {:f} sec = {:f} min".format( (endtime-starttime).seconds, (endtime-starttime).seconds/60. ))
myout("--------------------------------------------------------------------------------") 

