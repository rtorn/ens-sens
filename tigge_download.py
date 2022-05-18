from ecmwfapi import ECMWFDataServer
import sys, os
import argparse
import datetime as dt
import configparser
import numpy as np

'''
Program that retrieves forecast fields from a single initialization time from the 
TIGGE database.

From command line:

python run_AR_sens.py -init yyyymmddhh --param paramfile

where:

  -init is the initialization date in yyyymmddhh format
  -storm is TC name, where XXXXXXXNNB, where XXXXXXXX is the name, NN is the number, B is the basin (optional)
  -param is the parameter file path (optional, otherwise goes to default values in example.parm)
'''

#  Read the initialization time and storm from the command line
exp_parser = argparse.ArgumentParser()
exp_parser.add_argument('--init',  action='store', type=str, required=True)
exp_parser.add_argument('--storm', action='store', type=str)
exp_parser.add_argument('--param', action='store', type=str)

args = exp_parser.parse_args()

yyyymmddhh = args.init

#  Read the parameter file
if args.param:
   paramfile = args.param
else:
   paramfile = 'example.parm'

confin = configparser.ConfigParser()
confin.read(paramfile)

config = {}
config.update(confin['model'])
config.update(confin['locations'])

#  Modify work and output directory for specific case/time
if args.storm:
  storm = args.storm
  config['work_dir']   = '{0}/{1}.{2}'.format(config['work_dir'],yyyymmddhh,storm)
else:
  config['work_dir']   = '{0}/{1}'.format(config['work_dir'],yyyymmddhh)

#  Create appropriate directories
if not os.path.isdir(config['work_dir']):
   try:
      os.makedirs(config['work_dir'])
   except OSError as e:
      raise e

#  Set the TIGGE read command values
if not config['tigge_forecast_time']:
   raise "tigge_forecast_time is missing"

daystr = yyyymmddhh[0:4] + '-' + yyyymmddhh[4:6] + '-' + yyyymmddhh[6:8]
hhstr  = yyyymmddhh[8:10] + ":00:00"

if config['model_src'] == 'ECMWF':
  model_str = "ecmf"
  model_num = "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30/31/32/33/34/35/36/37/38/39/40/41/42/43/44/45/46/47/48/49/50"
elif config['model_src'] == 'GEFS':
  model_str = "kwbc"
  model_num = "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20/21/22/23/24/25/26/27/28/29/30"
#  model_num = "1/2/3/4/5/6/7/8/9/10/11/12/13/14/15/16/17/18/19/20"

atmres = config.get('tigge_forecast_grid_space','1.0/1.0')
sfcres = "0.25/0.25"

os.chdir(config['work_dir'])

server = ECMWFDataServer()

#  Retrieve pressure level data
server.retrieve({
    'origin'    : model_str,
    'levelist'  : "200/250/300/500/700/850/925/1000",
    'levtype'   : "pl",
    'expver'    : "prod",
    'parameter' : "130/131/132/133/156",
    'number'    : model_num,
    'dataset'   : "tigge",
    'step'      : config['tigge_forecast_time'],
    'grid'      : atmres,
    'time'      : hhstr,
    'type'      : "pf",
    'date'      : daystr,
    'class'     : "ti",
    'target'    : "tigge_output_pl_pf.grib"
})

server.retrieve({
    'origin'    : model_str,
    'levelist'  : "200/250/300/500/700/850/925/1000",
    'levtype'   : "pl",
    'expver'    : "prod",
    'parameter' : "130/131/132/133/156",
    'dataset'   : "tigge",
    'step'      : config['tigge_forecast_time'],
    'grid'      : atmres,
    'time'      : hhstr,
    'type'      : "cf",
    'date'      : daystr,
    'class'     : "ti",
    'target'    : "tigge_output_pl_cf.grib"
})

'''
#  Retrieve PV level data
server.retrieve({
    'origin'    : model_str,
    'levtype'   : "pv",
    'levelist'  : "2",
    'expver'    : "prod",
    'parameter' : "3/131/132",
    'number'    : model_num,
    'dataset'   : "tigge",
    'step'      : config['tigge_forecast_time'],
    'grid'      : atmres,
    'time'      : hhstr,
    'type'      : "pf",
    'date'      : daystr,
    'class'     : "ti",
    'target'    : "tigge_output_pv_pf.grib"
})

server.retrieve({
    'origin'    : model_str,
    'levtype'   : "pv",
    'levelist'  : "2",
    'expver'    : "prod",
    'parameter' : "3/131/132",
    'dataset'   : "tigge",
    'step'      : config['tigge_forecast_time'],
    'grid'      : atmres,
    'time'      : hhstr,
    'type'      : "cf",
    'date'      : daystr,
    'class'     : "ti",
    'target'    : "tigge_output_pv_cf.grib"
})

'''

#  Retrieve surface data
server.retrieve({
    'origin'    : model_str,
    'levtype'   : "sfc",
    'expver'    : "prod",
    'parameter' : "59/136/151/167/168",
    'number'    : model_num,
    'dataset'   : "tigge",
    'step'      : config['tigge_forecast_time'],
    'grid'      : atmres,
    'time'      : hhstr,
    'type'      : "pf",
    'date'      : daystr,
    'class'     : "ti",
    'target'    : "tigge_output_sfc_pf.grib"
})

server.retrieve({
    'origin'    : model_str,
    'levtype'   : "sfc",
    'expver'    : "prod",
    'parameter' : "59/136/151/167/168",
    'dataset'   : "tigge",
    'step'      : config['tigge_forecast_time'],
    'grid'      : atmres,
    'time'      : hhstr,
    'type'      : "cf",
    'date'      : daystr,
    'class'     : "ti",
    'target'    : "tigge_output_sfc_cf.grib"
})

#  Retrieve higher-resolution surface fields
server.retrieve({
    'origin'    : model_str,
    'levtype'   : "sfc",
    'expver'    : "prod",
    'parameter' : "165/166/228",
    'number'    : model_num,
    'dataset'   : "tigge",
    'step'      : config['tigge_forecast_time'],
    'grid'      : sfcres,
    'time'      : hhstr,
    'type'      : "pf",
    'date'      : daystr,
    'class'     : "ti",
    'target'    : "tigge_output_hrsfc_pf.grib"
})

server.retrieve({
    'origin'    : model_str,
    'levtype'   : "sfc",
    'expver'    : "prod",
    'parameter' : "165/166/228",
    'dataset'   : "tigge",
    'step'      : config['tigge_forecast_time'],
    'grid'      : sfcres,
    'time'      : hhstr,
    'type'      : "cf",
    'date'      : daystr,
    'class'     : "ti",
    'target'    : "tigge_output_hrsfc_cf.grib"
})

ftype = ("pl", "sfc", "hrsfc")

for i in range(len(ftype)):
   os.system("wgrib2 -s tigge_output_{0}_pf.grib >& grib_{0}_pf.out".format(ftype[i]))
   os.system("wgrib2 -s tigge_output_{0}_cf.grib >& grib_{0}_cf.out".format(ftype[i]))

hlist = config['tigge_forecast_time'].split("/")

init   = dt.datetime.strptime(yyyymmddhh, '%Y%m%d%H')
init_s = init.strftime("%m%d%H%M")

#  Loop over all times, create one file per time
for t in range(len(hlist)):

  datef   = init + dt.timedelta(hours=int(hlist[t]))
  datef_s = datef.strftime("%m%d%H%M")

  if int(hlist[t]) == 0:
     timestr = ":anl:"
  else:
     timestr = ":{0} hour fcst:".format(hlist[t])
  hhh = '%0.3i' % int(hlist[t])

  if np.remainder(int(hlist[t]),24) == 0:
     pcpstr = ":0-{0} day acc fcst:".format(round(float(hlist[t]) / 24.))
  else:
     pcpstr = ":0-{0} hour acc fcst:".format(hlist[t])

  if os.path.isfile("f{0}_fields.grb".format(hhh)):
     os.remove("f{0}_fields.grb".format(hhh))

  print(hlist[t],hhh,timestr,pcpstr)

  gribout = 'E1E{0}{1}1'.format(str(init_s), str(datef_s))  

  for i in range(len(ftype)):
     os.system("cat grib_{0}_cf.out | grep \"{1}\" | wgrib2 -fix_ncep -i -append tigge_output_{0}_cf.grib \
                   -grib {2} >& /dev/null".format(ftype[i],timestr,gribout))
     os.system("cat grib_{0}_pf.out | grep \"{1}\" | wgrib2 -fix_ncep -i -append tigge_output_{0}_pf.grib \
                   -grib {2} >& /dev/null".format(ftype[i],timestr,gribout))

  if int(hlist[t]) > 0:
     os.system("cat grib_hrsfc_cf.out | grep \"{0}\" | wgrib2 -fix_ncep -i -append tigge_output_hrsfc_cf.grib \
                   -grib {1} >& /dev/null".format(pcpstr,gribout))
     os.system("cat grib_hrsfc_pf.out | grep \"{0}\" | wgrib2 -fix_ncep -i -append tigge_output_hrsfc_pf.grib \
                   -grib {1} >& /dev/null".format(pcpstr,gribout))

for i in range(len(ftype)):
   os.remove("tigge_output_{0}_pf.grib".format(ftype[i]))
   os.remove("tigge_output_{0}_cf.grib".format(ftype[i]))
   os.remove("grib_{0}_pf.out".format(ftype[i]))
   os.remove("grib_{0}_cf.out".format(ftype[i]))

