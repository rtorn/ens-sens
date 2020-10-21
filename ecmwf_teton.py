import os
import time
import shutil
import sys
import cfgrib
import datetime as dt
import glob
import pandas as pd
import numpy as np
import xarray as xr

def stage_grib_files(datea, config):

    '''
    This is a generic class for copying or linking grib file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.
 
    This particular instance is for ECMWF data on the machine teton at UAlbany.  In this 
    case, the grib files are linked to files that are located in another directory on that
    machine.

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        config  (dict):  The dictionary with configuration information
    '''

    freq = config.get('fcst_hour_int', 12)
    fmax = config.get('fcst_hour_max', 120)

    #  Make the work directory if it does not exist
    if not os.path.isdir(config['work_dir']):
       try:
          os.makedirs(config['work_dir'])
       except OSError as e:
          raise e

    init   = dt.datetime.strptime(datea, '%Y%m%d%H')
    init_s = init.strftime("%m%d%H%M")

    #  Loop over all forecast times, link to the source file
    for fhr in range(0, int(fmax)+int(freq), int(freq)):

       datef   = init + dt.timedelta(hours=fhr)
       datef_s = datef.strftime("%m%d%H%M")

       grib_file = "E1E{0}{1}1".format(str(init_s), str(datef_s))
       infile    = config['model_dir'] + '/' + grib_file

       #  Only try to copy if the file is not there
       if ( not os.path.isfile(config['work_dir'] + '/' + grib_file) ):

          #  Wait for the source file to be present 
          while not os.path.exists(infile):
             time.sleep(20.1)

          #  Wait for the file to be finished being copied
          while ( (time.time() - os.path.getmtime(infile)) < 60 ):
             time.sleep(10)

          try:  #  Try to link from the source to the work directory
             os.symlink(infile, config['work_dir'] + '/' + grib_file)
          except Exception as err:
             print(err)


def stage_atcf_files(datea, bbnnyyyy, config):
    '''
    This is a generic class for copying or linking ATCF file data from a specific location
    to a directory where calculations are performed.  No matter where these calculations are
    carried out, this routine must exist.

    The result is a set of ATCF files in the work directory of the format atcf_NN.dat, 
    where NN is the ensemble member number.

    This particular instance is for ECMWF data on the machine teton at UAlbany.  In this 
    case, all ATCF data for a particular storm is in one file, so the code waits for this
    initialization time to exist, then uses sed to get the lines attributed to each 
    ensemble member and places that data in a seperate file.

    Attributes:
        datea (string):  The initialization time of the forecast (yyyymmddhh)
        config  (dict):  The dictionary with configuration information
    '''
    src  = config['atcf_dir'] + "/a" + bbnnyyyy + ".dat"
    nens = int(config['num_ens'])

    #  Wait for the source file to be present 
    while not os.path.exists(src):
       time.sleep(20.5)

    #  Wait for the ensemble ATCF information to be placed in the file
    while ( len(os.popen("sed -ne /" + datea + "/p " + src + " | sed -ne /EE/p").read()) == 0 ):
       time.sleep(20.7)

    #  Wait for the file to be finished being copied
    while ( (time.time() - os.path.getmtime(src)) < 60 ):
       time.sleep(10)

    for n in range(nens + 1):

       nn = '%0.2i' % n
       file_name = config['work_dir'] + "/atcf_" + nn + ".dat"

       #  If the specific member's ATCF file does not exist, copy from the source file with sed.
       if not os.path.isfile(file_name):

          fo = open(file_name,"w")
          fo.write(os.popen("sed -ne /" + datea + "/p " + src + " | sed -ne /EE" + nn + "/p").read())
          fo.close()


#  Class to read information from ensemble grib files
class ReadGribFiles:
    def __init__(self, datea, fhr, config):

        self.datea = datea
        self.src_path = config['work_dir']
        self.datea = dt.datetime.strptime(datea, '%Y%m%d%H')
        self.datea_s = self.datea.strftime("%m%d%H%M")
        self.datea_1 = self.datea + dt.timedelta(hours=fhr)
        self.datea_1 = self.datea_1.strftime("%m%d%H%M")
        self.grib_dict = self.__read_grib_dict()

        self.var_dict = {'zonal_wind': 'u', 'meridional_wind': 'v', 'geopotential_height': 'gh', 'temperature': 't', 
                         'relative_humidity': 'r', 'sea_level_pressure': 'msl'}

    ###  Function that creates a dictionary for each ensemble gribfile
    def __read_grib_dict(self):
        grib_file = "E1E{0}{1}1".format(str(self.datea_s), str(self.datea_1))
        file_name = os.path.join(self.src_path, grib_file)
        if self.__check_file_exists(file_name):
            ds = cfgrib.open_datasets(file_name)
            ds_dict = {}
            for d in ds:
                for tt in d:
                    if 'number' in d[tt].dims:
                        ds_dict.update({'{0}_pf'.format(tt): d[tt]})
                    else:
                        ds_dict.update({'{0}_cf'.format(tt): d[tt]})

            return ds_dict
        else:
            exit(2)

    ###  Function to read an ensemble array of fields based on input information in vdict
    def create_ens_array(self, varname, nens, vdict):

       vname = self.var_dict[varname] + '_cf'

       #  Determine latitude bounds.  If values are in vdict, use, otherwise, use entire array
       if 'latitude' in vdict:

          if float(self.grib_dict[vname].attrs['GRIB_latitudeOfFirstGridPointInDegrees']) > float(self.grib_dict[vname].attrs['GRIB_latitudeOfLastGridPointInDegrees']):
             slat1 = int(vdict['latitude'][1])
             slat2 = int(vdict['latitude'][0])
          else:
             slat1 = int(vdict['latitude'][0])
             slat2 = int(vdict['latitude'][1])

       else:

          latvec = list(self.grib_dict[vname].latitude.data)
          slat1  = latvec[0]
          slat2  = latvec[-1]

       #  Determine longitude bounds.  If values are in vdict, use, otherwise, use entire array
       if 'longitude' in vdict:

          lon1 = int(vdict['longitude'][0])
          lon2 = int(vdict['longitude'][1])

       else:

          lonvec = list(self.grib_dict[vname].longitude.data)
          lon1   = lonvec[0]
          lon2   = lonvec[-1]

       #  Create attributes based on what is in the file
       attrlist = {}
       if 'description' in vdict:
         attrlist['description'] = vdict['description']
       if 'units' in vdict:
         attrlist['units'] = vdict['units']
       if '_FillValue' in vdict:
         attrlist['_FillValue'] = vdict['_FillValue']

       #  Create a dummy array that can be used to copy data into
       lonvec = list(self.grib_dict[vname].sel(longitude=slice(lon1, lon2)).longitude.data)
       latvec = list(self.grib_dict[vname].sel(latitude=slice(slat1, slat2)).latitude.data)
       ensarr = xr.DataArray(name='ensemble_data', data=np.zeros([nens, len(latvec), len(lonvec)]),
                             dims=['ensemble', 'latitude', 'longitude'], attrs=attrlist, 
                             coords={'ensemble': [i for i in range(nens)], 'latitude': latvec, 'longitude': lonvec}) 

       return(ensarr)

    #  Function to read a single ensemble member's forecast field
    def read_grib_field(self, varname, member, vdict):

       vname = self.var_dict[varname] + '_cf'

       if 'latitude' in vdict:

          if float(self.grib_dict[vname].attrs['GRIB_latitudeOfFirstGridPointInDegrees']) > float(self.grib_dict[vname].attrs['GRIB_latitudeOfLastGridPointInDegrees']):
             slat1 = int(vdict['latitude'][1])
             slat2 = int(vdict['latitude'][0])
          else:
             slat1 = int(vdict['latitude'][0])
             slat2 = int(vdict['latitude'][1])

       else:

          latvec = list(self.grib_dict[vname].latitude.data)
          slat1  = latvec[0]
          slat2  = latvec[-1]

       if 'longitude' in vdict:

          lon1 = int(vdict['longitude'][0])
          lon2 = int(vdict['longitude'][1])

       else:

          lonvec = list(self.grib_dict[vname].longitude.data)
          lon1   = lonvec[0]
          lon2   = lonvec[-1]

       #  Read a single pressure level of data, if this is a variable that has pressure levels
       if 'isobaricInhPa' in vdict:

          if self.grib_dict[vname].isobaricInhPa[0] > self.grib_dict[vname].isobaricInhPa[1]:
            slev1 = int(vdict['isobaricInhPa'][1])
            slev2 = int(vdict['isobaricInhPa'][0])
          else:
            slev1 = int(vdict['isobaricInhPa'][0])
            slev2 = int(vdict['isobaricInhPa'][1])
 

          if member == 0:
             vname = self.var_dict[varname] + '_cf'
             vout  = self.grib_dict[vname].sel(latitude=slice(slat1, slat2), longitude=slice(lon1, lon2), isobaricInhPa=slice(slev1, slev2))

          else:
             vname = self.var_dict[varname] + '_pf'
             vout  = self.grib_dict[vname].sel(number=member, latitude=slice(slat1, slat2), longitude=slice(lon1, lon2), 
                                                isobaricInhPa=slice(slev1, slev2))

       #  Read the only level if it is a single level variable
       else:

          if member == 0:
             vname = self.var_dict[varname] + '_cf'
             vout  = (self.grib_dict[vname].sel(latitude=slice(slat1, slat2), longitude=slice(lon1, lon2)))

          else:

             vname = self.var_dict[varname] + '_pf'
             vout  = (self.grib_dict[vname].sel(number=member, latitude=slice(slat1, slat2), longitude=slice(lon1, lon2)))

       return(vout)


    @staticmethod
    def __check_file_exists(filename):
        isfile = False
        try:
            if os.path.isfile(filename):
                isfile = True
        except Exception as err:
            print(err)

        return isfile

if __name__ == '__main__':

    src1 = "/Users/parthpatwari/RA_Atmospheric_Science/Old_Code/atcf_data"
    grib_src = "/Users/parthpatwari/RA_Atmospheric_Science/GRIB_files"
    dest1 = "/Users/parthpatwari/RA_Atmospheric_Science/New_Code/atcf_data"
    atcf_src = "/Users/parthpatwari/RA_Atmospheric_Science/Old_Code/atcf_data"
    # c1 = CopyFiles(src, dest)
    # if c1.checkandcreatedir():
    #     c1.copy_filestowork()
    g1 = ReadGribFiles(grib_src, '2019082900', 180)
    print(g1.grib_dict)
    a1 = Readatcfdata(atcf_src)
    print(a1.atcf_array)
