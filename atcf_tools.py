import os
import sys
import datetime as dt
import glob
import pandas as pd
import numpy as np

class ReadATCFData:
    '''
    Class that reads ATCF-format data from an ensemble of files (assumes that each file contains one
    ensemble member at one initialization time.  The class itself contains a series of dictionaries that
    contain the TC track and intensity information.  The data is saved both with pandas and using a 
    dictionary.  The data can be accessed through a series of routines.

    This class is currently under development, especially the best track portions

    Attributes:
        infiles (string): list of ATCF ensemble member files

    '''

    def __init__(self, infiles):

        self.cols = ['basin', 'tcnum', 'datea', 'mm', 'ftype', 'fhr', 'lat', 'lon', 'wnd', 'mslp', \
                     'stype', 'rval', 'ord', 'rad1', 'rad2', 'rad3', 'rad4', 'a']
        self.atcf_files = {}
        self.atcf_array = {}
        self.no_atcf_files = 0
        self.missing = -9999.
        atcf_data = glob.glob(infiles)
        atcf_data = sorted(atcf_data)
        self.no_atcf_files = len(atcf_data)
        for f in range(self.no_atcf_files):

            #  Read entire file both through pandas and ascii reading routines
            file_name = "df_{0}".format(str((f + 1000))[1:])
            file_s = file_name
            self.atcf_files.update({file_s: file_name})
            self.atcf_array.update({f: {}})
            with open(atcf_data[f], 'r') as fo:
                data = np.array(fo.readlines())
                total_rows = data.shape[0]
                if total_rows == 0:
#                    if data[0][0:7] == "missing":
                  n_lines = 0
#                else:
#                  n_lines = total_rows
#                  n_cols = len(data[0])
                else:
                   n_lines = total_rows
                   n_cols = len(data[0])
            fo.close()

            if n_lines > 0:
              file_name = pd.read_csv(filepath_or_buffer=atcf_data[f], header=None)
              file_name.columns = self.cols[0:len(file_name.columns)]

            prfhr    = -1
            ntimes   = -1
            fhrvec   = []
            fmivec   = []
            latvec   = []
            lonvec   = []
            slpvec   = []
            wndvec   = []
            r34nevec = []
            r34sevec = []
            r34swvec = []
            r34nwvec = []

            #  Loop over lines in the ATCF file
            for line in range(n_lines):
                self.atcf_array.get(f).update({line: [0 for i in range(11)]})
                data[line] = data[line].strip()

                #  Parse data and create entries that have all information for each file line
                if str(data[line][38:39]) == "N":
                  self.atcf_array.get(f).get(line)[0] = float(data[line][35:38]) * 0.1
                else:
                  self.atcf_array.get(f).get(line)[0] = -float(data[line][35:38]) * 0.1
                if str(data[line][45:46]) == "E":
                  self.atcf_array.get(f).get(line)[1] = float(data[line][41:45]) * 0.1
                else:
                  self.atcf_array.get(f).get(line)[1] = -float(data[line][41:45]) * 0.1
                self.atcf_array.get(f).get(line)[2] = float(data[line][53:57])
                self.atcf_array.get(f).get(line)[3] = float(data[line][48:51])
                self.atcf_array.get(f).get(line)[4] = 0.0
                self.atcf_array.get(f).get(line)[5] = float(str(data[line][30:33]).strip())
                self.atcf_array.get(f).get(line)[10] = float(data[line][20:22])

                #  handle wind radii data, if it exists in the file
                if n_cols >= 64:
                    if str(data[line][63:66]) != "   ": 
                        if int(data[line][63:66]) == 34:
                            self.atcf_array.get(f).get(line)[6] = float(data[line][73:77])
                            self.atcf_array.get(f).get(line)[7] = float(data[line][79:83])
                            self.atcf_array.get(f).get(line)[8] = float(data[line][85:89])
                            self.atcf_array.get(f).get(line)[9] = float(data[line][91:95])

                datea = data[line][8:18]
                fctid = data[line][24:28]
                fhr = float(str(data[line][30:33]).strip())
                if fhr != prfhr:

                   #  append TC data to quantity-specific vectors, such as lat, lon, slp, etc.
                   ntimes = ntimes + 1
                   fhrvec.append(float(str(data[line][30:33]).strip()))
                   prfhr = float(str(data[line][30:33]).strip())

                   if str(data[line][38:39]) == "N":
                     latvec.append(float(data[line][35:38]) * 0.1)
                   elif str(data[line][38:39]) == "S":
                     latvec.append(-float(data[line][35:38]) * 0.1)

                   if str(data[line][45:46]) == "E":
                     lonvec.append(float(data[line][41:45]) * 0.1)
                   elif str(data[line][45:46]) == "W":
                     lonvec.append(-float(data[line][41:45]) * 0.1)

                   slpvec.append(float(data[line][53:57]))
                   wndvec.append(float(data[line][48:51]))

                   r34nevec.append(None)
                   r34sevec.append(None)
                   r34swvec.append(None)
                   r34nwvec.append(None)

                #  Add 34 knot wind radii data if available
                if n_cols >= 64:
                   if int(data[line][63:66]) == 34:
                      r34nevec[ntimes] = float(data[line][73:77])
                      r34sevec[ntimes] = float(data[line][79:83])
                      r34swvec[ntimes] = float(data[line][85:89])
                      r34nwvec[ntimes] = float(data[line][91:95])

            #  Append final vectors into the dictionary 
            self.atcf_array.get(f).update({'initialization_time': datea})
            self.atcf_array.get(f).update({'forecast_id': fctid})
            self.atcf_array.get(f).update({'num_lines': n_lines})
            self.atcf_array.get(f).update({'number_times': ntimes})
            self.atcf_array.get(f).update({'forecast_hour': fhrvec})
            self.atcf_array.get(f).update({'latitude': latvec}) 
            self.atcf_array.get(f).update({'longitude': lonvec})
            self.atcf_array.get(f).update({'sea_level_pressure': slpvec})
            self.atcf_array.get(f).update({'max_wind_speed': wndvec})

    def ens_lat_lon_time(self, fhr):
        '''
        Function that returns all ensemble member's latitude and longitude for a given 
        forecast hour.  The result is two vectors, one with the latitude and one with
        the ensemble TC longitude.
      
        Attributes:
            fhr (int):  forecast hour
        '''

        ens_lat = list(np.ones(len(self.atcf_files)) * self.missing)
        ens_lon = list(np.ones(len(self.atcf_files)) * self.missing)
        for n in range(len(self.atcf_files)):
           if fhr in list(self.atcf_array.get(n)['forecast_hour']):
              i = list(self.atcf_array.get(n)['forecast_hour']).index(fhr)
              ens_lat[n]=self.atcf_array.get(n)['latitude'][i]
              ens_lon[n]=self.atcf_array.get(n)['longitude'][i]

        return ens_lat, ens_lon

    def ens_intensity_time(self, fhr):
        '''
        Function that returns all ensemble member's minimum sea-level pressure and
        maximum wind speed for a given forecast hour.  The result is two vector arrays, 
        one with the ensemble SLP, and the other with the ensemble maximum wind.
      
        Attributes:
            fhr (int):  forecast hour
        '''

        ens_slp = list(np.ones(len(self.atcf_files)) * self.missing)
        ens_wnd = list(np.ones(len(self.atcf_files)) * self.missing)
        for n in range(len(self.atcf_files)):
           if fhr in list(self.atcf_array.get(n)['forecast_hour']):
              i = list(self.atcf_array.get(n)['forecast_hour']).index(fhr)
              ens_slp[n]=self.atcf_array.get(n)['sea_level_pressure'][i]
              ens_wnd[n]=self.atcf_array.get(n)['max_wind_speed'][i]

        return ens_slp, ens_wnd

    def get_best_data(self, bestfile):
        '''
        Function that reads the best track data from the specified best track file 
        and saves this within the current ATCF class.  This data will be accessed 
        through the current ATCF class. 
        ****** ROUTINE CURRENTLY UNDER CONSTRUCTION ******
      
        Attributes:
            bestfile (string):  best track file to read
        '''

#        bestpath = 'b' + storm + '.dat'
#        bestfile = os.path.join(self.src_path, bestpath)
        try:
           if os.path.isfile(bestfile):
              with open(bestfile, 'r') as fo:
                 best_data = np.array(fo.readlines())
                 return best_data
        except:
           print("{0} not found".format(bestfile))


if __name__ == '__main__':

    src1 = "/Users/parthpatwari/RA_Atmospheric_Science/Old_Code/atcf_data"
    grib_src = "/Users/parthpatwari/RA_Atmospheric_Science/GRIB_files"
    dest1 = "/Users/parthpatwari/RA_Atmospheric_Science/New_Code/atcf_data"
    atcf_src = "/Users/parthpatwari/RA_Atmospheric_Science/Old_Code/atcf_data"
    # c1 = CopyFiles(src, dest)
    # if c1.checkandcreatedir():
    #     c1.copy_filestowork()
    g1 = ReadGribFiles(grib_src, '2019082900', 180)
    a1 = Readatcfdata(atcf_src)
