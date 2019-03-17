# coding:utf-8

'''
 @author Yanqun Pan
 2018.02.10
'''

import numpy as np
import glob,os
import h5py
import matplotlib.pyplot as plt
import tables
from tables import *
from tables import Float32Atom
from tables import Filters
import pandas as pd
import gdal
import cv2


class Cbin():
    _fillvalue = -32767.0
    def __init__(self,l2dir=None,validrange=[0,10]):
        self.validRange = validrange
        self.l2dir = l2dir
        self.MONTH = [str(m).zfill(2) for m in range(1,13)]
        print(self.MONTH)
        self.SEASON = [str(m).zfill(2) for m in range(1,5)]
        self.LINES,self.PIXELS = 3801,2101
        # self.LINES,self.PIXELS = 353,384

        filter = 'COMS*GA*l2'
        self.l2files = glob.glob(os.path.join(l2dir,'2013',filter))
        data = gdal.Open(self.l2files[0])
        data_subds  = data.GetSubDatasets()

        self.LINES, self.PIXELS = gdal.Open(data_subds[0][0]).ReadAsArray().shape

        l2flags = gdal.Open(data_subds[30][0]).ReadAsArray()
        l2flags = l2flags.flatten()
        def func(item):
            if bin(item)[-2] == '1':
                return True
            else:
                return False
        # self.landmask = np.array(list(map(func, l2flags)),dtype=np.bool).reshape(3801,2101)[2057:2410,795:1179]
        self.landmask = np.array(list(map(func, l2flags)), dtype=np.bool).reshape(self.LINES, self.PIXELS)

        # self.landmask = (l2flags==2)
        # self.longitude = f['navigation_data/longitude'].value[2057:2410,795:1179]
        # self.latitude = f['navigation_data/latitude'].value[2057:2410,795:1179]
        self.longitude = gdal.Open(data_subds[0][0]).ReadAsArray()
        self.latitude =gdal.Open(data_subds[1][0]).ReadAsArray()


    def getMonthAve(self,year=None,productName="aph_443"):
        if year==None:
            year = "*"
        # files_stat = []
        for mon in self.MONTH:
            # l2files = glob.glob(os.path.join(self.l2dir,str(year),'COMS*%s%s*CONS.l2'%(year,mon)))
            # l2files = glob.glob(os.path.join(self.l2dir, str(year), 'COMS*%s%s*he5.l2' % (year, mon)))
            # l2files = glob.glob(os.path.join(self.l2dir, str(year), 'COMS*%s%s*adg.l2' % (year, mon)))
            l2files = glob.glob(os.path.join(self.l2dir, str(year),'COMS*%s%s*sert.l2' % (year, mon)))
            print (mon,len(l2files))
            # files_stat.append([year,mon,len(l2files)])
            if len(l2files) == 0:
                print('no data in %s %s!'%(year,mon))
                continue
            values = self.process(l2files,productName,self.validRange)
            values = cv2.blur(values,(5,5))
            if year=="*":
                l2binfile = os.path.join(self.l2dir,'COMS%s%s_%s_bin.l2'%("2012-2016",mon,productName.split('/')[-1]))
            else:
                l2binfile = os.path.join(self.l2dir,'%s/COMS%s%s_%s_bin.l2'%(year,year,mon,productName.split('/')[-1]))

            h5file_l2 = tables.open_file(l2binfile,'w')
            atom = Float32Atom()
            filters = Filters(complevel=5, complib='zlib')
            h5file_l2.root._v_attrs.title = 'L2 bin product(%s) produced by SKLEC'%(productName.split('/')[-1])
            h5file_l2.root._v_attrs.Scans =self.LINES
            h5file_l2.root._v_attrs.Pixels = self.PIXELS
            h5file_l2.root._v_attrs.AlgorithmName = 'ATC_MPL'
            h5file_l2.root._v_attrs.AlgorithmAuthor = 'Yanqun Pan, State Key Laboratory of Estuarine and Coastal Research'
            ca = h5file_l2.create_carray(h5file_l2.root, productName.split('/')[-1], atom, (self.LINES,self.PIXELS),filters=filters)
            ca._v_attrs._FillValue = self._fillvalue
            ca[:,:] = values

            ca = h5file_l2.create_carray(h5file_l2.root, 'longitude', atom, (self.LINES,self.PIXELS),filters=filters)
            ca._v_attrs._FillValue = self._fillvalue
            ca[:,:] = self.longitude

            ca = h5file_l2.create_carray(h5file_l2.root, 'latitude', atom, (self.LINES,self.PIXELS),filters=filters)
            ca._v_attrs._FillValue = self._fillvalue
            ca[:,:] = self.latitude

            h5file_l2.close()


            print(values.shape)
            # plt.imshow(values)
            # plt.show()
        # pd.DataFrame(data=files_stat,columns=['year','mon','files']).to_csv('goci_files_stat_%s.csv'%(year))
        return

    def getMultiYearMonAve(self):
        '''
        计算多年月平均
        :return:
        '''
        return

    def getSeasonAve(self,year=None,productName='aph_443'):
        seasonDic = {'Spring':[3,4,5],'Summer':[6,7,8],'Autumn':[9,10,11],'Winter':[12,1,2]}
        if year==None:
            year = "*"
        for season in seasonDic.keys():
            files = []
            for mon in seasonDic[season]:
                mon = str(mon).zfill(2)

                files_m = glob.glob(os.path.join(self.l2dir,str(year),'COMS*%s%s_%s_bin.l2'%(year,mon,productName)))
                if len(files) == 0:
                    print('no data in %s %s!'%(year,mon))
                    print(os.path.join(self.l2dir,'COMS*%s%s_%s_bin.l2'%(year,mon,productName)))
                    continue
                files.extend(files_m)

            values = self.process(files,productName,self.validRange)
            if year == "*":
                l2binfile = os.path.join(self.l2dir,'COMS%s%s_%s_season_bin.l2'%("2012-2016",season,productName))
            else:
                l2binfile = os.path.join(self.l2dir,'%s/COMS%s%s_%s_season_bin.l2'%(year,year,season,productName))

            h5file_l2 = tables.open_file(l2binfile,'w')
            atom = Float32Atom()
            filters = Filters(complevel=5, complib='zlib')
            h5file_l2.root._v_attrs.title = 'L2 bin product(%s) produced by SKLEC'%(productName)
            h5file_l2.root._v_attrs.Scans =self.LINES
            h5file_l2.root._v_attrs.Pixels = self.PIXELS
            h5file_l2.root._v_attrs.AlgorithmName = 'ATC_MPL'
            h5file_l2.root._v_attrs.AlgorithmAuthor = 'Yanqun Pan, State Key Laboratory of Estuarine and Coastal Research'
            ca = h5file_l2.create_carray(h5file_l2.root, productName.split('/')[-1], atom, (self.LINES,self.PIXELS),filters=filters)
            ca._v_attrs._FillValue = self._fillvalue
            ca[:,:] = values

            ca = h5file_l2.create_carray(h5file_l2.root, 'longitude', atom, (self.LINES,self.PIXELS),filters=filters)
            ca._v_attrs._FillValue = self._fillvalue
            ca[:,:] = self.longitude

            ca = h5file_l2.create_carray(h5file_l2.root, 'latitude', atom, (self.LINES,self.PIXELS),filters=filters)
            ca._v_attrs._FillValue = self._fillvalue
            ca[:,:] = self.latitude

            h5file_l2.close()

    def getYearAve(self,year=None,productName='aph_443'):
        seasonDic = {'Spring':[3,4,5],'Summer':[6,7,8],'Autumn':[9,10,11],'Winter':[12,1,2]}
        if year==None:
            year = "*"

        # files = glob.glob(os.path.join(self.l2dir,year,'COMS*%s*%s_season_bin.l2'%(year,productName)))
        files = glob.glob(os.path.join(self.l2dir, year, 'COMS*%s*%s_season_bin.l2' % (year, productName)))
        if len(files) == 0:
            print('no data in %s!'%(year))
            print(os.path.join(self.l2dir,'COMS*%s*%s_season_bin.l2'%(year,productName)))
            return

        values = self.process(files,productName,self.validRange)
        if year=="*":
            l2binfile = os.path.join(self.l2dir,'COMS%s%s_year_bin.l2'%("2012-2016",productName))
        else:
            l2binfile = os.path.join(self.l2dir,'%s/COMS%s%s_year_bin.l2'%(year,year,productName))
        h5file_l2 = tables.open_file(l2binfile,'w')
        atom = Float32Atom()
        filters = Filters(complevel=5, complib='zlib')
        h5file_l2.root._v_attrs.title = 'L2 bin product(%s) produced by SKLEC'%(productName)
        h5file_l2.root._v_attrs.Scans =self.LINES
        h5file_l2.root._v_attrs.Pixels = self.PIXELS
        h5file_l2.root._v_attrs.AlgorithmName = 'ATC_MPL'
        h5file_l2.root._v_attrs.AlgorithmAuthor = 'Yanqun Pan, State Key Laboratory of Estuarine and Coastal Research'
        ca = h5file_l2.create_carray(h5file_l2.root, productName.split('/')[-1], atom, (self.LINES,self.PIXELS),filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:,:] = values

        ca = h5file_l2.create_carray(h5file_l2.root, 'longitude', atom, (self.LINES,self.PIXELS),filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:,:] = self.longitude

        ca = h5file_l2.create_carray(h5file_l2.root, 'latitude', atom, (self.LINES,self.PIXELS),filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:,:] = self.latitude

        h5file_l2.close()
        return 0



    def process(self,files,productName,validRange):
        index = np.zeros((self.LINES,self.PIXELS))
        values = np.zeros((self.LINES,self.PIXELS))
        for fl in files:
            f = h5py.File(fl)
            try:
                value = f[productName].value
            except:
                print(fl)
            value[(value>validRange[1]) | (value<validRange[0])] = 0
            vc = value.copy()
            vc[~((vc>validRange[1]) | (vc<validRange[0]))] = 1
            index += vc
            values += value
            f.close()

        values[index>0] = values[index>0]/index[index>0]
        values[self.landmask] = self._fillvalue
        return values


    def getAve(self,year,filter,productName):
        l2files = glob.glob(os.path.join(self.l2dir, str(year),filter))
        # files_stat.append([year,mon,len(l2files)])
        if len(l2files) == 0:
            print('no data in %s!'%(year))
            return
        print(len(l2files))
        values = self.process(l2files,productName,[0,4000])
        values = cv2.blur(values,(5,5))
        if year=="*":
            l2binfile = os.path.join(self.l2dir,'COMS%s%s_bin.l2'%("2012-2016",productName.split('/')[-1]))
        else:
            l2binfile = os.path.join(self.l2dir,'%s/COMS%s%s_bin.l2'%(year,year,productName.split('/')[-1]))

        h5file_l2 = tables.open_file(l2binfile,'w')
        atom = Float32Atom()
        filters = Filters(complevel=5, complib='zlib')
        h5file_l2.root._v_attrs.title = 'L2 bin product(%s) produced by SKLEC'%(productName.split('/')[-1])
        h5file_l2.root._v_attrs.Scans =self.LINES
        h5file_l2.root._v_attrs.Pixels = self.PIXELS
        h5file_l2.root._v_attrs.AlgorithmName = 'ATC_MPL'
        h5file_l2.root._v_attrs.AlgorithmAuthor = 'Yanqun Pan, State Key Laboratory of Estuarine and Coastal Research'
        ca = h5file_l2.create_carray(h5file_l2.root, productName.split('/')[-1], atom, (self.LINES,self.PIXELS),filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:,:] = values

        ca = h5file_l2.create_carray(h5file_l2.root, 'longitude', atom, (self.LINES,self.PIXELS),filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:,:] = self.longitude

        ca = h5file_l2.create_carray(h5file_l2.root, 'latitude', atom, (self.LINES,self.PIXELS),filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:,:] = self.latitude

        h5file_l2.close()


        print(values.shape)

if __name__ == '__main__':

    # 有效值范围，只有在此范围内的才bin,否则mask
    validRange = [0,10]

    # cbin = Cbin(l2dir='E:\\Data\\GOCI')
    # for year in ['2012','2013','2014','2015','2016']:
    #     print(year)
    #     cbin.getMonthAve(year=year,productName='adg/adg_443')
    #     cbin.getSeasonAve(year=year,productName='adg_443')
    #     cbin.getYearAve(year=year,productName='adg_443')
    # cbin.getMonthAve(productName='adg/adg_443')
    # cbin.getSeasonAve(productName='adg_443')
    # cbin.getYearAve(productName='adg_443')


    cbin = Cbin(l2dir='E:\\GOCI L2')
    # cbin.getMonthAve(year=2016,productName='POC/POC')
    # cbin.getSeasonAve(year=2016, productName='POC')
    #
    # cbin.getMonthAve(year=2016, productName='Chla/Chla-YOC')
    # cbin.getSeasonAve(year=2016, productName='Chla-YOC')

    cbin.getMonthAve(year=2014, productName='SSC/SSC_SERT')
    cbin.getSeasonAve(year=2014, productName='SSC_SERT')

    # cbin.getAve(year='2015',filter='COMS*POC.l2',productName='Chla/Chla-YOC')

    # cbin.getAve(year='2012', filter='COMS*POC.l2', productName='POC/POC')
    # cbin.getAve(year='2012', filter='COMS*ssc*.l2', productName='SSC/SSC_SERT')

    # cbin.getMonthAve(productName='POC')
    # cbin.getSeasonAve(productName='Chla-YOC')
    # cbin.getSeasonAve(productName='SSC_SERT')


