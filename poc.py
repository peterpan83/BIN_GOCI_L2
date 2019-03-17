'''
 @author Yanqun Pan
 2018.02.10
'''

import gdal
import tables
from tables import *
from tables import Float32Atom
from tables import Filters
import glob,os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt

class Cprocess():

    def __init__(self,dir,filter):
        self.dir = dir
        self.filter = filter
        self.l2files = glob.glob(os.path.join(dir,filter))
        self.waves = [412,443,490,555,660,680,745,865]
        self.dic = {}
        data = gdal.Open(self.l2files[0])
        data_subds  = data.GetSubDatasets()
        for index,item in enumerate(data_subds):
            for wa in self.waves:
                if item[1].find('Rrs_%d'%(wa))>0:
                    self.dic['Rrs_%d'%(wa)] = index
                if item[1].find('ssc')>0:
                    self.dic['ssc'] = index
        self.lines,self.pixels = gdal.Open(data_subds[0][0]).ReadAsArray().shape
        self._fillvalue = -32767.0
        self.slope,self.intercept =  2.0E-6,0.05
        print(self.dic)

    def run(self,product='ssc'):
        for l2f in self.l2files:
            self.processSingle(l2f,product)


    def processSingle(self,l2f,product='ssc'):
        print(l2f)
        if product == 'ssc':
            self.__ssc(l2f)
        else:
            self.__poc(l2f)


    def __ssc(self,l2f):
        data = gdal.Open(l2f)
        ds  = data.GetSubDatasets()
        ssc_sert = gdal.Open(ds[self.dic['ssc']][0]).ReadAsArray()*1000
        sscfname = '%s_ssc_sert.l2' % (l2f.split('.')[0])
        if os.path.exists(sscfname):
            return
        h5file_l2 = tables.open_file(sscfname, 'w')
        shape = (self.lines, self.pixels)
        atom = Float32Atom()
        filters = Filters(complevel=5, complib='zlib')

        h5file_l2.root._v_attrs.title = 'GOCI SSC product produced by SKLEC,Yanqun Pan'

        grpChla = h5file_l2.create_group(h5file_l2.root, 'SSC', 'SSC')
        grpChla._v_attrs.Scans = self.lines
        grpChla._v_attrs.Pixels = self.pixels
        grpChla._v_attrs.AlgorithmName = 'ATC_MPL'

        ca = h5file_l2.create_carray(grpChla, 'SSC_SERT', atom, shape, filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:] = ssc_sert
        h5file_l2.close()

    def __poc(self,l2f):
        data = gdal.Open(l2f)
        ds  = data.GetSubDatasets()
        Rrs412 = gdal.Open(ds[self.dic['Rrs_412']][0]).ReadAsArray()*self.slope+self.intercept
        Rrs443 = gdal.Open(ds[self.dic['Rrs_443']][0]).ReadAsArray()*self.slope+self.intercept
        Rrs490 = gdal.Open(ds[self.dic['Rrs_490']][0]).ReadAsArray()*self.slope+self.intercept
        Rrs555 = gdal.Open(ds[self.dic['Rrs_555']][0]).ReadAsArray()*self.slope+self.intercept
        Rrs660 = gdal.Open(ds[self.dic['Rrs_660']][0]).ReadAsArray() * self.slope + self.intercept
        Rrs680 = gdal.Open(ds[self.dic['Rrs_680']][0]).ReadAsArray() * self.slope + self.intercept

        ssc_sert = gdal.Open(ds[self.dic['ssc']][0]).ReadAsArray()
        chla_oc3 = gdal.Open(ds[self.dic['chl_oc3']][0]).ReadAsArray()

        m412,m443,m490,m555,m660,m680 = Rrs412<0, Rrs443<0, Rrs490<0, Rrs555<0,Rrs660<0,Rrs680<0

        mask = m412|m443|m490|m555

        Rrs412_m = ma.array(Rrs412, mask=mask)
        Rrs443_m = ma.array(Rrs443, mask=mask)
        Rrs490_m = ma.array(Rrs490, mask=mask)
        Rrs555_m = ma.array(Rrs555, mask=mask)
        Rrs660_m = ma.array(Rrs660, mask=mask)
        Rrs680_m = ma.array(Rrs680, mask=mask)

        tempR = (Rrs443_m/Rrs555_m)*np.power(Rrs412_m/Rrs490_m,-1.012)
        temp  = 0.342-2.511*np.log10(tempR)-0.277*np.power(np.log10(tempR),2)
        chla = np.power(10,temp)
        print(chla.data)
        chla_data = chla.data
        chla_data[mask] = self._fillvalue
        # plt.imshow(chla_data)
        # plt.show()


        pocfname = '%s_POC.l2'%(l2f.split('.')[0])

        if os.path.exists(pocfname):
            return
        h5file_l2 = tables.open_file(pocfname, 'w')
        shape = (self.lines, self.pixels)
        atom = Float32Atom()
        filters = Filters(complevel=5, complib='zlib')

        h5file_l2.root._v_attrs.title = 'GOCI POC product produced by SKLEC,Yanqun Pan'

        grpChla = h5file_l2.create_group(h5file_l2.root, 'Chla', 'remote sensing reflectance')
        grpChla._v_attrs.Scans = self.lines
        grpChla._v_attrs.Pixels = self.pixels
        grpChla._v_attrs.AlgorithmName = 'ATC_MPL'

        ca = h5file_l2.create_carray(grpChla, 'Chla-OC3', atom, shape, filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:] = chla_oc3

        ca = h5file_l2.create_carray(grpChla, 'Chla-YOC', atom, shape, filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:] = chla_data

        mtemp1 =  Rrs660 > Rrs490
        mtemp2  = Rrs660 > Rrs680

        ratio = Rrs490_m / Rrs555_m
        POC = np.zeros((self.lines,self.pixels),dtype=np.float)
        POC[mtemp1] = ssc_sert[mtemp1]*1000*5.06+37.33
        POC[(~mtemp1) & (mtemp2)] = 87.3*np.power(ratio.data[(~mtemp1) & (mtemp2)],-2.04)
        POC[(~mtemp1) & (~mtemp2)] =  69.9*np.power(chla_data[(~mtemp1) & (~mtemp2)],0.63)

        POC[mask] = self._fillvalue

        h5file_l2.root._v_attrs.title = 'GOCI POC product produced by SKLEC,Yanqun Pan'

        grpPOC = h5file_l2.create_group(h5file_l2.root, 'POC', 'remote sensing reflectance')
        grpPOC._v_attrs.Scans = self.lines
        grpPOC._v_attrs.Pixels = self.pixels
        grpPOC._v_attrs.AlgorithmName = 'ATC_MPL'

        ca = h5file_l2.create_carray(grpPOC, 'POC', atom, shape, filters=filters)
        ca._v_attrs._FillValue = self._fillvalue
        ca[:] = POC

        h5file_l2.close()


if __name__ == '__main__':
    cp = Cprocess('E:\\GOCI L2\\2014',filter='COMS_GOCI_L1B_GA_*he5.GA.l2')
    # cp.run('POC')
    cp.run('ssc')



