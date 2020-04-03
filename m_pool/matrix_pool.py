#!/usr/bin/env python
# -*- coding: utf8 -*-

r"""
A wrapper for numpy arrays providing named axes, interpolation, iteration, disk persistence and numerical calcs


M_Pool wraps multidimensional numpy arrays to provide the following features::

    #. MatrixPool objects contain related Axis and Matrix objects
        - MP = MatrixPool(name='N2O4_MMH')
        
    #. Axis objects are added by name and interpolation transform (used to linearize interpolation)
        - epsAxis = Axis({'name':'eps', 'valueL':[10., 20., 30., 40., 50.], 'units':'', 'transform':'log10'})
        - pcAxis = Axis({'name':'pc', 'valueL':[100.,200.,300,400], 'units':'psia', 'transform':'log10'})
        - mrAxis = Axis({'name':'mr', 'valueL':[1,2,3], 'units':'', 'transform':''})
        
    #. Matrix objects added by name 
        - M = MP.add_matrix( name='cea_isp', units='sec', axisNameL=['eps','pc','mr'] )
        
    #. Find interpolated minimum or maximum
        - interpD, max_val = M.solve_interp_max( order=3, method='TNC', tol=1.0E-8)
            - where interpD and max_val look something like:
            - interpD = {'pc': 225.00641803120988, 'eps': 34.991495018803455, 'mr': 1.7499612975876655}
            - max_val = -0.000155216246295
    
    #. Disk-based persistence
        - Save to pickle or hdf5 file
            - MP.save_to_pickle() # saves MP to "N2O4_MMH_matrix.pool"
            
    #. Built-in statistics (standard deviation, median, mean/average, sum, minimum, maximum
        - M.get_range()
        - M.get_ave()
        - M.get_mean()
        - M.get_std()
        - M.get_median()
    
    #. Interpolation on axes with named values
        - interp_val = M.interp(order=2, pc=100, eps=20, mr=2.0)
        - Uses transformed axes to help linearize interpolation
        
    #. Iterate over matrix
        - for indeces,D,val in M.full_iter_items():
            - gives something like:
            - (0, 0, 0) {'pc': 100.0, 'eps': 10.0, 'mr': 1.0} 111.0
            - (0, 0, 1) {'pc': 100.0, 'eps': 10.0, 'mr': 2.0} 112.0
            - (0, 0, 2) {'pc': 100.0, 'eps': 10.0, 'mr': 3.0} 113.0
            - ...



M_Pool
Copyright (C) 2015  Charlie Taylor

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

-----------------------

"""
import os
here = os.path.abspath(os.path.dirname(__file__))


# for multi-file projects see LICENSE file for authorship info
# for single file projects, insert following information
__author__ = 'Charlie Taylor'
__copyright__ = 'Copyright (c) 2015 Charlie Taylor'
__license__ = 'GPL-3'
exec( open(os.path.join( here,'_version.py' )).read() )  # creates local __version__ variable
__email__ = "cet@appliedpython.com"
__status__ = "3 - Alpha" # "3 - Alpha", "4 - Beta", "5 - Production/Stable"

import pickle as pickle
import os
import numpy as np
try:
    import tables
except:
    print('... WARNING... No HDF5 Support Available (failed import tables)')
from m_pool.axis_obj import Axis
from m_pool.axis_pool import AxisPool
from m_pool.matrix_obj import Matrix

def read_mpool( name='' ):
    if not name:
        name = '%s_matrix.pool'%(name)
    if not name.endswith('_matrix.pool'):
        name += '_matrix.pool'
        
            
    if os.path.isfile(name):            
        MP = MatrixPool(name=name)
        MP.read_from_pickle(fname=name) # Opens "CHECK_matrix.pool"
        return MP
    else:
        print('ERROR... No Matrix File Found For:', name)

class MatrixPool(object):
    """
    A wrapper for numpy arrays providing named axes, interpolation, iteration, disk persistence and numerical calcs

    A MatrixPool object is a Collection of Axis and Matrix objects. 
    Like AxisPool, it is used to define a collection of Matrix objects
         
    *** Structured to easily pickle via a dictionary of named values for properties. ***
    """
    
    def __init__(self, name='Pool Data', descD=None):
        
        self.name = name
        self.descD = descD  # descriptive data in dictionary format
        
        self.axisPoolObj = AxisPool()
        
        self.matrixL = [] # list of Matrix objects
        self.matrixD = {} # cross ref by name
        
    def add_axis(self, axOrD):
        self.axisPoolObj.add_axis( axOrD )
    
    def rename_matrix(self, old_name, new_name ):
        M = self.matrixD[old_name]
        M.name = new_name
        del self.matrixD[old_name]
        self.matrixD[new_name] = M
    
    def add_matrix(self, name='matrixName', units='', axisNameL=None,
        matValArr=None):
        
        D={'name':name,  'units':units, 'matValArr':matValArr,
            'axisNameL':axisNameL, 'axisPoolObj':self.axisPoolObj}
        self.matrixL.append( Matrix(D) )
        
        self.matrixD[name] = self.matrixL[-1]
    
        return self.matrixL[-1]

    def get_matrix_by_name(self, matName):
        return self.matrixD.get(matName,None)

    def get_axis_by_name(self, axName):
        return self.axisPoolObj[axName]

    def __str__(self):
        sL = ['MatrixPool: %s'%self.name ]
        for M in self.matrixL:
            sL.append( '  Matrix:%s, shape=%s, units=%s, %%full=%i'%(M.name, M.shape(), M.units, M.iPercentFull()) )
        return '\n'.join( sL )
    
    def summ(self, use_long=False):
        sL = ['MatrixPool: %s'%self.name ]
        if self.descD is not None:
            keyL = sorted( list(self.descD.keys()), key=str.lower )
            try:
                m = max( [len(str(k)) for k in keyL] )
            except:
                m = 6
            fmt = '%' + '%is'%m
            for key in keyL:
                sL.append( '  ' + fmt%str(key) + ' : ' + str(self.descD[key]) )
        
        for M in self.matrixL:
            if use_long:
                s = M.long_summ()
            else:
                s = M.short_summ()
            ssL = s.split('\n')
            for s in ssL:
                sL.append( '    ' + s )
        return '\n'.join( sL )
        
    
    def __len__(self):
        return len( self.matrixL )
    
    def __iter__(self):
        for value in self.matrixL:
            yield value
    
    def __getitem__(self, i): # retrieve as: A[i]
        return self.matrixL[i]
        
    def save_to_hdf5(self, fname=None):
        if fname==None:
            fname = '%s_mpool.h5'%(self.name)
            
        h5file = tables.open_file(fname, mode='w')
        # Get the HDF5 root group
        root = h5file.root
        
        h5file.create_array(root, 'axes_name_list', [A.name for A in self.axisPoolObj], "String array")
        h5file.create_array(root, 'matrix_name_list', [M.name for M in self.matrixL], "String array")
        axes_group = h5file.create_group(root, 'axes', 'All the Axes used in Matrix Pool')
        mat_group = h5file.create_group(root, 'matrices', 'All the Matrices used in Matrix Pool')
        
        '''  {'name':self.name, 'valueL':self.valueL, 'units':self.units, 
            'transform':self.transform, 'roundDigits':self.roundDigits}'''
        for A in self.axisPoolObj:
            d = h5file.create_array(axes_group, A.name, A.valueArr) 
            h5file.create_array(axes_group, 'transform_%s'%A.name, A.transArr) 
            h5file.create_array(axes_group, 'desc_%s'%A.name, 
                [str(A.units), str(A.transform), str(A.roundDigits)], "String array")
            #d.transform_desc = A.transform
            #d.roundDigits_desc = A.roundDigits
            #d.units_desc = A.units
        
        # Just use last axis "A" for atom creation
        atom = tables.Atom.from_dtype(A.valueArr.dtype)
        filters = tables.Filters(complib='blosc', complevel=5)    
        
        for M in self.matrixL:    
            ds = h5file.create_carray(mat_group, M.name, atom, M.matValArr.shape, filters=filters)
            ds.attrs.units = M.units
            ds[:] = M.matValArr
                
            h5file.create_array(mat_group, M.name+'_axis_list', [aname for aname in M.axisNameL], "String array")
            
        h5file.flush()
        h5file.close()
    
    def read_from_hdf5(self, fname=None):
        if fname==None:
            fname = '%s_mpool.h5'%(self.name)
            print('Reading:', fname)
            
        if os.path.exists(fname):            
            h5file = tables.open_file(fname, mode='r')
               
            # reinit self
            self.axisPoolObj = AxisPool()
            self.matrixL = [] # list of Matrix objects
            self.matrixD = {} # cross ref by name
            self.descD = None
            
            root = h5file.root
            
            # First get the axes 
            axis_nameL = [_.decode('utf8') for _ in root.axes_name_list]
            print('axis_nameL =',axis_nameL)

            for aname in axis_nameL:
                #print('Getting axis_nameL, aname =', aname)
                val_arr = getattr( root.axes, aname ).read()
                dname = 'desc_' + aname
                desc = getattr( root.axes, dname ).read()
                units, trans, alen = desc
                #print('units, trans, alen =',units, trans, alen)
                units, trans, alen = units.decode('utf8'), trans.decode('utf8'), alen.decode('utf8')
                
                A = Axis({'name':aname, 'valueL':val_arr, 'units':units, 'transform':trans})
                self.add_axis( A )
            
            # Then get the matrices
            matrix_nameL = [_.decode('utf8') for _ in root.matrix_name_list]
            print('matrix_nameL =',matrix_nameL)

            for mname in matrix_nameL:
                m = getattr( root.matrices, mname ).read()
                units = getattr( root.matrices, mname ).attrs.units
                a_list = getattr( root.matrices, mname+'_axis_list' ).read()
                #print('units =', units, type(units))
                
                a_list = [_.decode('utf8') for _ in a_list]
                #print('a_list =', a_list)
                self.add_matrix( name=mname, units=units, axisNameL=list(a_list), matValArr=m )
            
            h5file.close()

    def save_to_pickle(self, fname=None):
        if fname==None:
            fname = '%s_matrix.pool'%(self.name)
        if not fname.endswith('_matrix.pool'):
            fname += '_matrix.pool'
        
        
        D = {}
        D['axisL'] = [A.get_pickleable_dict() for A in self.axisPoolObj]
        D['matrixL'] = [M.get_pickleable_dict() for M in self.matrixL]
        D['descD'] = self.descD
        
        fOut = open(fname, 'wb')
        pickle.dump( D, fOut )
        fOut.close()
        
    def read_from_pickle(self, fname=None):
        if fname==None:
            fname = '%s_matrix.pool'%(self.name)
        if not fname.endswith('_matrix.pool'):
            fname += '_matrix.pool'
            
        print('Reading:', fname)
            
        if os.path.exists(fname):            
            
            try:
                fInp = open(fname, 'rb')
                D = pickle.load( fInp )
                got_it = True
            except:
                got_it = False
            finally:
                fInp.close()
                
            if not got_it:
                fInp = open(fname, 'rb')
                #D = pickle.load( fInp, encoding='bytes' )
                D = pickle.load( fInp, encoding='latin1' )
                fInp.close()
            
            # reinit self
            self.axisPoolObj = AxisPool()
            self.matrixL = [] # list of Matrix objects
            self.matrixD = {} # cross ref by name
            
            self.descD = D.get('descD',{}) # if descriptive dictionary was saved, restore it
            
            for AD in D['axisL']:
                #print 'Adding',AD
                self.add_axis( AD )
            #print self.axisPoolObj
            #print
            #print self
                
            for MD in D['matrixL']:
                #print 'add Matrix',MD
                self.add_matrix( **MD )
        else:
            print('...WARNING... could not find:',fname)

if __name__=="__main__":

    MP = MatrixPool(name='N2O4_MMH')
    epsAxis = Axis({'name':'eps', 'valueL':[10., 20., 30., 40., 50.], 'units':'', 'transform':'log10'})
    pcAxis = Axis({'name':'pc', 'valueL':[100.,200.,300,400], 'units':'psia', 'transform':'log10'})
    mrAxis = Axis({'name':'mr', 'valueL':[1,2,3], 'units':'', 'transform':''})
    for A in [epsAxis, pcAxis, mrAxis]:
        MP.add_axis( A )
    
    M = MP.add_matrix( name='cea_isp', units='sec', axisNameL=['eps','pc','mr'] )
    for eps in epsAxis:
        for pc in pcAxis:
            for mr in mrAxis:
                M.setByName( pc=pc, eps=eps, mr=mr, val=eps+pc+mr )

    
    M = MP.add_matrix( name='SEA_isp', units='sec', axisNameL=['eps','pc','mr'] )
    for eps in epsAxis:
        for pc in pcAxis:
            for mr in mrAxis:
                M.setByName( pc=pc, eps=eps, mr=mr, val=eps+pc+mr+0.321 )
    

    print(MP)
    #MP.save_to_pickle()
    #MP.save_to_hdf5()
    print('='*55)
    print(MP.summ())
