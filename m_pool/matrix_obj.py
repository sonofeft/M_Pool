#!/usr/bin/env python
# -*- coding: utf8 -*-
import sys
import itertools
import copy
import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator
try:
    from scipy.optimize import minimize
except:
    print("...WARNING... scipy.optimize.minimize did NOT import...")
    print("   ... min/max functions are UNAVAILABLE ...")

from m_pool.axis_obj import Axis
from m_pool.axis_pool import AxisPool, axis_obj_dammit
from m_pool.InterpProp_scipy import InterpProp

try:
    import pylab
    got_pylab = True
except:
    got_pylab = False

class Matrix(object):
    '''An Matrix object holds data for N dimensional data
       There are N Axis objects for the data.
       
       The data is a single number indexed by the axes index values.
         
       *** Structured to easily pickle via a dictionary of named values for properties. ***
    '''

    def __init__(self, D={'name':'matrixName', 'matValArr':None, 'units':'', 
            'axisNameL':None, 'axisPoolObj':None} ):
        '''Initialize with a dictionary so that pickle files can easily save and read objects
        
           axisNameL holds the names of axes that are in the axisPoolObj.
           The Matrix is dimensioned by the size of the axes in the order specified.
           
           An Axis obj can be shared by many Matrix objects.
        '''
        
        self.name = D.get('name','UnkMatrix')
        self.matValArr = D.get('matValArr', None)
        self.units = D.get('units','')
        
        # Let it crash if axisNameL and axisPoolObj are not specified
        try:
            self.axisNameL = D.get('axisNameL')
            self.axisPoolObj = D.get('axisPoolObj')
        except:
            print('ERROR... both axisNameL and axisPoolObj MUST be specified in Matrix')
            sys.exit()
        
        self.axisL = [self.axisPoolObj.axisD[name] for name in self.axisNameL]
        self.NumAxes = len( self.axisL )
        shape = [len(A) for A in self.axisL]
            
        # Init to Zeros if axes specified, but data not specified
        if self.matValArr is None and shape:
            self.matValArr = np.zeros( shape )
            
        self.axisPoolObj.connectMatrixToAxes(self, self.axisNameL)
        
        # temporary list of numpy matrices used for interpolation
        self.terp_mL = [self.matValArr] # list of matrices used for interpolation
        
        self.terp_reg_grid = None # will be initialized if used
        self.terp_reg_grid_shape = None

    def solve_interp_min(self, order=3, method='TNC', tol=1.0E-8): # method can be: SLSQP, TNC
        return self.solve_interp_minmax( find_min=True, order=order, method=method, tol=tol)

    def solve_interp_max(self, order=3, method='TNC', tol=1.0E-8): # method can be: SLSQP, TNC
        return self.solve_interp_minmax( find_min=False, order=order, method=method, tol=tol)

    def solve_interp_minmax(self, find_min=False, 
        order=3, method='TNC', tol=1.0E-8): # method can be: SLSQP, TNC
        
        boundsL = []
        startValL = []
        axisNameL = []
        mn,mx = self.get_min_max()
        range = mx - mn
        interpD = {} # dictionary of axis values
        
        if find_min:
            iminmax = self.get_minima_indeces()
        else:
            iminmax = self.get_peak_indeces()
        
        for i,im in enumerate( iminmax ):
            #print 'minmax value at %s=%g'%(self.axisL[i].name, self.axisL[i][im])
            #EPS=1.0E-10*abs(self.axisL[i][-1] - self.axisL[i][0])
            boundsL.append( (self.axisL[i][0],self.axisL[i][-1]) )
            startValL.append( self.axisL[i][im] )
            axisNameL.append( self.axisL[i].name )
            interpD[self.axisL[i].name] = self.axisL[i][im]
        #print 'minmax value =',self.matValArr[ iminmax ],'   Min =',mn,'  Max =',mx
        #print 'boundsL =',boundsL
        #print 'startValL =',startValL
        #print 'axisNameL =',axisNameL
        #print 'interpD =',interpD
        
        def fun( row ): # row is in axis-order from self.axisL
            for i,val in enumerate(row):
                interpD[ axisNameL[i] ] = val
            mval = self.interp(order=order, **interpD )
            
            norm_val = float( (mval-mn)/range ) # normalize to help convergence
            if find_min:
                return norm_val
            else:
                return -norm_val
            
        res = minimize(fun, tuple(startValL), method=method, 
            bounds=tuple(boundsL), tol=tol, options={'disp':False})
        print(res)
        
        fun( res.x )# make sure interpD is set
        
        minmax_val = float( self.interp( **interpD ) )
        
        return interpD, minmax_val
        

    def interp(self, order=1, **kwds): # kwds contains axis names... returns interpolated val
        if order>1:
            return self.interp_higher_order( order=order, **kwds )
        else:
            return self.interp_linear( order=order, **kwds )
    
    def interp_linear(self, **kwds ):
        
        if (self.terp_reg_grid is None) or (self.terp_reg_grid_shape != self.matValArr.shape ):
            
            self.terp_reg_grid_shape = self.matValArr.shape
            
            axis_valL = [ A.get_trans_valueL() for A in self.axisL ]
            self.terp_reg_grid = RegularGridInterpolator(axis_valL, self.matValArr)
        
        ptArr = np.array( [A.transObj( kwds[ A.name ] ) for A in self.axisL] )
        
        ans = self.terp_reg_grid( ptArr )
        #print( 'ans=',ans )
        return ans[0]
        
    def interp_higher_order(self, order=3, **kwds): # kwds contains axis names... returns interpolated val
        '''
        Call as: M.interp(order=3, pc=100, eps=20, mr=2.0)
        
        Uses scipy.interpolate.interp1d
        '''
        
        
        # Only generate list of temporary matrices if 1st time, or if shape change
        if (len(self.terp_mL)==1) or (self.terp_mL[0].shape != self.matValArr.shape):
            #print 'orig shape =',self.matValArr.shape
            self.terp_mL[0] = self.matValArr # list of matrices used for interpolation
            
            #remove first dimension from each subsequent member of self.terp_mL
            next_shape = list( self.matValArr.shape )[1:]
            #print 'next_shape =',next_shape
            while len(next_shape)>0:
                self.terp_mL.append( np.zeros( next_shape ) )
                next_shape = next_shape[1:]
                #print 'next_shape =',next_shape
        else:
            self.terp_mL[0] = self.matValArr # verify 1st matrix is current
        
        #  interp from previous matrix  for next matrix
        for ia,m in enumerate(self.terp_mL[1:]): # ia is index into self.axisL for current axis
            A = self.axisL[ia]
            xval = A.transObj( kwds[ A.name ] )
            kind = min(len(A)-1, order)
            #print '... interpolating into',A.name,' xval=',xval,A
            for mindeces in itertools.product(*(list(range(s)) for s in m.shape)):
                # mindeces is a tuple index into m
                # indeces is index into last m
                
                yL = []
                #print 'mindeces =',mindeces
                mindecesL = list( mindeces )
                for iv,vax in enumerate( A ):
                    indeces = tuple( [iv] + mindecesL )
                    val = self.terp_mL[ia][indeces]
                    #print indeces, val
                    yL.append( val )
                #print 'xL=',A.transArr
                #print 'yL=',yL
                try:
                    # do NOT set , fill_value="extrapolate" 
                    #     ... let if fail so Extrapolating logic is used.
                    m[mindeces] = interp1d( A.transArr , yL, kind=kind, fill_value="extrapolate")(xval)
                except:
                    #print('Extrapolating',A.name,'axis =',A.transArr,'  xval=',xval)
                    print('Extrapolating',A.name,'axis  %s='%A.name, kwds[ A.name ])
                    #print('           yL =',yL)
                    if xval>=A.transArr[-2]:
                        m[mindeces] = yL[-1] # assume out of bounds at high end
                    else:
                        m[mindeces] = yL[0] # assume out of bounds at low end
        
        #print 'Last matrix(array) =',self.terp_mL[-1]
        A = self.axisL[-1]
        kind = min(len(A)-1, order)
        xval = A.transObj( kwds[ A.name ] )
        m = self.terp_mL[-1]
        #print 'm =',m
        #print 'axis =',A,'  xval=',xval
        try:
            result = interp1d( A.transArr, m, kind=kind, fill_value="extrapolate")( xval )
        except:
            print('Extrapolating','axis =',A,'  xval=',xval)
            print('           m =',m)
            if xval>=A.transArr[-2]:
                result = m[-1] # assume out of bounds at high end
            else:
                result = m[0] # assume out of bounds at low end
        
        #print 'type(result)=',type(result), result.shape
        #return result
        return  float( result )
        
        
    def numNonZero(self):
        return np.count_nonzero( self.matValArr )
        
    def iPercentFull(self): # return an integer percent full
        ntotal = 1
        for i in self.matValArr.shape:
            ntotal *= i
        nfull = np.count_nonzero( self.matValArr )
        return (100*nfull) / ntotal
    
    def get_pickleable_dict(self):
        '''Note that matrix_pool supplies axisPoolObj for pickled Matrix'''
        return {'name':self.name, 'matValArr':self.matValArr, 'units':self.units,
            'axisNameL':self.axisNameL}
            
    def insert_dimension(self, iaxis,i ):
        newMat = np.insert( self.matValArr, i, 0.0, axis=iaxis )
        self.matValArr = newMat

    def long_summ(self):
        sL = [self.short_summ()]
        
        sL.append( 'get_range   = %s'%( self.get_range(), ))
        sL.append( 'get_ave     = %s'%( self.get_ave(), ))
        sL.append( 'get_mean    = %s'%( self.get_mean(), ))
        sL.append( 'get_std     = %s'%( self.get_std(), ))
        sL.append( 'get_median  = %s'%( self.get_median(), ))
        sL.append( 'get_min_max = %s'%( self.get_min_max(), ))

        return '\n'.join( sL )

    def short_summ(self):
        if self.matValArr is None:
            sL = ['Matrix %s (shape=%s) %s (units=%s)'%(self.name, str(self.matValArr),self.name, self.units)]
        else:
            sL = ['Matrix %s (shape=%s) %s (units=%s)'%(self.name, str(self.matValArr.shape),self.name, self.units)]
        for A in self.axisL:
            s = str(A)
            ssL = s.split('\n')
            for s in ssL:
                sL.append( '    ' + s )
            #sL.append( str(A) )
        
        return '\n'.join( sL )

    def __str__(self):
        s = self.short_summ()
        sL = s.split('\n')
        #if self.matValArr is None:
        #    sL = ['Matrix %s (shape=%s) %s (units=%s)'%(self.name, str(self.matValArr),self.name, self.units)]
        #else:
        #    sL = ['Matrix %s (shape=%s) %s (units=%s)'%(self.name, str(self.matValArr.shape),self.name, self.units)]
        #for A in self.axisL:
        #    sL.append( str(A) )
        
        sL.append( str(self.matValArr) )
        return '\n'.join( sL )
        
    def __getitem__(self, iL):
        return self.matValArr[ tuple(iL) ]

        
    def __setitem__(self, iL, val): # use as M[(i,j,k)] = val
        if val is None:
            print('ERROR... illegal value for "val" in Matrix.set.  val =',val)
        else:
            self.matValArr[tuple(iL)] = float(val)

    def setByName(self, **kwds): # kwds contains axis names and "val"
        '''Call as: M.setByName(pc=100, eps=20, mr=2.0, val=29.23)'''
        
        iL = [] # list of indeces into matrix array
        for A in self.axisL:
            iL.append( A.getExactIndex( kwds[A.name] ) )
        self.matValArr[tuple(iL)] = float( kwds['val'] )

    def getByName(self, **kwds): # kwds contains axis names... returns val
        '''Call as: M.getByName(pc=100, eps=20, mr=2.0)'''
        
        iL = [] # list of indeces into matrix array
        for A in self.axisL:
            iL.append( A.getExactIndex( kwds[A.name] ) )
        return self.matValArr[tuple(iL)]
    
    def get_list_of_peak_indeces(self):
        """Returns a list of all occurances of max value"""
        
        max_val = self.get_max()
        return np.argwhere( self.matValArr == max_val )
        
    
    def get_peak_indeces(self):
        """Returns 1st occurance of max value"""
        imax = np.unravel_index(self.matValArr.argmax(), self.matValArr.shape)
        return imax

    def get_peak_dict(self):
        """Returns 1st occurance of max value"""
        imax = np.unravel_index(self.matValArr.argmax(), self.matValArr.shape)
        D = {}
        for i,im in enumerate( imax ):
            D[self.axisL[i].name] = self.axisL[i][im]
        return D
    
    def get_minima_indeces(self):
        """Returns 1st occurance of min value"""
        imin = np.unravel_index(self.matValArr.argmin(), self.matValArr.shape)
        return imin

    def get_minima_dict(self):
        """Returns 1st occurance of min value"""
        imin = np.unravel_index(self.matValArr.argmin(), self.matValArr.shape)
        D = {}
        for i,im in enumerate( imin ):
            D[self.axisL[i].name] = self.axisL[i][im]
        return D
    
    def get_min_max(self):
        return np.nanmin(self.matValArr), np.nanmax(self.matValArr)
        
    def get_min(self):
        return np.nanmin(self.matValArr)
        
    def get_max(self):
        return np.nanmax(self.matValArr)
        
    def get_sum(self):
        return np.nansum(self.matValArr)
        
    def get_ave(self):
        return np.average(self.matValArr)
        
    def get_mean(self):
        return np.mean(self.matValArr)
        
    def get_std(self):
        return np.std(self.matValArr)
        
    def get_median(self):
        return np.median(self.matValArr)

    def get_range(self): # returns max - min
        return np.ptp(self.matValArr) # peak to peak
    
    def __len__(self):
        return len( self.matValArr )
        
    def shape(self):
        return self.matValArr.shape
    
    def size(self):
        return np.prod( self.matValArr.shape )
        
    def iter_indeces(self): # an iterator over the indeces of the matrix
        for indeces in itertools.product(*(list(range(s)) for s in self.matValArr.shape)):
            yield indeces
        
    def iter_items(self): # iterator returns indeces and value at location
        for indeces in itertools.product(*(list(range(s)) for s in self.matValArr.shape)):
            val = self.matValArr[indeces]
            yield indeces,val
        
    def full_iter_items(self): # iterator returns indeces, value and axes value dictionary
        self.axisNameL
        for indeces in itertools.product(*(list(range(s)) for s in self.matValArr.shape)):
            val = self.matValArr[indeces]
            D={}
            for i,axname in enumerate( self.axisNameL ):
                D[axname] = self.axisL[ i ][indeces[i]]
            yield indeces,D,val

    def clone(self):
        Mclone = copy.deepcopy( self )
        Mclone.name = self.name + '(clone)'
        return Mclone
    
    def __neg__(self):
        Mclone = self.clone()
        Mclone.matValArr = np.negative( Mclone.matValArr )
        return Mclone
        
        #return self * (-1.0)
        
    def __abs__(self):
        Mclone = self.clone()
        Mclone.matValArr = abs(self.matValArr)
        return Mclone
    
    def __add__(self, other):
        Mclone = self.clone()
        if isinstance(other, Matrix):
            Mclone.name = self.name + ' + %s'%other.name
            Mclone.matValArr = self.matValArr + other.matValArr
        else:
            Mclone.name = self.name + ' + %s'%other
            Mclone.matValArr = self.matValArr + other
        return Mclone
        
    def __radd__(self, other):
        return self.__add__(other)
        
    def __iadd__(self, other):
        if isinstance(other, Matrix):
            self.name = self.name + ' + %s'%other.name
            self.matValArr = self.matValArr + other.matValArr
        else:
            self.name = self.name + ' + %s'%other
            self.matValArr = self.matValArr + other
        return self
        
    def __sub__(self, other):
        Mclone = self.clone()
        if isinstance(other, Matrix):
            Mclone.name = self.name + ' - %s'%other.name
            Mclone.matValArr = self.matValArr - other.matValArr
        else:
            Mclone.name = self.name + ' - %s'%other
            Mclone.matValArr = self.matValArr - other
        return Mclone
        
    def __rsub__(self, other):
        Mclone = self.clone()
        Mclone.matValArr = np.negative( Mclone.matValArr )
        return Mclone + other
        
    def __isub__(self, other):
        
        if isinstance(other, Matrix):
            self.name = self.name + ' - %s'%other.name
            self.matValArr = self.matValArr - other.matValArr
        else:
            self.name = self.name + ' - %s'%other
            self.matValArr = self.matValArr - other
        return self
        
    def __mul__(self, other):
        Mclone = self.clone()
        if isinstance(other, Matrix):
            Mclone.name = self.name + ' * %s'%other.name
            Mclone.matValArr = self.matValArr * other.matValArr
        else:
            Mclone.name = self.name + ' * %s'%other
            Mclone.matValArr = self.matValArr * other
        return Mclone

    def __rmul__(self, other):
        return self * other

    def __imul__(self, other):
        if isinstance(other, Matrix):
            self.name = self.name + ' * %s'%other.name
            self.matValArr = self.matValArr * other.matValArr
        else:
            self.name = self.name + ' * %s'%other
            self.matValArr = self.matValArr * other
        return self
        
    def __div__(self, other):
        Mclone = self.clone()
        if isinstance(other, Matrix):
            Mclone.name = self.name + ' / %s'%other.name
            Mclone.matValArr = self.matValArr / other.matValArr
        else:
            Mclone.name = self.name + ' / %s'%other
            Mclone.matValArr = self.matValArr / other
        return Mclone

    def __rdiv__(self, other):
        Mclone = self.clone()
        Mclone.matValArr = np.reciprocal( Mclone.matValArr )
        return Mclone * other

    def __idiv__(self, other):
        #print ' plain div'
        if isinstance(other, Matrix):
            self.name = self.name + ' / %s'%other.name
            self.matValArr = self.matValArr / other.matValArr
        else:
            self.name = self.name + ' / %s'%other
            self.matValArr = self.matValArr / other
        return self

    def __truediv__(self, other): # assumes from __future__ import division
        return self.__div__(other)
        
    def __rtruediv__(self, other): # assumes from __future__ import division
        return self.__rdiv__(other)

    def __itruediv__(self, other): # assumes from __future__ import division
        #print 'truediv'
        return self.__idiv__(other)
        
    def __pow__(self, other):
        Mclone = self.clone()
        if isinstance(other, Matrix):
            Mclone.name = self.name + ' ** %s'%other.name
            Mclone.matValArr = self.matValArr ** other.matValArr
        else:
            Mclone.name = self.name + ' ** %s'%other
            Mclone.matValArr = self.matValArr ** other
        return Mclone

    def __rpow__(self, other):
        Mclone = self.clone()
        Mclone.matValArr = (Mclone.matValArr*0.0) + other
        return Mclone**self

    def __ipow__(self, other):
        #print ' plain div'
        if isinstance(other, Matrix):
            self.name = self.name + ' ** %s'%other.name
            self.matValArr = self.matValArr ** other.matValArr
        else:
            self.name = self.name + ' ** %s'%other
            self.matValArr = self.matValArr ** other
        return self
        
    def get_sub_matrix(self, **kwds): # kwds contains axis names... returns val
        '''Call as: M.get_sub_matrix(pc=100, eps=20, mr=2.0)
           Return a smaller Matrix at specified values in kwds'''
        
        is_in_cutL=[0 for axname in self.axisNameL] # set to 1 if axname is a cut plane
        orig_indexL = is_in_cutL[:] # hold index into axis for input axis value
        newAxisNameL = [] # smaller list of axis names in new, smaller Matrix
        for ia,axname in enumerate(self.axisNameL):
            if axname in kwds:
                is_in_cutL[ia]=1
                
                # Also hold const index in cut axis
                orig_indexL[ia] = self.axisL[ia].getExactIndex( kwds[axname] ) 
            else:
                newAxisNameL.append( axname )
        
        #print 'is a slice plane =',is_in_cutL
        #print 'Index of slice plane =',orig_indexL
        
        new_name = self.name +'_'+ '_'.join( ['%s=%s'%(n,v) for n,v in list(kwds.items())] )
        M = Matrix( {'name':new_name,  'units':self.units, 
            'axisNameL':newAxisNameL, 'axisPoolObj':self.axisPoolObj} )
        
        # TODO: change to faster numpy slicing method.
        for new_indeces in M.iter_indeces():
            inew = 0
            for i,ia in enumerate(is_in_cutL):
                if ia==0: # if axis in new Matrix, iterate indeces
                    orig_indexL[i] = new_indeces[inew]
                    inew += 1
            M[ tuple(new_indeces) ] =  self.matValArr[ tuple(orig_indexL) ]
        
        return M
    
    def values_in_range(self, **kwds):
        for k,v in kwds.items():
            A = self.get_axis_by_name( k )
            if not A.value_in_range( v ):
                return False
        return True
    
    def get_axis_by_name(self, aname):
        for a in self.axisL:
            if a.name == aname:
                return a
        return None
    
    def matrix_axis_name_list(self):
        return [a.name for a in self.axisL]
    
    def is_axis_name(self, axis_name):
        return axis_name in self.matrix_axis_name_list()

    
    def get_indeces_where(self, if_gt=0.0, if_lt=None):
        """Return indeces of values in range. Ignore if set to None"""
        if if_lt is None:
            return np.argwhere( self.matValArr > if_gt )
        elif if_gt is None:
            return np.argwhere( self.matValArr < if_lt )
        else:
            return np.argwhere( self.matValArr < if_lt and self.matValArr > if_gt )
    
    def get_dict_of_indeces(self, indeces):
        D={}
        for i,axname in enumerate( self.axisNameL ):
            D[axname] = self.axisL[ i ][indeces[i]]
        return D
    
    def fill_missing_from_good_neighbors(self, bad_if_lt=0.0, bad_if_gt=None,
                                         good_if_gt=0.0, good_if_lt=None):
                                             
        badL = self.get_indeces_where( if_gt=bad_if_gt, if_lt=bad_if_lt)
        
        for iL in badL:
            good_ivL = self.get_nearest_good_neighbors(iL, good_if_gt=good_if_gt, good_if_lt=good_if_lt)
            
            sum_wts = 0.0 # sum of data pt weights
            sum_val_x_wts = 0.0 # sum of value * wt
            
            for good_iv in good_ivL:
                good_indeces = good_iv[0]
                dist = sum( [ abs(i1-i2) for (i1,i2) in zip(iL,good_indeces) ] )
                #print(dist, iL, good_indeces)
                #print(iL, good_indeces, dist)
                if dist > 0:
                    wt = 1.0/float(dist)
                    sum_wts += wt
                    sum_val_x_wts += wt * self[ good_indeces ]
            if sum_wts > 0.0:
                new_val = sum_val_x_wts / sum_wts
                
                #iD = self.get_dict_of_indeces(iL)
                #iD['val'] = new_val
                #self.setByName( **iD )
                self[ iL ] = new_val
                
                #print(iL,'new_val',new_val, self.get_dict_of_indeces(iL), self[iL])
                #for good_iv in good_ivL:
                #    print('    ',good_iv, self.get_dict_of_indeces(good_iv[0]), self[good_iv[0]])
    
    def get_nearest_good_neighbors(self, iL, good_if_gt=0.0, good_if_lt=None ):
        """Return the indeces of nearest good neighbors"""
        def is_good_val( val ):
            
            if good_if_gt is None:
                return val < good_if_lt
            elif good_if_lt is None:
                return val > good_if_gt
            else:
                return val > good_if_gt and val < good_if_lt
        
        
        iL = list( iL ) # makes a list copy
        good_ivL = [] # list of tuples  (indeces, val)
        for ia, i in enumerate( iL ):
            a = self.axisL[ia]
            
            itestL = iL[:]
            j = i+1
            while j < len( a ):
                itestL[ia] = j
                if is_good_val( self[ itestL ] ):
                    good_ivL.append( (itestL, self[itestL]) )
                    j += len(a)
                j += 1
            
            itestL = iL[:]
            j = i-1
            while j >= 0:
                itestL[ia] = j
                if is_good_val( self[ itestL ] ):
                    good_ivL.append( (itestL, self[itestL]) )
                    j -= len(a)
                j -= 1
            
        return good_ivL

    
    def interp_missing_from_good_neighbors(self, bad_if_lt=0.0, bad_if_gt=None,
                                           good_if_gt=0.0, good_if_lt=None):
                                             
        badL = self.get_indeces_where( if_gt=bad_if_gt, if_lt=bad_if_lt)
        print( "Replacing %i bad values from Nearest Neighbors in"%len(badL), self.name )
        
        for iL in badL:
            good_ivL = self.get_nearest_good_neighbors(iL, good_if_gt=good_if_gt, good_if_lt=good_if_lt)
            
            sum_wts = 0.0 # sum of data pt weights
            sum_val_x_wts = 0.0 # sum of value * wt
            
            for good_iv in good_ivL:
                good_indeces = good_iv[0]
                dist = sum( [ abs(i1-i2) for (i1,i2) in zip(iL,good_indeces) ] )
                #print(dist, iL, good_indeces)
                #print(iL, good_indeces, dist)
                if dist > 0:
                    wt = 1.0/float(dist)
                    sum_wts += wt
                    sum_val_x_wts += wt * self[ good_indeces ]
            if sum_wts > 0.0:
                new_val = sum_val_x_wts / sum_wts
                
                #iD = self.get_dict_of_indeces(iL)
                #iD['val'] = new_val
                #self.setByName( **iD )
                self[ iL ] = new_val
                
                #print(iL,'new_val',new_val, self.get_dict_of_indeces(iL), self[iL])
                #for good_iv in good_ivL:
                #    print('    ',good_iv, self.get_dict_of_indeces(good_iv[0]), self[good_iv[0]])
    
    def get_1d_interp_fill_value(self, i_centerL, good_if_gt=0.0, good_if_lt=None):
        """
        Given the indeces, i_centerL, of a point in the matrix, M, return all
        of the 1D matrices with GOOD values as defined by good_if_gt and good_if_lt.
        """
        
        valueL = [] # list of interpolated values (will take average at end)
        
        for ia,a in enumerate(self.axisL):

            # start with a fresh center list
            cL = list( i_centerL )
            aL = [] # good axis value list
            vL = [] # good value list
            
            for i in range( len(a) ):
                cL[ia] = i
                val = self[ cL ]
            
                if good_if_gt is None:
                    if val < good_if_lt:
                        aL.append( a.transObj( a.valueL[i] ) )
                        vL.append( val )
                elif good_if_lt is None:
                    if val > good_if_gt:
                        aL.append( a.transObj( a.valueL[i] ) )
                        vL.append( val )
                else:
                    if val > good_if_gt and val < good_if_lt:
                        aL.append( a.transObj( a.valueL[i] ) )
                        vL.append( val )
                        
            if aL:
                terp = InterpProp(aL, vL, extrapOK=1, linear=1)
                valueL.append( terp( a.transObj( a.valueL[ i_centerL[ia] ] ) ) )
                #print('  val:',val,'  terpVal:',valueL[-1], 'aL:',aL,'  vL:',vL)

            if valueL:
                val = sum(valueL) / len(valueL)
            else:
                val = self[ i_centerL ]
        return val
    
    
    def fill_missing_from_1d_interp(self, bad_if_lt=0.0, bad_if_gt=None,
                                         good_if_gt=0.0, good_if_lt=None):
                                             
        badL = self.get_indeces_where( if_gt=bad_if_gt, if_lt=bad_if_lt)
        print( "1D Interpolating %i bad values in"%len(badL), self.name )
        
        new_valD = {} # index:bad_indeces, value:val
        for iL in badL:
            val = self.get_1d_interp_fill_value( iL, good_if_gt=good_if_gt, good_if_lt=good_if_lt)
            new_valD[ tuple(iL) ] = val
            
        for k,v in new_valD.items():
            self[ k ] = v
        
        # Just in case interpolation fails, use good neighbors
        self.interp_missing_from_good_neighbors( bad_if_lt=bad_if_lt, bad_if_gt=bad_if_gt,
                                           good_if_gt=good_if_gt, good_if_lt=good_if_lt)
    
    def plot_x_param(self, xVar='', param='', fixedD=None, 
                     interp_pts=0, interp_order=2,
                     is_semilog=False, marker='o', markersize=0,
                     rev_param_order=False, show_grid=True,
                     min_val=float('-inf'), max_val=float('inf')):
        """
        Make a plot of xVar vs matrix value, parameterized by param.
        If left blank, use names of 1st two axes.
        Set any other axis values based on fixedD.
        If not in fixedD, then use median value of axis.
        If interp_pts>0, insert interpolated points between axis pts
        """
        #self.axisL self.matValArr
        if len( self.axisL ) < 2:
            print('ERROR... can not make plot_x_param with less than 2 axes.')
            return
            
        if not got_pylab:
            print('ERROR... pylab FAILED to import.')
            return
        
        # if xVar not input, set it to one of 1st 2 axes
        if not self.is_axis_name(xVar):
            if param != self.axisL[0].name:
                xVar  = self.axisL[0].name
            else:
                xVar  = self.axisL[1].name
                
        # if param not input, set it to one of 1st 2 axes
        if not self.is_axis_name(param):
            if xVar  != self.axisL[0].name:
                param = self.axisL[0].name
            else:
                param = self.axisL[1].name        
        #print('xVar=%s,  param=%s'%(xVar, param))
        
        xAxis = self.get_axis_by_name( xVar )
        pAxis = self.get_axis_by_name( param )
        #print( 'xAxis =',xAxis )
        #print( 'pAxis =',pAxis )
        changing_paramL = [xVar, param]
        
        # prepare fixedD of constant values
        fixed_paramL = []
        if fixedD is None:
            D = {}
        else:
            D = fixedD.copy()
            
        sL = [] # making title string of fixed values
        fixedD = {}
        for a in self.axisL:
            if a.name not in D:
                D[a.name] = a.get_middle_value()
            if a.name not in changing_paramL:
                fixed_paramL.append( a.name )
                sL.append( '%s=%g %s'%(a.name, D[a.name], a.units) )
                fixedD[a.name] = D[a.name]
        fixed_s = ', '.join(sL)
        #print( "D=", D,  '   fixedD',fixedD )
        #print( 'fixed_paramL',fixed_paramL, '   fixed_s',fixed_s )
        #print( 'changing_paramL',changing_paramL )
        
        # .......... get sub-matrix to speed things up ..................
        
        SP = self.get_sub_matrix( **fixedD )

        # ================= start making plot ========================
        if rev_param_order:
            paramL = reversed( pAxis.valueL )
        else:
            paramL = pAxis.valueL
        
        
        pylab.figure()
        markerD = {} # matplotlib options
        if marker:
            markerD['marker'] = '.'
            markerD['markevery'] = 1 + interp_pts
        if markersize:
            markerD['markersize'] = markersize
        
        # .... begin iterating over param and xVar
        for p in paramL:
            fL = []
            xL = []
            for x in xAxis.valueL:
                D[ xVar ] = x
                D[ param ] = p
                
                if interp_pts:
                    if x in xAxis.valueL:
                        f = SP.getByName( **D )
                    else:
                        f = SP.interp( order=interp_order, **D)
                        
                else:
                    f = SP.getByName( **D )
                    
                if f is not None and ( min_val <= f <= max_val):
                    fL.append( f )
                    xL.append( x )
                    
            if xL:
                if interp_pts:
                    # make a transformed list of x's for interpolation
                    xtL = [ xAxis.transObj(x) for x in xL]
                    
                    # make full xvarL for interpolation
                    xvarL = xAxis.valueL[:] # make a copy... it will be modified
                    
                    f = 1.0/(1.0 + interp_pts)
                    for i in range( len(xL) - 1 ):
                        for j in range( interp_pts ):
                            xvarL.append( xL[i] + f*(j+1) * (xL[i+1] - xL[i]) )
                    xL = sorted( xvarL )
                    
                    fL = [ interp1d( xtL , fL, kind=interp_order, fill_value="extrapolate")\
                           ( xAxis.transObj(x) ) for x in xL]
                
                    
                
                
                # Assume there are some interpolated points... plot twice.
                if is_semilog:
                    lastplot = pylab.semilogx(xL, fL,  label='%s=%g'%(param, p), **markerD)
                    c = lastplot[0].get_color()
                    pylab.semilogx(xL, fL, linestyle='None', marker='|', color=c)
                else:
                    lastplot = pylab.plot(xL, fL, label='%s=%g'%(param, p), **markerD)
                    c = lastplot[0].get_color()
                    pylab.plot(xL, fL, linestyle='None', marker='|', color=c)

        pylab.title( '%s: %s'%(self.name, fixed_s) )
        pylab.legend(loc='best', framealpha=0.3)
        def axis_label( a ):
            if a.units:
                return '%s (%s)'%(a.name, a.units)
            else:
                return a.name
        
        if show_grid:
            pylab.grid()
        pylab.xlabel( axis_label( xAxis ) )
        pylab.ylabel( self.name )
        
        

if __name__=="__main__":

    epsAxis = Axis({'name':'eps', 'valueL':[10., 20., 30., 40.], 
            'units':'', 'transform':'log10'})
    
    # Just a dict, not an Axis obj
    pcAxis = {'name':'pc', 'valueL':[100.,200.,300], 'units':'psia', 'transform':'log10'}

    mrAxis = Axis({'name':'mr', 'valueL':[1,2,3,4,5], 
            'units':'', 'transform':''})

    axesDefL = [epsAxis, pcAxis, mrAxis]
    
    AP = AxisPool( {'axesDefL':axesDefL} )
        
    axisNameL = ['eps','pc','mr']
    shape = [len(AP.axisD[name])  for name in axisNameL]
    print('shape =',shape)
    matValArr = np.zeros( shape )
    n0,n1,n2 = axisNameL
    for i0,v0 in enumerate(AP.axisD[n0]):
        for i1,v1 in enumerate(AP.axisD[n1]):
            for i2,v2 in enumerate(AP.axisD[n2]):
                matValArr[i0,i1,i2] = v0+v1+v2

    M = Matrix( {'name':'isp_ode', 'matValArr':matValArr, 'units':'', 
        'axisNameL':axisNameL, 'axisPoolObj':AP} )
        
    #print M.axisL
    print(M)
    #print type( M.axisL[0] ) == Axis
    #print type( {1:1} ) == dict
    print(M[(0,0,0)],M[3,2,4],'__getitem__ examples')
    print('_'*55)
    print(mrAxis.matrixConnectionL)
    #epsAxis.add_value( 16.0 )
    j = AP.add_value_to_Axis('pc', 250.0)
    print(M)
    print('   ...Added new axis value.  Matrix expands to accomodate')
    print('_'*55)
    for i in range( len(epsAxis) ):
        for k in range( len(mrAxis) ):
            M[(i,j,k)] = 7777.0
    print(M)
    print('   ...Set inserted value to 7777.  Use index from axis value insert.')
    print('_'*55)
    
    print()
    pc = 250.0
    for eps in epsAxis:
        for mr in mrAxis:
            M.setByName( pc=pc, eps=eps, mr=mr, val=9999.0 )
    print(M)
    print('   ...change 7777 to 9999 using dictionary indexing pc=pc.')
    print('_'*55)
    