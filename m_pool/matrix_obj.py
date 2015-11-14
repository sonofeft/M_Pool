import sys
import itertools
import copy
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize

from axis_obj import Axis
from axis_pool import AxisPool, axis_obj_dammit


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
            print 'ERROR... both axisNameL and axisPoolObj MUST be specified in Matrix'
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
        print res
        
        fun( res.x )# make sure interpD is set
        
        minmax_val = float( self.interp( **interpD ) )
        
        return interpD, minmax_val
        

    def interp(self, order=3, **kwds): # kwds contains axis names... returns interpolated val
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
            for mindeces in itertools.product(*(range(s) for s in m.shape)):
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
                    m[mindeces] = interp1d( A.transArr , yL, kind=kind)(xval)
                except:
                    print 'Extrapolating',A.name,'axis =',A.transArr,'  xval=',xval
                    print '           yL =',yL
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
            result = interp1d( A.transArr, m, kind=kind)( xval )
        except:
            print 'Extrapolating','axis =',A,'  xval=',xval
            print '           m =',m
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

    def short_summ(self):
        if self.matValArr is None:
            sL = ['Matrix %s (shape=%s) %s (units=%s)'%(self.name, str(self.matValArr),self.name, self.units)]
        else:
            sL = ['Matrix %s (shape=%s) %s (units=%s)'%(self.name, str(self.matValArr.shape),self.name, self.units)]
        for A in self.axisL:
            sL.append( str(A) )
        
        return '\n'.join( sL )

    def __str__(self):
        if self.matValArr is None:
            sL = ['Matrix %s (shape=%s) %s (units=%s)'%(self.name, str(self.matValArr),self.name, self.units)]
        else:
            sL = ['Matrix %s (shape=%s) %s (units=%s)'%(self.name, str(self.matValArr.shape),self.name, self.units)]
        for A in self.axisL:
            sL.append( str(A) )
        
        sL.append( str(self.matValArr) )
        return '\n'.join( sL )
        
    def __getitem__(self, iL):
        return self.matValArr[ tuple(iL) ]

        
    def __setitem__(self, iL, val): # use as M[(i,j,k)] = val
        if val is None:
            print 'ERROR... illegal value for "val" in Matrix.set.  val =',val
        else:
            self.matValArr[iL] = float(val)

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
        
    def get_peak_indeces(self):
        imax = np.unravel_index(self.matValArr.argmax(), self.matValArr.shape)
        return imax
    
    def get_minima_indeces(self):
        imin = np.unravel_index(self.matValArr.argmin(), self.matValArr.shape)
        return imin
    
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
        for indeces in itertools.product(*(range(s) for s in self.matValArr.shape)):
            yield indeces
        
    def iter_items(self): # iterator returns indeces and value at location
        for indeces in itertools.product(*(range(s) for s in self.matValArr.shape)):
            val = self.matValArr[indeces]
            yield indeces,val
        
    def full_iter_items(self): # iterator returns indeces, value and axes value dictionary
        self.axisNameL
        for indeces in itertools.product(*(range(s) for s in self.matValArr.shape)):
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
        
        return self * (-1.0)
        
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
            if kwds.has_key(axname):
                is_in_cutL[ia]=1
                
                # Also hold const index in cut axis
                orig_indexL[ia] = self.axisL[ia].getExactIndex( kwds[axname] ) 
            else:
                newAxisNameL.append( axname )
        
        #print 'is a slice plane =',is_in_cutL
        #print 'Index of slice plane =',orig_indexL
        
        new_name = self.name +'_'+ '_'.join( ['%s=%s'%(n,v) for n,v in kwds.items()] )
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
    print 'shape =',shape
    matValArr = np.zeros( shape )
    n0,n1,n2 = axisNameL
    for i0,v0 in enumerate(AP.axisD[n0]):
        for i1,v1 in enumerate(AP.axisD[n1]):
            for i2,v2 in enumerate(AP.axisD[n2]):
                matValArr[i0,i1,i2] = v0+v1+v2

    M = Matrix( {'name':'isp_ode', 'matValArr':matValArr, 'units':'', 
        'axisNameL':axisNameL, 'axisPoolObj':AP} )
        
    #print M.axisL
    print M
    #print type( M.axisL[0] ) == Axis
    #print type( {1:1} ) == dict
    print M[(0,0,0)],M[3,2,4],'__getitem__ examples'
    print '_'*55
    print mrAxis.matrixConnectionL
    #epsAxis.add_value( 16.0 )
    j = AP.add_value_to_Axis('pc', 250.0)
    print M
    print '   ...Added new axis value.  Matrix expands to accomodate'
    print '_'*55
    for i in xrange( len(epsAxis) ):
        for k in xrange( len(mrAxis) ):
            M[(i,j,k)] = 7777.0
    print M
    print '   ...Set inserted value to 7777.  Use index from axis value insert.'
    print '_'*55
    
    print
    pc = 250.0
    for eps in epsAxis:
        for mr in mrAxis:
            M.setByName( pc=pc, eps=eps, mr=mr, val=9999.0 )
    print M
    print '   ...change 7777 to 9999 using dictionary indexing pc=pc.'
    print '_'*55
    