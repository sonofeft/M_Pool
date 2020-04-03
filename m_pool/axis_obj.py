#!/usr/bin/env python
# -*- coding: utf8 -*-
import bisect
import numpy as np
from m_pool.transform import Transform

class Axis(object):
    '''An Axis object is used to define the axes of Matrix objects
    
       For interpolation, an axis can transform its values to (hopefully)
         make them more linear.  (e.g. log10 for an axis spanning orders of magnitude)
         
       *** Structured to easily pickle via a dictionary of named values for properties. ***
    '''

    def __init__(self, D={'name':'axisName', 'valueL':None, 'units':'', 
            'transform':None, 'roundDigits':4} ):
        '''Initialize with a dictionary so that pickle files can easily save and read objects'''
        
        self.name = D.get('name','Unknown')
        self.valueL = D.get('valueL', [0.0,1.0])
        self.units = D.get('units','')
        self.transform = D.get('transform','')
        self.transObj = Transform( self.transform )
        self.roundDigits = D.get('roundDigits', 4)
        
        self.valueL = [float(v) for v in self.valueL] # make sure it's all float
        
        self.matrixConnectionL = [] # holds list of tuple objects (iaxis, matrixObj)
        
        self.make_supt_objects()
    
    def get_pickleable_dict(self):
        return {'name':self.name, 'valueL':self.valueL, 'units':self.units, 
                'transform':self.transform, 'roundDigits':self.roundDigits}
    
    def get_valueL(self):
        return self.valueL[:] # return a copy
        
    def get_trans_valueL(self):
        return [ self.transObj(v) for v in self.valueL]
    
    def value_in_range(self, val):
        return val>=self.valueL[0] and val<=self.valueL[-1]
    
    def get_middle_value(self):
        i = int( len(self.valueL) / 2 )
        return self.valueL[i]
    
    def transform_str(self):
        if self.transform:
            return str(self.transform)
        else:
            return 'NONE'
    
    def __str__(self):
        sL = ['Axis: %s'%self.name + ' ' + str(self.valueL) + ' transform:' + self.transform_str() ]
        return '\n'.join(sL)
    
    def __len__(self):
        return len( self.valueL )
    
    def __iter__(self):
        for value in self.valueL:
            yield value
    
    def __getitem__(self, i): # retrieve as: A[i]
        return self.valueL[i]
    
    def connectMatrix(self, iaxis, matrixObj):
        '''Save all connections to Matrix objects'''
        # Save a tuple of index position in Matrix and Matrix itself
        self.matrixConnectionL.append( (iaxis,matrixObj) )
        
    
    def make_supt_objects(self):
        self.valueArr = np.array( self.valueL )
        # Keep a transformed array for later linear interpolation 
        self.transArr = self.transObj( self.valueArr )
        
        self.indexD = {} # quickly look up index from value
        for i,val in enumerate(self.valueL):
            self.indexD[val] = i
            self.indexD[ self.valueArr[i] ] = i # works for both float and np.float
        
    def getQuadterpIndex(self, val):
        '''Index used for quad_terp interpolation'''
        val = float(val)
        i = bisect.bisect_left(self.valueArr, val) - 1
        
        if i<0:
            return 0
        elif i>self.valueArr.size-3:
            return self.valueArr.size-3
        return i
            
    def getExactIndex(self, val):
        '''Get Index into self.valueL from the val input'''
        val = float(val) # make sure it's a float
        try:
            i = self.indexD[val] # if it's in the cross-index, use it
            return i
        except:
            pass

        # Otherwise, see if val is within roundoff error of an included value
        i = bisect.bisect_left( self.valueL, val )
        try:
            #print '   Round A i=%i'%i,round(val, self.roundDigits),round(self.valueL[i], self.roundDigits)
            if round(val, self.roundDigits)==round(self.valueL[i], self.roundDigits):
                return i
        except:
            pass
            
        try:
            #print '   Round B i=%i'%i,round(val, self.roundDigits),round(self.valueL[i-1], self.roundDigits)
            if round(val, self.roundDigits)==round(self.valueL[i-1], self.roundDigits):
                return i-1
        except:
            pass
            
        return None # Can't find the index... return None
        
        
    def add_value(self, val):
        val = float(val)
        i = bisect.bisect_left( self.valueL, val )
        if i>= len(self.valueL):
            self.valueL.append( val )
        
        else:
            if self.valueL[i]==val:
                print(val,'is already in Axis... skipping insert')
                return None
            
            #print 'i=',i,'for',val,'inserted into',self.valueL,'v[i]=',self.valueL[i]
            self.valueL.insert(i,val)
            #print 'new vL=',self.valueL
        self.make_supt_objects()
        
        for iaxis,M in self.matrixConnectionL:
            M.insert_dimension( iaxis,i )
        
        return i # calling routine might need insertion point
        

if __name__=="__main__":

    epsAxis = Axis({'name':'eps', 'valueL':[1.,2.,3.,5.,10.,20.,40.], 'units':'', 'transform':'log10'})
    print(epsAxis.get_pickleable_dict())
    print()
    print('epsAxis.transArr =',epsAxis.transArr)
    print()
    epsAxis.add_value( 1.6 )
    epsAxis.add_value( 1.5 )
    epsAxis.add_value( 7 )
    epsAxis.add_value( 40.0001 )
    epsAxis.add_value( 0.0 )
    print(epsAxis.valueArr)
    print('len(epsAxis) =',len(epsAxis))
    print()
    for val in [0, -1.1, 0.5, 1, 1.25, 4,40, 40.0001, 40.00001, 40.000001]:
    
        i = epsAxis.getExactIndex(val)
        if i==None:
            print('for val=',val,' i=',i)
        else:
            print('for val=',val,' i=%i'%i, ' valueL[i]=',epsAxis.valueL[i])
    
    