import numpy as np
from axis_obj import Axis

def axis_obj_dammit( axOrD ):
    '''No Matter What, return an Axis obj'''
    if type( axOrD ) == dict:
        return Axis( axOrD ) # make an axis from a dictionary

    if type( axOrD ) == Axis:
        return axOrD # it's already an Axis
        
    return Axis() # a default Axis


class AxisPool(object):
    '''An AxisPool object is a Collection of Axis objects. 
       It is used to define a collection of Matrix objects
         
       *** Structured to easily pickle via a dictionary of named values for properties. ***
    '''

    def __init__(self, D={'axesDefL':None} ):
        '''Initialize with a dictionary so that pickle files can easily save and read objects'''
        
        axesDefL = D.get('axesDefL', None)
        
        if axesDefL:
            self.axisL = [axis_obj_dammit(axOrD) for axOrD in axesDefL]
        else:
            self.axisL = []
        
        # uses self.axisD to get an Axis object by name
        self.make_supt_objects()
    
    def add_axis(self, axOrD):
        self.axisL.append( axis_obj_dammit(axOrD) )
            
        self.make_supt_objects()

    def make_supt_objects(self):
        # use axisD to get an Axis object by name
        self.axisD = {}
        for A in self.axisL:
            self.axisD[A.name] = A
        
    def connectMatrixToAxes(self, matrixObj, axisNameL):
        '''Be sure that axisNameL is the correct order for the axes'''
        for i,name in enumerate(axisNameL):
            self.axisD[name].connectMatrix( i, matrixObj )
                
    def add_value_to_Axis(self, axisName, val):
        A = self.axisD[ axisName ]
        i = A.add_value( val )
        return i
        
    def __str__(self):
        sL = ['AxisPool']
        for A in self.axisL:
            sL.append( str(A) )
        return '\n'.join( sL )

    def __len__(self):
        return len( self.axisL )
    
    def __iter__(self):
        for value in self.axisL:
            yield value
                
    def __getitem__(self, axisName):
        return self.axisD.get( axisName, None)

if __name__=="__main__":

    epsAxis = Axis({'name':'eps', 'valueL':[10., 20., 30., 40., 50.], 
            'units':'', 'transform':'log10'})
    
    # Just a dict, not an Axis obj
    pcAxis = {'name':'pc', 'valueL':[100.,200.,300,400], 'units':'psia', 'transform':'log10'}

    mrAxis = Axis({'name':'mr', 'valueL':[1,2,3], 
            'units':'', 'transform':''})

    axesDefL = [epsAxis, pcAxis]
    
    AP = AxisPool( {'axesDefL':axesDefL} )
        
    print AP
    #print AP.axisD
    print '_'*20,'Add another axis called "mr"','_'*20
    AP.add_axis( mrAxis )
    print AP
    print '_'*20,'Add value 2.5 to axis "mr"','_'*20
    i =AP.add_value_to_Axis( 'mr', 2.5 )
    print 'Added 2.5 at position',i,'in "mr"'
    print AP
