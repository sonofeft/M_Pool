
import numpy as np
from math import log10, log, sqrt, pow, exp

def square(x):
    return x*x
    
def pow10(x):
    return pow(10.0, x)

def nppow10(x):
    return np.power(10.0, x)

class Transform(object):
    '''Used to transform a data array into a more linear data array'''
    def __init__(self, name='log10'):
        
        self.name = name
        
        if name=='log10':
            self.func = log10
            self.npfunc = np.log10
            self.r_func = pow10
            self.r_npfunc = nppow10
        elif name=='log':
            self.func = log
            self.npfunc = np.log
            self.r_func = exp
            self.r_npfunc = np.exp
        elif name=='sqrt':
            self.func = sqrt
            self.npfunc = np.sqrt
            self.r_func = square
            self.r_npfunc = square
        elif name in ['',None]:
            self.func = None
            self.npfunc = None
            self.r_func = None
            self.r_npfunc = None
        else:
            print 'WARNING... Unknown Transform "%s"'%name
            self.func = None
            self.npfunc = None
            self.r_func = None
            self.r_npfunc = None
            
    def inverse(self, value):
        '''reverse transform either float value or numpy array and return same type.'''
        
        if self.r_func==None:
            return value
            
        
        if hasattr(value,'__len__'):
            return self.r_npfunc( value )
        else:
            return self.r_func( value )
        
        
    def __call__(self, value):
        '''transform either float value or numpy array and return same type.'''
        
        if self.func==None:
            return value
            
        
        if hasattr(value,'__len__'):
            return self.npfunc( value )
        else:
            return self.func( value )
            
if __name__=="__main__":

    T = Transform( 'log10' )
    
    print 'T(2.2)=',T(2.2), type(T(2.2))
    
    arr = np.array( [1.1, 2.2, 3.3, 4.4] )
    print 'T(arr)=',T(arr)
    
    print
    T2 = Transform('What???')
    print 'T2(2.2)=',T2(2.2)
    print 'T2(arr)=',T2(arr)

