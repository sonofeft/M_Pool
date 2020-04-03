from m_pool.matrix_pool import MatrixPool
from m_pool.matrix_obj import Matrix
from m_pool.axis_obj import Axis
from m_pool.axis_pool import AxisPool
from math import log10

MP = MatrixPool(name='CHECK')
epsAxis = Axis({'name':'eps', 'valueL':[10., 20., 30., 40., 50.], 'units':'', 'transform':'log10'})
pcAxis = Axis({'name':'pc', 'valueL':[100.,200.,300,400], 'units':'psia', 'transform':'log10'})
mrAxis = Axis({'name':'mr', 'valueL':[1,2,3], 'units':'', 'transform':''})
for A in [epsAxis, pcAxis, mrAxis]:
    MP.add_axis( A )

def chkfunc(pc, eps, mr):
    return log10(eps)**2 + 2.0*log10(pc)**2 + 0.1*mr 

M = MP.add_matrix( name='cea_isp', units='sec', axisNameL=['eps','pc','mr'] )
for eps in epsAxis:
    for pc in pcAxis:
        for mr in mrAxis:
            val = chkfunc(pc, eps, mr)
            M.setByName( pc=pc, eps=eps, mr=mr, val=val )

middleD = {}
for a in M.axisL:
    middleD[a.name] = a.get_middle_value()
print('middleD:', middleD)
# for each axis, make a small 1D matrix and test interp
for a in M.axisL:
    small_avalL = [a.valueL[0], a.valueL[-1]]
    print('small_avalL:',small_avalL)
    Asmall = Axis( {'name':a.name, 'valueL':small_avalL, 
                    'units':a.units, 'transform':a.transform} )
    print( Asmall )
    print('1'*55)
    axesDefL = [Asmall]
    AP = AxisPool( {'axesDefL':axesDefL} )
    print(AP)
    print('2'*55)
    Msmall = Matrix( {'axisNameL':[a.name], 'axisPoolObj':AP, 'name':'Msmall'} )

    for aval in small_avalL:    
        mD = middleD.copy()
        mD[a.name] = aval
        print('mD:',mD)
        val = M.getByName( **mD )
        Msmall.setByName( **{a.name:aval, 'val':val} )
    print(Msmall)
    
    # ... do interp
    for aval in a.valueL[1:-1]:
        iD = {a.name:aval, 'order':1}
        vsmall = Msmall.interp(**iD)
        mD = middleD.copy()
        mD[a.name] = aval
        vorig = M.interp(**middleD)
        print( 'at:%g, Msmall=%g,  M=%g'%(aval, vsmall, vorig  ) )
    
    print('#'*55)
