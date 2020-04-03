from m_pool.matrix_pool import MatrixPool
from m_pool.matrix_obj import Matrix
from m_pool.axis_obj import Axis
from m_pool.axis_pool import AxisPool

from math import log10
from copy import deepcopy

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
# ===================================================================
Mbad = deepcopy( M )
Mbad.setByName( eps=30, pc=200, mr=2, val=-1 )

Mbad.fill_missing_from_1d_interp( bad_if_lt=0.0, bad_if_gt=None,
                                  good_if_gt=0.0, good_if_lt=None)

print(Mbad)




#badL = Mbad.get_indeces_where( if_gt=None, if_lt=0.0)
#print( 'badL:',badL )
#print(Mbad)

#for bad in badL:
#    print('for bad:', bad)
#    val_est = Mbad.get_1d_interp_fill_value( bad, good_if_gt=0.0, good_if_lt=None)
#    print('    val_est:',val_est,'   val:',M[bad])
    