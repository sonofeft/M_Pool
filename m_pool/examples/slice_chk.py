from m_pool.matrix_pool import MatrixPool
from m_pool.axis_obj import Axis
import itertools

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

print(M)

print('_'*55)

newM = M.get_sub_matrix( pc=200. )
print(newM)

print('_'*55)
newM = M.get_sub_matrix( pc=300., eps=50.0 )
print(newM)
