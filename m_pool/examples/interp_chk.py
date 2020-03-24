from m_pool.matrix_pool import MatrixPool
from m_pool.axis_obj import Axis
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
            
MP.save_to_pickle()  # Saves to "CHECK_matrix.pool"

D = {'pc':250., 'eps':35.0, 'mr':1.5}
#val = M.interp(order=2, **{'pc':250., 'eps':35.0, 'mr':1.5})

#print M
print('len(M.shape()) =',len(M.shape()))

val = M.interp(order=1, **D)
print('type(val) =',type(val))
print('interp_1  =',val,'   func =',chkfunc(D['pc'],D['eps'],D['mr']))

val = M.interp(order=2, **D)
print('interp_2  =',val,'   func =',chkfunc(D['pc'],D['eps'],D['mr']))

val = M.interp(order=3, **D)
print('interp_3  =',val,'   func =',chkfunc(D['pc'],D['eps'],D['mr']))

