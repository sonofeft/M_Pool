from m_pool.matrix_pool import MatrixPool
from m_pool.axis_obj import Axis

MP = MatrixPool(name='CHECK')
epsAxis = Axis({'name':'eps', 'valueL':[10., 20., 30., 40., 50.], 'units':'', 'transform':''})
pcAxis = Axis({'name':'pc', 'valueL':[100.,200.,300,400], 'units':'psia', 'transform':''})
mrAxis = Axis({'name':'mr', 'valueL':[1,2,3,4], 'units':'', 'transform':''})
for A in [epsAxis, pcAxis, mrAxis]:
    MP.add_axis( A )

def chkfunc(pc, eps, mr):
    #return eps + 2.0*pc + 0.1*mr 
    #return -(((eps-25.0)**2) + 2.0*((pc-225.)**2) + 333.333*(mr-1.75)**2)
    return (62267.0-(((eps-25.0)**2) + 2.0*((pc-225.)**2) + 333.333*(mr-1.75)**2))*(300.0/62267.0) + 30.0

M = MP.add_matrix( name='cea_isp', units='sec', axisNameL=['eps','pc','mr'] )
for eps in epsAxis:
    for pc in pcAxis:
        for mr in mrAxis:
            val = chkfunc(pc, eps, mr)
            M.setByName( pc=pc, eps=eps, mr=mr, val=val )



#print M
print('len(M.shape()) =',len(M.shape()))
interpD, max_val = M.solve_interp_max( order=3, method='TNC', tol=1.0E-8)
print('interpD =',interpD)
print('max_val =',max_val)
print('range =',M.get_range())