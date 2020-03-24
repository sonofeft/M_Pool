
from m_pool.matrix_pool import MatrixPool

MP = MatrixPool(name='CHECK')
MP.read_from_pickle() # Opens "CHECK_matrix.pool"
print(MP)
print(MP.matrixD)

M = MP.get_matrix_by_name( 'cea_isp' )

print('len(M.shape()) =',len(M.shape()))

D = {'pc':250., 'eps':35.0, 'mr':1.5}
val = M.interp(order=1, **D)
print('type(val) =',type(val))
print('interp_1  =',val)

val = M.interp(order=2, pc=250, eps=35, mr=1.5)
print('interp_2  =',val)

val = M.interp(order=3, **D)
print('interp_3  =',val)

