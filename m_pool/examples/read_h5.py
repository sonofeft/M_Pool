from m_pool.matrix_pool import MatrixPool

MP = MatrixPool()
MP.read_from_hdf5('Pool Data_mpool.h5')
print(MP.summ())

MP2 = MatrixPool()
MP2.read_from_pickle('CHECK_matrix.pool')
print(MP2.summ())
