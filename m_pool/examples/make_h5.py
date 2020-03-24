from m_pool.matrix_pool import MatrixPool

MP = MatrixPool()
MP.read_from_pickle('CHECK_matrix.pool')
#MP.read_from_pickle(r'D:\win7_CDrive\PROGRAMS\ESP_for_TomF\py_ESP\py_esp\pool_data\N2O4_MMH_matrix.pool')

print(MP.summ())

MP.save_to_hdf5()
