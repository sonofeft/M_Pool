
.. image:: https://travis-ci.org/sonofeft/M_Pool.svg?branch=master

.. image:: https://img.shields.io/pypi/v/M_Pool.svg
        
.. image:: https://img.shields.io/badge/python-3.6|3.7-blue

.. image:: https://img.shields.io/pypi/l/M_Pool.svg

**A Wrapper For Numpy Arrays**

Provides Named Axes, Interpolation, Iteration, Disk Persistence And Numerical Calcs


See the Code at: `<https://github.com/sonofeft/M_Pool>`_

See the Docs at: `<http://m_pool.readthedocs.org/en/latest/>`_

See PyPI page at:`<https://pypi.python.org/pypi/m_pool>`_


M_Pool wraps multidimensional numpy arrays to provide the following features::

    #. MatrixPool objects contain related Axis and Matrix objects
        - MP = MatrixPool(name='N2O4_MMH')
        
    #. Axis objects are added by name and interpolation transform (used to linearize interpolation)
        - epsAxis = Axis({'name':'eps', 'valueL':[10., 20., 30., 40., 50.], 'units':'', 'transform':'log10'})
        - pcAxis = Axis({'name':'pc', 'valueL':[100.,200.,300,400], 'units':'psia', 'transform':'log10'})
        - mrAxis = Axis({'name':'mr', 'valueL':[1,2,3], 'units':'', 'transform':''})
        
    #. Matrix objects added by name 
        - M = MP.add_matrix( name='cea_isp', units='sec', axisNameL=['eps','pc','mr'] )
        
    #. Find interpolated minimum or maximum
        - interpD, max_val = M.solve_interp_max( order=3, method='TNC', tol=1.0E-8)
            - where interpD and max_val look something like:
            - interpD = {'pc': 225.00641803120988, 'eps': 34.991495018803455, 'mr': 1.7499612975876655}
            - max_val = -0.000155216246295
    
    #. Disk-based persistence
        - Save to pickle or hdf5 file
            - MP.save_to_pickle() # saves MP to "N2O4_MMH_matrix.pool"
            
    #. Built-in statistics (standard deviation, median, mean/average, sum, minimum, maximum
        - M.get_range()
        - M.get_ave()
        - M.get_mean()
        - M.get_std()
        - M.get_median()
    
    #. Interpolation on axes with named values
        - interp_val = M.interp(order=2, pc=100, eps=20, mr=2.0)
        - Uses transformed axes to help linearize interpolation
        
    #. Iterate over matrix
        - for indeces,D,val in M.full_iter_items():
            - gives something like:
            - (0, 0, 0) {'pc': 100.0, 'eps': 10.0, 'mr': 1.0} 111.0
            - (0, 0, 1) {'pc': 100.0, 'eps': 10.0, 'mr': 2.0} 112.0
            - (0, 0, 2) {'pc': 100.0, 'eps': 10.0, 'mr': 3.0} 113.0
            - ...

