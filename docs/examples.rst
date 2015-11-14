
.. examples


.. _internal_examples:

Example Scripts
===============

There are a number of example scripts in the examples subdirectory. 
Every example starts with either creating or reading from disc a MatrixPool.


Build a MatrixPool
------------------


The normal steps to create a MatrixPool are::

    1) Make a MatrixPool object
    2) Make some Axis objects
    3) Add the Axis objects to the pool
    4) Create Matrix objects that use the Axis objects
    5) Set the values of all elements of each Matrix object 
    
Typical creation steps are shown below.
This example creates 3 Axis objects, two of which use a ``log10`` transform to linearize their interpolation.
There is a dummy function used to set the initial values.

Initializing the values in the Matrix object uses named parameters for convenience.

The entire MatrixPool is saved to a file called "CHECK_matrix.pool".

.. code:: python

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
    
    
Read From File
--------------

The above file "CHECK_matrix.pool" can be opened with code like the following.

.. code:: python

    from m_pool.matrix_pool import MatrixPool

    MP = MatrixPool(name='CHECK')
    MP.read_from_pickle() # Opens "CHECK_matrix.pool"
    print MP
    print MP.matrixD
    
Which will generate the following output::

    Reading: CHECK_matrix.pool
    MatrixPool: CHECK
      Matrix:cea_isp, shape=(5, 4, 3), units=sec, %full=100
    {'cea_isp': <m_pool.matrix_obj.Matrix object at 0x052D0D50>}

    