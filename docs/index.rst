.. ASC-bla documentation master file, created by
   sphinx-quickstart on Tue Aug 29 06:39:02 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to triadtitans-bla's documentation!
===================================

bla is a C++ library for basic linear algebra operations.
The library provides template classes **Vector** and **Matrix**.

Installation is via git-clone:

..  code-block::
    
    git clone https://github.com/TUWien-ASC/ASC-bla.git


To configure and build some tests do

..  code-block::

    cd ASC-bla
    mkdir build
    cd build
    cmake ..
    make
    

To use bla in your code, set the compiler include path properly, and include the header files

..  code-block::

    #include <vector.h>
    #include <matrix.h>

All objects are implemented in the namespace bla. To use them with less typing, you can use

..  code-block::
    
    using namespace bla;

    

You can create vectors and compute with vectors like:

..  code-block:: cpp
                 
   Vector<double> x(5), y(5), z(5);
   for (int i = 0; i < x.Size(); i++)
      x(i) = i;
   y = 5.0
   z = x+3*y;
   cout << "z = " << z << endl;


For matrices you can choose between row-major (`RowMajor`) or column-major (`ColMajor`) storage,
default is row-major.

..  code-block:: cpp

   Matrix<double,RowMajor> m1(5,3), m2(3,3);
   for (int i = 0; i < m1.Height(); i++)
     for (int j = 0; j < m1.Width(); j++)
       m1(i,j) = i+j;
   m2 = 3.7;
   Matrix product = m1 * m2;
   
You can extract a rows or a columns from a matrix:

..  code-block:: cpp

   Vector col1 = product.Col(1);

Python bindings
===============

It is possible to use the library in python to benefit from better performance
when working with matrices. To install the library run:

..   code-block::

   pip install git+https://github.com/triadtitans/bla.git

It then possible to use the library as follows:

..   code-block:: python

   import ASCsoft.bla as bla

   m = bla.Matrix(3,3)
   for p in [(i,j,i+j) for i in range(3) for j in range(3)]:
     (r,c,v) = p
     m[r,c] = v

   print(m)
   print(3*m)
   print(m*m)

Pickling support
----------------

The python objects defined by bla can efficiently be serialized using pickle.
The following example shows how this can be done:

..   code-block:: Python
   <TODO>

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
