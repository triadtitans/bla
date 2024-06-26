.. ASC-bla documentation master file, created by
   sphinx-quickstart on Tue Aug 29 06:39:02 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to triadtitans-bla's documentation!
===================================

TriadTitans Basic Linear Algebra define a set of fundamental operations on vectors and matrices which can be used to create optimized higher-level linear algebra functionality.
Written in C++ and bound to python, TriadTitans BLA library combines the best of both sides, as it unleashes the full power of your PC while still giving you the comfort of python programming.


The library provides basic template classes **Vector** and **Matrix** and three levels of basic operations:

=========== ===============================================================
**Level 1** Vector operations, e.g. :math:`y = \alpha*x + y`
**Level 2** Matrix-vector operations, e.g. :math:`y = \alpha*A*x + \alpha*beta*y`
**Level 3** Matrix-matrix operations, e.g. :math:`C = \alpha*A*B + C`
=========== ===============================================================

Based on the beasic opertations and together with lapack even more complex operations are implementes possible.
With TriadTitans BLA library you can do

=========== ===============================================================
**Level 4** Inverting, e.g. :math:`C = \alpha*Inv(A)*B` (will soon be available for python)
**Lavel 5** Transposing e.g. :math:`C = \alpha*Transpose(A)*B (will soon be available for python)
**Level 6** LU Factorisation, e.g. :math:`y = \alpha*x + y` (will soon be available for python)
**Level 7** Calculating Eigenvalues, e.g. :math:`eigv = Eigenv(A)` (will soon be available for python)
=========== ===============================================================

In the near future our library will also be able to solve partial and ordinary differential equations!


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

..   code-block:: python

   import pickle

   # write matrix to file
   f = open("file.txt", 'wb')
   pickle.dump(matrix, f)
   del f

   # load matrix from file
   f2 = open("file.txt", 'rb')
   val = pickle.load(f2)
   del f2

   # print loaded matrix
   print (val)

Buffer protocol
---------------

In case conversion of bla objects to numpy arrays is necessary, the
Matrix class is compliant with pythons buffer protocol. If you initialise a
numpy array like this:

.. code-block:: python

   array = np.array(simple_matrix)

The buffer protocol will be used to efficiently pass the matrix to numpy.

.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
