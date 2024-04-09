==========================
Installation on Linux & Mac OSX 
==========================

To Install TEMPy you need to follow these steps:
 
**Download**

Download TEMPy-1.0 from `download`_
  
**Requirements**

You need to install the following software before installing TEMPy:
 
| 1. `Python`_ 
|    We recommend that you use 2.7.
| 2. `NumPy`_ 
|    We recommend that you use the latest version of NumPy, but v1.6 and later are supported.
| 3. `SciPy`_
|    We recommend that you use the latest version of Scipy, but v0.10.0 and later are supported.
| 4. `Biopython`_ 
| 	We recommend that you use the latest version of BioPython, but v1.58 and later are supported.
 
Default installation of these packages requires root (Linux & Mac OSX) privileges. 

An appropriate version of Python is very likely to be available by default on recent Linux and Mac OSX systems. 
`NumPy`_ , `SciPy`_ and `Biopython`_  modules should be installed in user directories. 
Please read the installation instructions for these packages carefully. 

**Installation**

| tar -xzf TEMPy-1.0.tar.gz	
| cd TEMPy-1.0 
| python setup.py build
| sudo python setup.py install

If you don't have root privileges see alternative installation scheme in `Installing Python Modules`_ .

**Recommendations**

We recommend that you also install `matplotlib`_ for plotting.

.. _NumPy:
   http://www.numpy.org/
.. _SciPy:
   http://www.scipy.org/
.. _Biopython:
   http://biopython.org/
.. _matplotlib:
   http://matplotlib.org/ 
.. _Python:
   http://www.python.org
.. _download:
   http://tempy.ismb.lon.ac.uk/download/TEMPy-1.0.tar.gz
   
.. _Installing Python Modules:
   http://docs.python.org/2/install/
