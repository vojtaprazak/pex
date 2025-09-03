Written by Daven Vasishtan, contributions by Vojta Prazak

Dependencies:

NumPy
sciPy
scikit-image
lxml
matplotlib
mrcfile

add these to your PATH and PYTHONPATH

export PR=/path/to/your/pex
export PYTHONPATH=${PYTHONPATH}:${PR}/TEMPy_py3/build/lib/TEMPy:${PR}/TEMPy_extensions_py3:${PR}/TEMPy_py3/build/lib:${PR}/user_py_libs
export PATH=${PR}/TEMPy_extensions_py3/bin:${PATH}
