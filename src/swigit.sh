export CFLAGS="-I/Users/marshalltekell/anaconda3/lib/python3.7/site-packages/numpy/core/include $CFLAGS"
swig -python cfunctions.i
python setup.py build_ext --inplace