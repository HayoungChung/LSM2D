rm -rf build/lib*
rm lsm_classes.cpp
python setup_class_test.py build
cp build/lib*/lsm_classes.so ./
