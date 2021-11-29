echo "Cleaning up the CMake output files..."
find . \( -iwholename '*cmake*' -o -name 'Makefile' -o -name 'lib' -o -iwholename '*exec_*' \) -not -name CMakeLists.txt -delete
