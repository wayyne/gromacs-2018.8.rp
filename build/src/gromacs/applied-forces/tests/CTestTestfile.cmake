# CMake generated Testfile for 
# Source directory: /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/applied-forces/tests
# Build directory: /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/applied-forces/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(AppliedForcesUnitTest "/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/bin/applied-forces-test" "--gtest_output=xml:/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/Testing/Temporary/AppliedForcesUnitTest.xml")
set_tests_properties(AppliedForcesUnitTest PROPERTIES  LABELS "GTest;UnitTest" TIMEOUT "30" _BACKTRACE_TRIPLES "/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/TestMacros.cmake;126;add_test;/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/TestMacros.cmake;136;gmx_register_gtest_test;/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/applied-forces/tests/CMakeLists.txt;35;gmx_add_unit_test;/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/applied-forces/tests/CMakeLists.txt;0;")
