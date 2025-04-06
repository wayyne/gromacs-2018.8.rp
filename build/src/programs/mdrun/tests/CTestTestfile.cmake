# CMake generated Testfile for 
# Source directory: /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/programs/mdrun/tests
# Build directory: /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/programs/mdrun/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(MdrunTests "/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/bin/mdrun-test" "--gtest_output=xml:/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/Testing/Temporary/MdrunTests.xml")
set_tests_properties(MdrunTests PROPERTIES  LABELS "GTest;IntegrationTest" TIMEOUT "120" _BACKTRACE_TRIPLES "/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/TestMacros.cmake;126;add_test;/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/programs/mdrun/tests/CMakeLists.txt;65;gmx_register_gtest_test;/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/programs/mdrun/tests/CMakeLists.txt;0;")
add_test(MdrunMpiTests "/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/bin/mdrun-mpi-test" "-ntmpi" "2" "--gtest_output=xml:/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/Testing/Temporary/MdrunMpiTests.xml")
set_tests_properties(MdrunMpiTests PROPERTIES  LABELS "GTest;IntegrationTest;MpiTest" TIMEOUT "120" _BACKTRACE_TRIPLES "/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/TestMacros.cmake;126;add_test;/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/programs/mdrun/tests/CMakeLists.txt;82;gmx_register_gtest_test;/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/programs/mdrun/tests/CMakeLists.txt;0;")
