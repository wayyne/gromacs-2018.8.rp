# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /apps/cmake/3.16.4/bin/cmake

# The command to remove a file.
RM = /apps/cmake/3.16.4/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build

# Include any dependencies generated for this target.
include src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/depend.make

# Include the progress variables for this target.
include src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/progress.make

# Include the compile flags for this target's objects.
include src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/flags.make

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/confio.cpp.o: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/flags.make
src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/confio.cpp.o: ../src/gromacs/fileio/tests/confio.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/confio.cpp.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fileio-test.dir/confio.cpp.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/confio.cpp

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/confio.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fileio-test.dir/confio.cpp.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/confio.cpp > CMakeFiles/fileio-test.dir/confio.cpp.i

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/confio.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fileio-test.dir/confio.cpp.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/confio.cpp -o CMakeFiles/fileio-test.dir/confio.cpp.s

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/readinp.cpp.o: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/flags.make
src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/readinp.cpp.o: ../src/gromacs/fileio/tests/readinp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/readinp.cpp.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fileio-test.dir/readinp.cpp.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/readinp.cpp

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/readinp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fileio-test.dir/readinp.cpp.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/readinp.cpp > CMakeFiles/fileio-test.dir/readinp.cpp.i

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/readinp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fileio-test.dir/readinp.cpp.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/readinp.cpp -o CMakeFiles/fileio-test.dir/readinp.cpp.s

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/tngio.cpp.o: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/flags.make
src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/tngio.cpp.o: ../src/gromacs/fileio/tests/tngio.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/tngio.cpp.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fileio-test.dir/tngio.cpp.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/tngio.cpp

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/tngio.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fileio-test.dir/tngio.cpp.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/tngio.cpp > CMakeFiles/fileio-test.dir/tngio.cpp.i

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/tngio.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fileio-test.dir/tngio.cpp.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests/tngio.cpp -o CMakeFiles/fileio-test.dir/tngio.cpp.s

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.o: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/flags.make
src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.o: ../src/testutils/unittest_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/unittest_main.cpp

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/unittest_main.cpp > CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.i

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/unittest_main.cpp -o CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.s

# Object files for target fileio-test
fileio__test_OBJECTS = \
"CMakeFiles/fileio-test.dir/confio.cpp.o" \
"CMakeFiles/fileio-test.dir/readinp.cpp.o" \
"CMakeFiles/fileio-test.dir/tngio.cpp.o" \
"CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.o"

# External object files for target fileio-test
fileio__test_EXTERNAL_OBJECTS =

bin/fileio-test: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/confio.cpp.o
bin/fileio-test: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/readinp.cpp.o
bin/fileio-test: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/tngio.cpp.o
bin/fileio-test: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/__/__/__/testutils/unittest_main.cpp.o
bin/fileio-test: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/build.make
bin/fileio-test: lib/libtestutils.a
bin/fileio-test: lib/libgromacs_rp.so.3.5.0
bin/fileio-test: lib/libgmock.a
bin/fileio-test: src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable ../../../../bin/fileio-test"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fileio-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/build: bin/fileio-test

.PHONY : src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/build

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/clean:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests && $(CMAKE_COMMAND) -P CMakeFiles/fileio-test.dir/cmake_clean.cmake
.PHONY : src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/clean

src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/depend:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++ /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/fileio/tests /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/gromacs/fileio/tests/CMakeFiles/fileio-test.dir/depend

