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
include src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/depend.make

# Include the progress variables for this target.
include src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/progress.make

# Include the compile flags for this target's objects.
include src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/flags.make

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/gputest.cpp.o: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/flags.make
src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/gputest.cpp.o: ../src/gromacs/gpu_utils/tests/gputest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/gputest.cpp.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gpu_utils-test.dir/gputest.cpp.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/gputest.cpp

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/gputest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gpu_utils-test.dir/gputest.cpp.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/gputest.cpp > CMakeFiles/gpu_utils-test.dir/gputest.cpp.i

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/gputest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gpu_utils-test.dir/gputest.cpp.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/gputest.cpp -o CMakeFiles/gpu_utils-test.dir/gputest.cpp.s

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.o: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/flags.make
src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.o: ../src/gromacs/gpu_utils/tests/hostallocator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/hostallocator.cpp

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/hostallocator.cpp > CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.i

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/hostallocator.cpp -o CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.s

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.o: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/flags.make
src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.o: ../src/gromacs/gpu_utils/tests/pinnedmemorychecker.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/pinnedmemorychecker.cpp

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/pinnedmemorychecker.cpp > CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.i

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests/pinnedmemorychecker.cpp -o CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.s

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.o: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/flags.make
src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.o: ../src/testutils/unittest_main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/unittest_main.cpp

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/unittest_main.cpp > CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.i

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/testutils/unittest_main.cpp -o CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.s

# Object files for target gpu_utils-test
gpu_utils__test_OBJECTS = \
"CMakeFiles/gpu_utils-test.dir/gputest.cpp.o" \
"CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.o" \
"CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.o" \
"CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.o"

# External object files for target gpu_utils-test
gpu_utils__test_EXTERNAL_OBJECTS =

bin/gpu_utils-test: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/gputest.cpp.o
bin/gpu_utils-test: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/hostallocator.cpp.o
bin/gpu_utils-test: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/pinnedmemorychecker.cpp.o
bin/gpu_utils-test: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/__/__/__/testutils/unittest_main.cpp.o
bin/gpu_utils-test: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/build.make
bin/gpu_utils-test: lib/libtestutils.a
bin/gpu_utils-test: lib/libgromacs_rp.so.3.5.0
bin/gpu_utils-test: lib/libgmock.a
bin/gpu_utils-test: lib/libgpu_utilstest_cuda.so
bin/gpu_utils-test: src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable ../../../../bin/gpu_utils-test"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gpu_utils-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/build: bin/gpu_utils-test

.PHONY : src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/build

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/clean:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests && $(CMAKE_COMMAND) -P CMakeFiles/gpu_utils-test.dir/cmake_clean.cmake
.PHONY : src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/clean

src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/depend:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++ /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/gromacs/gpu_utils/tests /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/gromacs/gpu_utils/tests/CMakeFiles/gpu_utils-test.dir/depend

