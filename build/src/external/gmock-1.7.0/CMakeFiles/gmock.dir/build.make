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
include src/external/gmock-1.7.0/CMakeFiles/gmock.dir/depend.make

# Include the progress variables for this target.
include src/external/gmock-1.7.0/CMakeFiles/gmock.dir/progress.make

# Include the compile flags for this target's objects.
include src/external/gmock-1.7.0/CMakeFiles/gmock.dir/flags.make

src/external/gmock-1.7.0/CMakeFiles/gmock.dir/src/gmock-all.cc.o: src/external/gmock-1.7.0/CMakeFiles/gmock.dir/flags.make
src/external/gmock-1.7.0/CMakeFiles/gmock.dir/src/gmock-all.cc.o: ../src/external/gmock-1.7.0/src/gmock-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/external/gmock-1.7.0/CMakeFiles/gmock.dir/src/gmock-all.cc.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gmock.dir/src/gmock-all.cc.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/external/gmock-1.7.0/src/gmock-all.cc

src/external/gmock-1.7.0/CMakeFiles/gmock.dir/src/gmock-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock.dir/src/gmock-all.cc.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/external/gmock-1.7.0/src/gmock-all.cc > CMakeFiles/gmock.dir/src/gmock-all.cc.i

src/external/gmock-1.7.0/CMakeFiles/gmock.dir/src/gmock-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock.dir/src/gmock-all.cc.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/external/gmock-1.7.0/src/gmock-all.cc -o CMakeFiles/gmock.dir/src/gmock-all.cc.s

src/external/gmock-1.7.0/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o: src/external/gmock-1.7.0/CMakeFiles/gmock.dir/flags.make
src/external/gmock-1.7.0/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o: ../src/external/gmock-1.7.0/gtest/src/gtest-all.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object src/external/gmock-1.7.0/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -Wno-unused-variable -o CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o -c /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/external/gmock-1.7.0/gtest/src/gtest-all.cc

src/external/gmock-1.7.0/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.i"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -Wno-unused-variable -E /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/external/gmock-1.7.0/gtest/src/gtest-all.cc > CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.i

src/external/gmock-1.7.0/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.s"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -Wno-unused-variable -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/external/gmock-1.7.0/gtest/src/gtest-all.cc -o CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.s

# Object files for target gmock
gmock_OBJECTS = \
"CMakeFiles/gmock.dir/src/gmock-all.cc.o" \
"CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o"

# External object files for target gmock
gmock_EXTERNAL_OBJECTS =

lib/libgmock.a: src/external/gmock-1.7.0/CMakeFiles/gmock.dir/src/gmock-all.cc.o
lib/libgmock.a: src/external/gmock-1.7.0/CMakeFiles/gmock.dir/gtest/src/gtest-all.cc.o
lib/libgmock.a: src/external/gmock-1.7.0/CMakeFiles/gmock.dir/build.make
lib/libgmock.a: src/external/gmock-1.7.0/CMakeFiles/gmock.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX static library ../../../lib/libgmock.a"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && $(CMAKE_COMMAND) -P CMakeFiles/gmock.dir/cmake_clean_target.cmake
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gmock.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/external/gmock-1.7.0/CMakeFiles/gmock.dir/build: lib/libgmock.a

.PHONY : src/external/gmock-1.7.0/CMakeFiles/gmock.dir/build

src/external/gmock-1.7.0/CMakeFiles/gmock.dir/clean:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 && $(CMAKE_COMMAND) -P CMakeFiles/gmock.dir/cmake_clean.cmake
.PHONY : src/external/gmock-1.7.0/CMakeFiles/gmock.dir/clean

src/external/gmock-1.7.0/CMakeFiles/gmock.dir/depend:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++ /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/src/external/gmock-1.7.0 /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0 /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/src/external/gmock-1.7.0/CMakeFiles/gmock.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/external/gmock-1.7.0/CMakeFiles/gmock.dir/depend

