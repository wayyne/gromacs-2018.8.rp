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

# Utility rule file for doxygen-lib.

# Include the progress variables for this target.
include docs/doxygen/CMakeFiles/doxygen-lib.dir/progress.make

docs/doxygen/CMakeFiles/doxygen-lib: docs/doxygen/doxygen-lib-timestamp.txt


docs/doxygen/doxygen-lib-timestamp.txt: docs/doxygen/doxygen-source-timestamp.txt
docs/doxygen/doxygen-lib-timestamp.txt: docs/doxygen/Doxyfile-version
docs/doxygen/doxygen-lib-timestamp.txt: docs/doxygen/Doxyfile-lib
docs/doxygen/doxygen-lib-timestamp.txt: docs/doxygen/dep-graphs-dot-timestamp.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating library documentation with Doxygen"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E make_directory /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/depgraphs
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -DDOCTYPE=lib -P RunDoxygen.cmake
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E touch /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/doxygen-lib-timestamp.txt

docs/doxygen/Doxyfile-version: ../docs/doxygen/Doxyfile-version.cmakein
docs/doxygen/Doxyfile-version: VersionInfo.cmake
docs/doxygen/Doxyfile-version: ../cmake/gmxConfigureVersionInfo.cmake
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating Doxyfile-version"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -D VERSION_VARIABLES=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/VersionInfo.cmake -D VERSION_CMAKEIN=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/docs/doxygen/Doxyfile-version.cmakein -D VERSION_OUT=Doxyfile-version -P /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/cmake/gmxConfigureVersionInfo.cmake
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E touch Doxyfile-version

docs/doxygen/dep-graphs-dot-timestamp.txt: docs/doxygen/doxygen-xml-timestamp.txt
docs/doxygen/dep-graphs-dot-timestamp.txt: ../docs/doxygen/doxygenxml.py
docs/doxygen/dep-graphs-dot-timestamp.txt: ../docs/doxygen/gmxtree.py
docs/doxygen/dep-graphs-dot-timestamp.txt: ../docs/doxygen/graphbuilder.py
docs/doxygen/dep-graphs-dot-timestamp.txt: ../docs/doxygen/cycle-suppressions.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating include dependency graphs for dot"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /bin/python2.7 /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/docs/doxygen/graphbuilder.py -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++ -B /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build --ignore-cycles /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/docs/doxygen/cycle-suppressions.txt -o /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/depgraphs
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E touch /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/dep-graphs-dot-timestamp.txt

docs/doxygen/doxygen-xml-timestamp.txt: docs/doxygen/doxygen-source-timestamp.txt
docs/doxygen/doxygen-xml-timestamp.txt: docs/doxygen/Doxyfile-version
docs/doxygen/doxygen-xml-timestamp.txt: docs/doxygen/Doxyfile-xml
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Extracting Doxygen documentation to XML"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E make_directory /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/depgraphs
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -DDOCTYPE=xml -P RunDoxygen.cmake
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E touch /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/doxygen-xml-timestamp.txt

doxygen-lib: docs/doxygen/CMakeFiles/doxygen-lib
doxygen-lib: docs/doxygen/doxygen-lib-timestamp.txt
doxygen-lib: docs/doxygen/Doxyfile-version
doxygen-lib: docs/doxygen/dep-graphs-dot-timestamp.txt
doxygen-lib: docs/doxygen/doxygen-xml-timestamp.txt
doxygen-lib: docs/doxygen/CMakeFiles/doxygen-lib.dir/build.make

.PHONY : doxygen-lib

# Rule to build all files generated by this target.
docs/doxygen/CMakeFiles/doxygen-lib.dir/build: doxygen-lib

.PHONY : docs/doxygen/CMakeFiles/doxygen-lib.dir/build

docs/doxygen/CMakeFiles/doxygen-lib.dir/clean:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && $(CMAKE_COMMAND) -P CMakeFiles/doxygen-lib.dir/cmake_clean.cmake
.PHONY : docs/doxygen/CMakeFiles/doxygen-lib.dir/clean

docs/doxygen/CMakeFiles/doxygen-lib.dir/depend:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++ /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/docs/doxygen /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/CMakeFiles/doxygen-lib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : docs/doxygen/CMakeFiles/doxygen-lib.dir/depend

