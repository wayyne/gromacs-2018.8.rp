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

# Utility rule file for dep-graphs-dot.

# Include the progress variables for this target.
include docs/doxygen/CMakeFiles/dep-graphs-dot.dir/progress.make

docs/doxygen/CMakeFiles/dep-graphs-dot: docs/doxygen/dep-graphs-dot-timestamp.txt


docs/doxygen/dep-graphs-dot-timestamp.txt: docs/doxygen/doxygen-xml-timestamp.txt
docs/doxygen/dep-graphs-dot-timestamp.txt: ../docs/doxygen/doxygenxml.py
docs/doxygen/dep-graphs-dot-timestamp.txt: ../docs/doxygen/gmxtree.py
docs/doxygen/dep-graphs-dot-timestamp.txt: ../docs/doxygen/graphbuilder.py
docs/doxygen/dep-graphs-dot-timestamp.txt: ../docs/doxygen/cycle-suppressions.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating include dependency graphs for dot"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /bin/python2.7 /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/docs/doxygen/graphbuilder.py -S /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++ -B /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build --ignore-cycles /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/docs/doxygen/cycle-suppressions.txt -o /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/depgraphs
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E touch /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/dep-graphs-dot-timestamp.txt

docs/doxygen/doxygen-xml-timestamp.txt: docs/doxygen/doxygen-source-timestamp.txt
docs/doxygen/doxygen-xml-timestamp.txt: docs/doxygen/Doxyfile-version
docs/doxygen/doxygen-xml-timestamp.txt: docs/doxygen/Doxyfile-xml
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Extracting Doxygen documentation to XML"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E make_directory /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/depgraphs
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -DDOCTYPE=xml -P RunDoxygen.cmake
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E touch /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/doxygen-xml-timestamp.txt

docs/doxygen/Doxyfile-version: ../docs/doxygen/Doxyfile-version.cmakein
docs/doxygen/Doxyfile-version: VersionInfo.cmake
docs/doxygen/Doxyfile-version: ../cmake/gmxConfigureVersionInfo.cmake
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating Doxyfile-version"
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -D VERSION_VARIABLES=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/VersionInfo.cmake -D VERSION_CMAKEIN=/work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/docs/doxygen/Doxyfile-version.cmakein -D VERSION_OUT=Doxyfile-version -P /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/cmake/gmxConfigureVersionInfo.cmake
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && /apps/cmake/3.16.4/bin/cmake -E touch Doxyfile-version

dep-graphs-dot: docs/doxygen/CMakeFiles/dep-graphs-dot
dep-graphs-dot: docs/doxygen/dep-graphs-dot-timestamp.txt
dep-graphs-dot: docs/doxygen/doxygen-xml-timestamp.txt
dep-graphs-dot: docs/doxygen/Doxyfile-version
dep-graphs-dot: docs/doxygen/CMakeFiles/dep-graphs-dot.dir/build.make

.PHONY : dep-graphs-dot

# Rule to build all files generated by this target.
docs/doxygen/CMakeFiles/dep-graphs-dot.dir/build: dep-graphs-dot

.PHONY : docs/doxygen/CMakeFiles/dep-graphs-dot.dir/build

docs/doxygen/CMakeFiles/dep-graphs-dot.dir/clean:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen && $(CMAKE_COMMAND) -P CMakeFiles/dep-graphs-dot.dir/cmake_clean.cmake
.PHONY : docs/doxygen/CMakeFiles/dep-graphs-dot.dir/clean

docs/doxygen/CMakeFiles/dep-graphs-dot.dir/depend:
	cd /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++ /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/docs/doxygen /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen /work_bgfs/g/gdayhoff/grayson/lab/gmx+/gromacs-2018.8++/build/docs/doxygen/CMakeFiles/dep-graphs-dot.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : docs/doxygen/CMakeFiles/dep-graphs-dot.dir/depend

