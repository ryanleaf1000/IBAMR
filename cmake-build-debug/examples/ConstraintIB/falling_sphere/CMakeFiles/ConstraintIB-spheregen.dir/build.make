# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake

# The command to remove a file.
RM = /Applications/CLion.app/Contents/bin/cmake/mac/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/qisun/sfw/ibamr/IBAMR

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug

# Include any dependencies generated for this target.
include examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/progress.make

# Include the compile flags for this target's objects.
include examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/flags.make

examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o: examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/flags.make
examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o: ../examples/ConstraintIB/falling_sphere/sphereGen3d.cpp
examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o: examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/examples/ConstraintIB/falling_sphere && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o -MF CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o.d -o CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o -c /Users/qisun/sfw/ibamr/IBAMR/examples/ConstraintIB/falling_sphere/sphereGen3d.cpp

examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.i"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/examples/ConstraintIB/falling_sphere && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/qisun/sfw/ibamr/IBAMR/examples/ConstraintIB/falling_sphere/sphereGen3d.cpp > CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.i

examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.s"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/examples/ConstraintIB/falling_sphere && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/qisun/sfw/ibamr/IBAMR/examples/ConstraintIB/falling_sphere/sphereGen3d.cpp -o CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.s

# Object files for target ConstraintIB-spheregen
ConstraintIB__spheregen_OBJECTS = \
"CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o"

# External object files for target ConstraintIB-spheregen
ConstraintIB__spheregen_EXTERNAL_OBJECTS =

examples/ConstraintIB/falling_sphere/sphereGen3d: examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/sphereGen3d.cpp.o
examples/ConstraintIB/falling_sphere/sphereGen3d: examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/build.make
examples/ConstraintIB/falling_sphere/sphereGen3d: examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable sphereGen3d"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/examples/ConstraintIB/falling_sphere && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ConstraintIB-spheregen.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/build: examples/ConstraintIB/falling_sphere/sphereGen3d
.PHONY : examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/build

examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/clean:
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/examples/ConstraintIB/falling_sphere && $(CMAKE_COMMAND) -P CMakeFiles/ConstraintIB-spheregen.dir/cmake_clean.cmake
.PHONY : examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/clean

examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/depend:
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/qisun/sfw/ibamr/IBAMR /Users/qisun/sfw/ibamr/IBAMR/examples/ConstraintIB/falling_sphere /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/examples/ConstraintIB/falling_sphere /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/ConstraintIB/falling_sphere/CMakeFiles/ConstraintIB-spheregen.dir/depend

