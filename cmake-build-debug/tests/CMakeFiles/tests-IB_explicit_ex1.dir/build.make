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
include tests/CMakeFiles/tests-IB_explicit_ex1.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/tests-IB_explicit_ex1.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/tests-IB_explicit_ex1.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/tests-IB_explicit_ex1.dir/flags.make

tests/CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o: tests/CMakeFiles/tests-IB_explicit_ex1.dir/flags.make
tests/CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o: ../tests/IB/explicit_ex1.cpp
tests/CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o: tests/CMakeFiles/tests-IB_explicit_ex1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests/CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o -MF CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o.d -o CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o -c /Users/qisun/sfw/ibamr/IBAMR/tests/IB/explicit_ex1.cpp

tests/CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.i"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/qisun/sfw/ibamr/IBAMR/tests/IB/explicit_ex1.cpp > CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.i

tests/CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.s"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/qisun/sfw/ibamr/IBAMR/tests/IB/explicit_ex1.cpp -o CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.s

# Object files for target tests-IB_explicit_ex1
tests__IB_explicit_ex1_OBJECTS = \
"CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o"

# External object files for target tests-IB_explicit_ex1
tests__IB_explicit_ex1_EXTERNAL_OBJECTS =

tests/IB/explicit_ex1: tests/CMakeFiles/tests-IB_explicit_ex1.dir/IB/explicit_ex1.cpp.o
tests/IB/explicit_ex1: tests/CMakeFiles/tests-IB_explicit_ex1.dir/build.make
tests/IB/explicit_ex1: lib/libIBAMR2d.dylib
tests/IB/explicit_ex1: lib/libIBTK2d.dylib
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_algs.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_appu.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_geom.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_hier.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_math_std.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_mesh.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_pdat_std.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_solv.a
tests/IB/explicit_ex1: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI2d_xfer.a
tests/IB/explicit_ex1: /usr/local/Cellar/hdf5/1.12.2_2/lib/libhdf5.dylib
tests/IB/explicit_ex1: /usr/local/lib/libsz.dylib
tests/IB/explicit_ex1: /Users/qisun/anaconda3/lib/libz.dylib
tests/IB/explicit_ex1: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib/libdl.tbd
tests/IB/explicit_ex1: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib/libm.tbd
tests/IB/explicit_ex1: lib/libBUNDLED_MUPARSER.dylib
tests/IB/explicit_ex1: /Users/qisun/sfw/petsc/3.13.4/mac-debug/lib/libpetsc.dylib
tests/IB/explicit_ex1: /Users/qisun/sfw/libmesh/1.6.2/1.6.2-debug/lib/libmesh_dbg.dylib
tests/IB/explicit_ex1: /Users/qisun/sfw/libmesh/1.6.2/1.6.2-debug/lib/libmesh_dbg.dylib
tests/IB/explicit_ex1: tests/CMakeFiles/tests-IB_explicit_ex1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable IB/explicit_ex1"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tests-IB_explicit_ex1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/tests-IB_explicit_ex1.dir/build: tests/IB/explicit_ex1
.PHONY : tests/CMakeFiles/tests-IB_explicit_ex1.dir/build

tests/CMakeFiles/tests-IB_explicit_ex1.dir/clean:
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && $(CMAKE_COMMAND) -P CMakeFiles/tests-IB_explicit_ex1.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/tests-IB_explicit_ex1.dir/clean

tests/CMakeFiles/tests-IB_explicit_ex1.dir/depend:
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/qisun/sfw/ibamr/IBAMR /Users/qisun/sfw/ibamr/IBAMR/tests /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests/CMakeFiles/tests-IB_explicit_ex1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/tests-IB_explicit_ex1.dir/depend

