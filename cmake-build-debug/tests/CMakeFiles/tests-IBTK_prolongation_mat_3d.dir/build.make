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
include tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/flags.make

tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o: tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/flags.make
tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o: ../tests/IBTK/prolongation_mat.cpp
tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o: tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o -MF CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o.d -o CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o -c /Users/qisun/sfw/ibamr/IBAMR/tests/IBTK/prolongation_mat.cpp

tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.i"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/qisun/sfw/ibamr/IBAMR/tests/IBTK/prolongation_mat.cpp > CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.i

tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.s"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && /usr/local/bin/mpicxx $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/qisun/sfw/ibamr/IBAMR/tests/IBTK/prolongation_mat.cpp -o CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.s

# Object files for target tests-IBTK_prolongation_mat_3d
tests__IBTK_prolongation_mat_3d_OBJECTS = \
"CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o"

# External object files for target tests-IBTK_prolongation_mat_3d
tests__IBTK_prolongation_mat_3d_EXTERNAL_OBJECTS =

tests/IBTK/prolongation_mat_3d: tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/IBTK/prolongation_mat.cpp.o
tests/IBTK/prolongation_mat_3d: tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/build.make
tests/IBTK/prolongation_mat_3d: lib/libIBAMR3d.dylib
tests/IBTK/prolongation_mat_3d: lib/libIBTK3d.dylib
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_algs.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_appu.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_geom.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_hier.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_math_std.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_mesh.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_pdat_std.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_solv.a
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/samrai/2.4.4/mac-g++-debug/lib/libSAMRAI3d_xfer.a
tests/IBTK/prolongation_mat_3d: /usr/local/Cellar/hdf5/1.12.2_2/lib/libhdf5.dylib
tests/IBTK/prolongation_mat_3d: /usr/local/lib/libsz.dylib
tests/IBTK/prolongation_mat_3d: /Users/qisun/anaconda3/lib/libz.dylib
tests/IBTK/prolongation_mat_3d: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib/libdl.tbd
tests/IBTK/prolongation_mat_3d: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/lib/libm.tbd
tests/IBTK/prolongation_mat_3d: lib/libBUNDLED_MUPARSER.dylib
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/petsc/3.13.4/mac-debug/lib/libpetsc.dylib
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/libmesh/1.6.2/1.6.2-debug/lib/libmesh_dbg.dylib
tests/IBTK/prolongation_mat_3d: /Users/qisun/sfw/libmesh/1.6.2/1.6.2-debug/lib/libmesh_dbg.dylib
tests/IBTK/prolongation_mat_3d: tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable IBTK/prolongation_mat_3d"
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/build: tests/IBTK/prolongation_mat_3d
.PHONY : tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/build

tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/clean:
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && $(CMAKE_COMMAND) -P CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/clean

tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/depend:
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/qisun/sfw/ibamr/IBAMR /Users/qisun/sfw/ibamr/IBAMR/tests /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/tests-IBTK_prolongation_mat_3d.dir/depend

