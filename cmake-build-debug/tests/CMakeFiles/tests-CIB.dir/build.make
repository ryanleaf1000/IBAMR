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

# Utility rule file for tests-CIB.

# Include any custom commands dependencies for this target.
include tests/CMakeFiles/tests-CIB.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/tests-CIB.dir/progress.make

tests-CIB: tests/CMakeFiles/tests-CIB.dir/build.make
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && bash /Users/qisun/sfw/ibamr/IBAMR/tests/link-test-files.sh /Users/qisun/sfw/ibamr/IBAMR/tests/CIB /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests/CIB
.PHONY : tests-CIB

# Rule to build all files generated by this target.
tests/CMakeFiles/tests-CIB.dir/build: tests-CIB
.PHONY : tests/CMakeFiles/tests-CIB.dir/build

tests/CMakeFiles/tests-CIB.dir/clean:
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests && $(CMAKE_COMMAND) -P CMakeFiles/tests-CIB.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/tests-CIB.dir/clean

tests/CMakeFiles/tests-CIB.dir/depend:
	cd /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/qisun/sfw/ibamr/IBAMR /Users/qisun/sfw/ibamr/IBAMR/tests /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests/CMakeFiles/tests-CIB.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/tests-CIB.dir/depend

