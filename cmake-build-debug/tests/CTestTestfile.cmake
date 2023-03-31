# CMake generated Testfile for 
# Source directory: /Users/qisun/sfw/ibamr/IBAMR/tests
# Build directory: /Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(autotests "/Users/qisun/sfw/ibamr/IBAMR/attest" "-j2" "--verbose" "--test-timeout=120" "-E" "mpirun=[3-9]|explicit_ex2_3d.nodal_quadrature.input|explicit_ex5_3d.mpirun=2.input|explicit_ex8_2d.input|explicit_ex8_2d.scratch_hier.input|nwt_cylinder|cib_double_shell.input|cib_double_shell.cholesky.input")
set_tests_properties(autotests PROPERTIES  ENVIRONMENT "PYTHONUNBUFFERED=yes" TIMEOUT "1800" WORKING_DIRECTORY "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug" _BACKTRACE_TRIPLES "/Users/qisun/sfw/ibamr/IBAMR/tests/CMakeLists.txt;336;ADD_TEST;/Users/qisun/sfw/ibamr/IBAMR/tests/CMakeLists.txt;0;")
