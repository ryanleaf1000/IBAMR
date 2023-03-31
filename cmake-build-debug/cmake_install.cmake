# Install script for directory: /Users/qisun/sfw/ibamr/IBAMR

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/qisun/sfw/ibamr/0.12.0/0.12.0-debug")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Users/qisun/anaconda3/bin/llvm-objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/contrib/muparser" TYPE FILE MESSAGE_LAZY FILES
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParser.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserBase.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserBytecode.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserCallback.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserDLL.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserDef.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserError.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserFixes.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserInt.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserStack.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserTemplateMagic.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserTest.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserToken.h"
    "/Users/qisun/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserTokenReader.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY MESSAGE_LAZY FILES "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/lib/libBUNDLED_MUPARSER.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.dylib")
    execute_process(COMMAND /Users/qisun/anaconda3/bin/install_name_tool
      -delete_rpath "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/lib/"
      -add_rpath "/Users/qisun/sfw/ibamr/0.12.0/0.12.0-debug/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Users/qisun/anaconda3/bin/llvm-strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ibtk" TYPE FILE MESSAGE_LAZY FILES "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/ibtk/include/ibtk/config.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY MESSAGE_LAZY FILES "/Users/qisun/sfw/ibamr/IBAMR/include/ibamr")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY MESSAGE_LAZY FILES "/Users/qisun/sfw/ibamr/IBAMR/ibtk/include/ibtk")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr/IBAMRTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr/IBAMRTargets.cmake"
         "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles/Export/lib/cmake/ibamr/IBAMRTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr/IBAMRTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr/IBAMRTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr" TYPE FILE MESSAGE_LAZY FILES "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles/Export/lib/cmake/ibamr/IBAMRTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr" TYPE FILE MESSAGE_LAZY FILES "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles/Export/lib/cmake/ibamr/IBAMRTargets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr" TYPE FILE MESSAGE_LAZY FILES
    "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/cmake/IBAMRConfig.cmake"
    "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/cmake/IBAMRConfigVersion.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/ibtk/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/src/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/tests/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/examples/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/qisun/sfw/ibamr/IBAMR/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
