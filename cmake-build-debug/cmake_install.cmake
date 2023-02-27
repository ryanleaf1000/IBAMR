# Install script for directory: /home/qi/sfw/ibamr/IBAMR

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
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

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/contrib/eigen" TYPE DIRECTORY MESSAGE_LAZY FILES "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/eigen/Eigen")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/contrib/eigen" TYPE DIRECTORY MESSAGE_LAZY FILES "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/eigen/unsupported")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/contrib/muparser" TYPE FILE MESSAGE_LAZY FILES
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParser.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserBase.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserBytecode.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserCallback.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserDLL.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserDef.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserError.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserFixes.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserInt.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserStack.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserTemplateMagic.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserTest.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserToken.h"
    "/home/qi/sfw/ibamr/IBAMR/ibtk/contrib/muparser/include/muParserTokenReader.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.so"
         RPATH "/usr/local/lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY MESSAGE_LAZY FILES "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/lib/libBUNDLED_MUPARSER.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.so"
         OLD_RPATH "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/lib/:"
         NEW_RPATH "/usr/local/lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libBUNDLED_MUPARSER.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibraryx" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ibtk" TYPE FILE MESSAGE_LAZY FILES "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/ibtk/include/ibtk/config.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY MESSAGE_LAZY FILES "/home/qi/sfw/ibamr/IBAMR/include/ibamr")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY MESSAGE_LAZY FILES "/home/qi/sfw/ibamr/IBAMR/ibtk/include/ibtk")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr/IBAMRTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr/IBAMRTargets.cmake"
         "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles/Export/lib/cmake/ibamr/IBAMRTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr/IBAMRTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr/IBAMRTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr" TYPE FILE MESSAGE_LAZY FILES "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles/Export/lib/cmake/ibamr/IBAMRTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr" TYPE FILE MESSAGE_LAZY FILES "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/CMakeFiles/Export/lib/cmake/ibamr/IBAMRTargets-debug.cmake")
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ibamr" TYPE FILE MESSAGE_LAZY FILES
    "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/cmake/IBAMRConfig.cmake"
    "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/cmake/IBAMRConfigVersion.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/ibtk/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/src/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/tests/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/examples/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/qi/sfw/ibamr/IBAMR/cmake-build-debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
