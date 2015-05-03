cmake_minimum_required(VERSION 2.6)

project(TA C CXX)# Fortran)

# add a NativeRelease build type
SET( CMAKE_CXX_FLAGS_NATIVERELEASE "-O3 -DNDEBUG -march=native -ftree-vectorize" CACHE STRING
    "Flags used by the C++ compiler during native builds."
    FORCE )
SET( CMAKE_C_FLAGS_NATIVERELEASE "-O3 -DNDEBUG -march=native -ftree-vectorize" CACHE STRING
    "Flags used by the C compiler during native builds."
    FORCE )
SET( CMAKE_EXE_LINKER_FLAGS_NATIVERELEASE "" CACHE STRING
    "Flags used for linking binaries during native builds."
    FORCE )
SET( CMAKE_SHARED_LINKER_FLAGS_NATIVERELEASE "" CACHE STRING
    "Flags used by the shared libraries linker during native builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_NATIVERELEASE
    CMAKE_C_FLAGS_NATIVERELEASE
    CMAKE_EXE_LINKER_FLAGS_NATIVERELEASE
    CMAKE_SHARED_LINKER_FLAGS_NATIVERELEASE )

if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES NativeRelease)
  list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES "${CMAKE_CONFIGURATION_TYPES}" CACHE STRING
     "Add the native configuration"
     FORCE)
endif()

## default build type is native
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE NativeRelease CACHE STRING
      "Choose the type of build. Options are: None Debug NativeRelease Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)
MESSAGE(STATUS "Build type: " ${CMAKE_BUILD_TYPE})




# handy stuff for eclipse
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# definitions for visual studio
add_definitions(-D_CRT_SECURE_NO_WARNINGS)

# my own includes/libs
# #include_directories(~/local/include)
#link_directories(~/local/lib)









# search for boost library
# SET(BOOST_ROOT ~/local) # default
# SET(Boost_NO_SYSTEM_PATHS TRUE) # needed with cmake on lenny to prevent it grabbing an old system version...
find_package(Boost  COMPONENTS serialization thread program_options regex filesystem chrono REQUIRED)

message(STATUS "Boost found " ${Boost_LIBRARY_DIRS})
link_directories(${Boost_LIBRARY_DIRS})
include_directories(${Boost_INCLUDE_DIRS})
message(STATUS "Boost_LIBRARY_DIRS in  " ${Boost_LIBRARY_DIRS})
message(STATUS "Boost_INCLUDE_DIRS in  " ${Boost_INCLUDE_DIRS})


# search for GSL
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${TA_SOURCE_DIR})
# find_package(GSL REQUIRED)
# include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})



find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

################################## for finding MLPACK ###########################



# find_package(Armadillo 3.6.0 REQUIRED)
# include_directories(${ARMADILLO_INCLUDE_DIRS})
# 
# find_package(LibXml2 2.6.0 REQUIRED)
# include_directories(${LIBXML2_INCLUDE_DIR})



################################## ################################## ################################## 

# include directories
include_directories(${TA_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} )

# install directories
SET(LIB_INSTALL_DIR lib)
SET(INCLUDE_INSTALL_DIR include)

# go
#configure_file("${TA_SOURCE_DIR}/log4j.properties" "${CMAKE_CURRENT_BINARY_DIR}/")
add_subdirectory(util)
add_subdirectory(taLib)
add_subdirectory(tests)
add_subdirectory(tools)





