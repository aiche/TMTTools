### example CMakeLists.txt to develop programs using OpenMS
project("TMTTools")
cmake_minimum_required(VERSION 2.8)

## list all your executables here (a corresponding .C file should exist, e.g. Main.C)
set(my_executables
	TMTExporter
)

## find OpenMS configuration and register target "OpenMS" (our library)
find_package(OpenMS)

# check whether the OpenMS package was found
if (OpenMS_FOUND)

  ## include directories for OpenMS headers (and contrib)
  include_directories(${OPENMS_INCLUDE_DIRS})

  ## append precompiler macros and compiler flags specific to OpenMS
  ## Warning: this could be harmful to your project. Check this if problems occur.
  ## Also, use this to add your own compiler flags, e.g. for OpenMP support.
  ## e.g. for Visual Studio use /openmp
  ## set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS} /openmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENMS_ADDCXX_FLAGS}")
  add_definitions(${OPENMS_DEFINITIONS})


  ## add targets for the executables
  foreach(i ${my_executables})
    add_executable(${i} ${i}.C)
    ## link executables against OpenMS
  	target_link_libraries(${i} OpenMS)
  endforeach(i)

else(OpenMS_FOUND)
  message(FATAL_ERROR "OpenMSConfig.cmake file not found!")
endif(OpenMS_FOUND)
