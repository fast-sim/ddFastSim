#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.15)
project(DDFastSim)

# Set up C++ Standard
# ``-DCMAKE_CXX_STANDARD=<standard>`` when invoking CMake
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

# Prevent CMake falls back to the latest standard the compiler does support
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Disables the use of compiler-specific extensions, hence makes sure the code
# works for a wider range of compilers
set(CMAKE_CXX_EXTENSIONS OFF)

#---------------------------------------------------
# Find Geant4 package

find_package(Geant4 REQUIRED)

#--------------------------------------------------
# Find DD4hep package

find_package(DD4hep REQUIRED COMPONENTS DDG4 DDParsers DDCore DDDetectors)

#--------------------------------------------------
set(DD4HEP_SET_RPATH ON)
dd4hep_set_compiler_flags()
dd4hep_use_python_executable()

file(GLOB sources src/*.cc)

set(headers  ${PROJECT_SOURCE_DIR}/include/ddfastsim/*.h)

if(DD4HEP_USE_PYROOT)
  message("dd4hep use pyroot")
  ROOT_GENERATE_DICTIONARY(G__DDFastSim ${headers} LINKDEF include/ROOT/LinkDef.h)
  list(APPEND sources G__DDFastSim.cxx)
endif()

add_dd4hep_plugin(${PROJECT_NAME} SHARED ${sources})

target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR}/include  )

target_link_libraries(${PROJECT_NAME} PUBLIC  ${Geant4_LIBRARIES} DD4hep::DDCore DD4hep::DDG4 )

install(TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}Targets)

dd4hep_instantiate_package(${PROJECT_NAME})

add_library(${PROJECT_NAME}detector SHARED Detector/src/SimpleCylinder_geo.cpp Detector/src/SandwichCylinders_geo.cpp)
target_link_libraries(${PROJECT_NAME}detector PUBLIC DD4hep::DDCore DD4hep::DDRec)
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(${PROJECT_NAME}detector PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
dd4hep_set_version(${PROJECT_NAME}detector MAJOR 1 MINOR 0 PATCH 0)
dd4hep_generate_rootmap(${PROJECT_NAME}detector)
install(TARGETS ${PROJECT_NAME}detector)

install(DIRECTORY Detector/xml DESTINATION Detector)
install(DIRECTORY options DESTINATION .)
