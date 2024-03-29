# Setup Project for CMAKE and find Components
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(reactorBay_sourceFiles)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})

# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
if(ROOT_FOUND)
add_definitions(-DG4ANALYSIS_USE_ROOT)
include(${Geant4_USE_FILE})
include_directories(${Geant4_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include)
else()
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
endif()

# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add the executable, and link it to the Geant4 libraries
#
add_executable(reactorBay main.cc ${sources} ${headers})
if(ROOT_FOUND)
target_link_libraries(reactorBay ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
else()
target_link_libraries(reactorBay ${Geant4_LIBRARIES})
endif()

# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(reactorBay_sourceFiles DEPENDS reactorBay)

# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS reactorBay DESTINATION bin)