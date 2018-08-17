cmake_minimum_required(VERSION 2.8)

####################################################################################################
## General settings
####################################################################################################

project("apt")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")#"-ftime-report")

if (${CMAKE_VERSION} VERSION_LESS 3.1)
    set(CMAKE_CXX_FLAGS "-std=c++17 ${CMAKE_CXX_FLAGS}")
else ()
    set(CMAKE_CXX_STANDARD 17)
endif ()

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")

####################################################################################################
## Self dependancies
####################################################################################################

set(CMAKE_PREFIX_PATH "../packages/")
find_package(tools REQUIRED)
message(STATUS "Tools found at " ${tools_DIR})
message(STATUS "Tools includes at " ${tools_INCLUDE_DIRS})
message(STATUS "Tools libs are " ${tools_LIBRARIES})
include_directories(${tools_INCLUDE_DIRS})

####################################################################################################
## Static libraries
####################################################################################################

add_library(apt-core STATIC
    "src/core/crossover.hpp"
    "src/core/phylogenictree.hpp"
    "src/core/ptreeconfig.h"
    "src/core/ptreeconfig.cpp"
)
target_link_libraries(apt-core ${tools_LIBRARIES})

####################################################################################################
## GUI management
####################################################################################################
if (NOT CLUSTER_BUILD)

    ################################################################################################
    ## Qt library
    ################################################################################################

    #set(CMAKE_PREFIX_PATH $(find $HOME/Qt* -path "*/lib/cmake"))

    find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui)
    message(STATUS "Qt found at " ${Qt5_DIR})
    set(QT_LIBS Qt5::Core Qt5::Widgets)

    ################################################################################################
    ## GUI sources files
    ################################################################################################

    add_library(apt-gui STATIC
        "src/core/crossover.hpp"
        "src/core/phylogenictree.hpp"
        "src/visu/phylogenyviewer.h"
        "src/visu/phylogenyviewer.cpp"
        "src/visu/graphicsviewzoom.h"
        "src/visu/graphicsviewzoom.cpp"
        "src/visu/graphicutils.h"
        "src/visu/graphicutils.cpp"
        "src/visu/standaloneviewer.hpp"
    )
    target_link_libraries(apt-gui apt-core ${tools_LIBRARIES} ${QT_LIBS})
    set_target_properties(apt-gui PROPERTIES AUTOMOC ON)
    set_property(TARGET apt-gui PROPERTY POSITION_INDEPENDENT_CODE ON)
    set(GUI_LIB "$<TARGET_FILE:apt-gui>")

endif()

####################################################################################################
## Options
####################################################################################################

message(STATUS "")

option(CLUSTER_BUILD "Whether or not building on a cluster (i-e no gui and local linkage)" ON)
message(STATUS "Cluster building mode is " ${CLUSTER_BUILD})
if(CLUSTER_BUILD)
    add_definitions(-DCLUSTER_BUILD)
endif()

####################################################################################################
## Package info generation
####################################################################################################

set(CONFIG "${PROJECT_SOURCE_DIR}/../packages/${PROJECT_NAME}-config.cmake")
file(GENERATE
    OUTPUT ${CONFIG}
    CONTENT
"# CMake configuration settings for project ${PROJECT_NAME}

set(${PROJECT_NAME}_INCLUDE_DIRS
    ${tools_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/src)

set(${PROJECT_NAME}_LIBRARIES
    $<JOIN:${tools_LIBRARIES};$<TARGET_FILE:apt-core>;${GUI_LIB},\n    >)\n"
)
