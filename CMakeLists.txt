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

find_package(tools REQUIRED)
include_directories(${tools_INCLUDE_DIRS})

####################################################################################################
## Static libraries
####################################################################################################

add_library(apt-core STATIC
    "src/core/crossover.hpp"
    "src/core/crossconfig.h"
    "src/core/crossconfig.cpp"
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
        "src/visu/phylogenyviewer.h"
        "src/visu/ptreeintrospecter.hpp"
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
    set(MAYBE_GUI "apt-gui;apt-basicviewer")
    set(GUI_LIB "$<TARGET_FILE_NAME:apt-gui>")

    add_executable(
        apt-basicviewer
        src/tests/basicviewer.cpp
    )
    target_link_libraries(apt-basicviewer apt-gui apt-core ${tools_LIBRARIES} ${QT_LIBS})
    set_property(TARGET apt-basicviewer PROPERTY POSITION_INDEPENDENT_CODE ON)

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

install(TARGETS apt-core ${MAYBE_GUI}
    ARCHIVE DESTINATION lib/kgd/apt
    RUNTIME DESTINATION bin/kgd/apt
)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/src/"
        DESTINATION include/kgd/apt
        FILES_MATCHING
            PATTERN "*.h"
        PERMISSIONS
            OWNER_EXECUTE OWNER_READ
            GROUP_EXECUTE GROUP_READ
            PATTERN "*.hpp"
        PERMISSIONS
            OWNER_EXECUTE OWNER_READ
            GROUP_EXECUTE GROUP_READ
)

set(CONFIG "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake")
file(GENERATE
    OUTPUT ${CONFIG}
    CONTENT
"# CMake configuration settings for project ${PROJECT_NAME} (${CMAKE_BUILD_TYPE})

set(${PROJECT_NAME}_INCLUDE_DIRS
    ${tools_INCLUDE_DIRS}
    ${CMAKE_INSTALL_PREFIX}/include)

set(${PROJECT_NAME}_LIBRARIES
    $<JOIN:${tools_LIBRARIES};$<TARGET_FILE_NAME:apt-core>;${GUI_LIB},
    ${CMAKE_INSTALL_PREFIX}/lib/kgd/apt/>)\n"
)
install(FILES ${CONFIG} DESTINATION cmake)
