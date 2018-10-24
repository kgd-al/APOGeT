cmake_minimum_required(VERSION 3.4)

################################################################################
## General settings
################################################################################


project("apt")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")#"-ftime-report")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(LIB_BASE "${CMAKE_INSTALL_PREFIX}/lib/kgd/apt")


################################################################################
## Self dependancies
################################################################################

find_package(tools REQUIRED)
include_directories(${tools_INCLUDE_DIRS})


####################################################################################################
## Managing uneven support of std 17 filesystem
####################################################################################################

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  link_directories(/usr/local/opt/llvm/lib/)
  message(STATUS "> Linking against manual clang installation")
endif()


################################################################################
## Static library
################################################################################

add_library(apt-core STATIC
    "src/core/crossover.h"
    "src/core/crossover.cpp"
    "src/core/ptreeconfig.h"
    "src/core/ptreeconfig.cpp"

    "src/core/tree/enumvector.hpp"
    "src/core/tree/treetypes.h"
    "src/core/tree/treetypes.cpp"
    "src/core/tree/enveloppecriteria.cpp"
    "src/core/tree/callbacks.hpp"
    "src/core/tree/speciesdata.hpp"
    "src/core/tree/speciescontributors.cpp"
    "src/core/tree/speciescontributors.h"
    "src/core/tree/node.hpp"
    "src/core/tree/phylogenictree.hpp"
)
target_link_libraries(apt-core ${tools_LIBRARIES})

################################################################################
## GUI management
################################################################################

if (NOT CLUSTER_BUILD)

    ############################################################################
    ## Qt library
    ############################################################################

    #set(CMAKE_PREFIX_PATH $(find $HOME/Qt* -path "*/lib/cmake"))

    find_package(Qt5 REQUIRED COMPONENTS Core Widgets Gui)
    message(STATUS "Qt found at " ${Qt5_DIR})
    set(QT_LIBS Qt5::Core Qt5::Widgets)

    ############################################################################
    ## GUI sources files
    ############################################################################

    add_library(apt-gui STATIC
        "src/visu/phylogenyviewer.h"
        "src/visu/ptgraphbuilder.h"
        "src/visu/ptgraphbuilder.cpp"
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
    set(GUI_LIB "${LIB_BASE}/$<TARGET_FILE_NAME:apt-gui>")

    add_executable(
        apt-basicviewer
        src/tests/basicviewer.cpp
    )
    target_link_libraries(apt-basicviewer apt-gui apt-core ${tools_LIBRARIES} ${QT_LIBS})
    set_property(TARGET apt-basicviewer PROPERTY POSITION_INDEPENDENT_CODE ON)

endif()

################################################################################
## Options
################################################################################


message(STATUS "")

option(CLUSTER_BUILD "Whether or not building on a cluster (i-e no gui and local linkage)" ON)
message(STATUS "Cluster building mode is " ${CLUSTER_BUILD})
if(CLUSTER_BUILD)
    add_definitions(-DCLUSTER_BUILD)
endif()

option(WITH_DEBUG_INFO "Sets whether to maximize debug data collection" OFF)
message(STATUS "With debug info " ${WITH_DEBUG_INFO})
if(WITH_DEBUG_INFO)
    add_definitions(-DWITH_DEBUG_INFO)
    string(APPEND CMAKE_CXX_FLAGS " -g") # Enable debug info

    # Enable small memory error detector (fuse=gold fixes linker errors)
#    set(ASAN "-fsanitize=thread")
    set(ASAN "-fsanitize=address -fsanitize=undefined")
    string(APPEND CMAKE_CXX_FLAGS " ${ASAN} -fno-omit-frame-pointer")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        string(APPEND CMAKE_CXX_FLAGS " -fuse-ld=gold")
    endif()
endif()

option(BUILD_TESTS "Sets whether to build the tests executables" OFF)
message(STATUS "Build tests " ${BUILD_TESTS})

################################################################################
## Make documentation
################################################################################

option(GENERATE_HELP "Whether or not to build the help through doxygen" OFF)
message(STATUS "Help generation is " ${GENERATE_HELP})

if (GENERATE_HELP)
    cmake_policy(SET CMP0057 NEW)
    find_package(Doxygen REQUIRED dot)

    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/docs")

    # Doxygen needs to see the included files as well (e.g. macros)
    set(DOXYGEN_INCLUDE_PATH "${tools_INCLUDE_DIRS}")

    # Do expand macros for I may have used some...
    set(DOXYGEN_MACRO_EXPANSION YES)

    # Ignore generated moc files
    set(DOXYGEN_EXCLUDE_PATTERNS "*moc_*"
                                 "*.cpp")

    # Dot not look into the abyss
    set(DOXYGEN_EXCLUDE_SYMBOLS "_details::deepmagic::*"
                                "_*_link_checker"
                                "_*_metadata"
                                "CFILE")

    # Qt help
    set(DOXYGEN_GENERATE_QHP YES)
    set(DOXYGEN_QCH_FILE "../apt.qch")
    set(DOXYGEN_QHP_NAMESPACE "kgd.apt")
    set(DOXYGEN_QHP_VIRTUAL_FOLDER "kgd_docs")
    set(DOXYGEN_QHG_LOCATION "qhelpgenerator")

    # Debugging
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        message(STATUS "Enabling full documentation")
        set(DOXYGEN_ENABLED_SECTIONS "internal")
        set(DOXYGEN_EXTRACT_PRIVATE YES)
    endif()

    doxygen_add_docs(docs ${SOURCES} COMMENT "Generating documentation")
endif()


################################################################################
## Package info generation / installation
################################################################################

install(TARGETS apt-core ${MAYBE_GUI}
    ARCHIVE DESTINATION lib/kgd/apt
    RUNTIME DESTINATION bin/kgd/apt
)

install(
    DIRECTORY "${CMAKE_SOURCE_DIR}/src/"
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

if (GENERATE_HELP)
    install(
        DIRECTORY "${CMAKE_SOURCE_DIR}/docs/"
        DESTINATION docs/apt
    )
endif()

set(CONFIG "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.cmake")
file(GENERATE
    OUTPUT ${CONFIG}
    CONTENT
"# CMake configuration settings for project ${PROJECT_NAME} (${CMAKE_BUILD_TYPE})

set(${PROJECT_NAME}_INCLUDE_DIRS
    ${tools_INCLUDE_DIRS}
    ${CMAKE_INSTALL_PREFIX}/include)

set(${PROJECT_NAME}_LIBRARIES
    $<JOIN:${tools_LIBRARIES};${LIB_BASE}/$<TARGET_FILE_NAME:apt-core>;${GUI_LIB},
    >)\n"
)
install(FILES ${CONFIG} DESTINATION cmake)
