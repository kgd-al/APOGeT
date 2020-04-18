cmake_minimum_required(VERSION 3.4)

################################################################################
## General settings
################################################################################


project("APOGeT")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")#"-ftime-report")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")


################################################################################
## Multi-configuration installation path
################################################################################

if (${CMAKE_INSTALL_PREFIX} MATCHES "^$ENV{HOME}")
    string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
    string(APPEND CMAKE_INSTALL_PREFIX "/${BUILD_TYPE}")
    set(CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})
endif()
set(LIB_BASE "${CMAKE_INSTALL_PREFIX}/lib/kgd/apt")
message("Install path is ${CMAKE_INSTALL_PREFIX}/")


################################################################################
## Self dependancies
################################################################################

find_package(tools CONFIG REQUIRED)
include_directories(${tools_INCLUDE_DIRS})


####################################################################################################
## Managing uneven support of std 17 filesystem
####################################################################################################

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  link_directories(/usr/local/opt/llvm/lib/)
  message("> Linking against manual clang installation")
endif()


################################################################################
## Static library
################################################################################

FUNCTION(PREPEND output prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${output} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

set(CORE_SRC
    "crossover.h"
    "crossover.cpp"
    "ptreeconfig.h"
    "ptreeconfig.cpp"
)
PREPEND(CORE_SRC "src/core" ${CORE_SRC})

set(TREE_SRC
    "enumvector.hpp"
    "treetypes.h"
    "treetypes.cpp"
    "enveloppecriteria.cpp"
    "callbacks.hpp"
    "speciesdata.hpp"
    "speciescontributors.cpp"
    "speciescontributors.h"
    "node.hpp"
    "phylogenetictree.hpp"
)
PREPEND(TREE_SRC "src/core/tree" ${TREE_SRC})

add_library(apt-core STATIC ${CORE_SRC} ${TREE_SRC})
target_link_libraries(apt-core ${tools_LIBRARIES})

message("${CORE_SRC}")

################################################################################
## GUI management
################################################################################

if (NOT CLUSTER_BUILD)

    ############################################################################
    ## Qt library
    ############################################################################

    #set(CMAKE_PREFIX_PATH $(find $HOME/Qt* -path "*/lib/cmake"))

    set(QT_COMPONENTS Core Widgets Gui)
    set(QT_LIBS Qt5::Core Qt5::Widgets)

    if (NOT S{NO_PRINTER})
        set(QT_COMPONENTS ${QT_COMPONENTS} PrintSupport)
        set(QT_LIBS ${QT_LIBS} Qt5::PrintSupport)
    endif()

    if (NOT S{NO_SVG})
        set(QT_COMPONENTS ${QT_COMPONENTS} Svg)
        set(QT_LIBS ${QT_LIBS} Qt5::Svg)
    endif()

    find_package(Qt5 REQUIRED COMPONENTS ${QT_COMPONENTS})
    message("Qt found at " ${Qt5_DIR})


    ############################################################################
    ## GUI sources files
    ############################################################################

    set(VISU_SRC
        "phylogenyviewer.h"
        "layer.hpp"
        "ptgraphbuilder.h"
        "ptgraphbuilder.cpp"
        "phylogenyviewer.cpp"
        "speciestracking.h"
        "speciestracking.cpp"
        "graphicsviewzoom.h"
        "graphicsviewzoom.cpp"
        "graphicutils.h"
        "graphicutils.cpp"
        "pviewerconfig.h"
        "pviewerconfig.cpp"

        "standaloneviewer.hpp"
    )
    PREPEND(VISU_SRC "src/visu" ${VISU_SRC})
    add_library(apt-gui STATIC ${VISU_SRC})
    target_link_libraries(apt-gui apt-core ${tools_LIBRARIES} ${QT_LIBS})
    set_target_properties(apt-gui PROPERTIES AUTOMOC ON)
    set_property(TARGET apt-gui PROPERTY POSITION_INDEPENDENT_CODE ON)
    set(MAYBE_GUI "apt-gui")
    set(GUI_LIB "${LIB_BASE}/$<TARGET_FILE_NAME:apt-gui> ${QT_LIBS}")

    option(BUILD_GENERIC_VIEWER "Whether or not to build the generic viewer" ON)
    message("Building generic viewer " ${BUILD_GENERIC_VIEWER})
    if (BUILD_GENERIC_VIEWER)
        add_executable(
            apt-basicviewer
            src/tests/basicviewer.cpp
        )
        target_link_libraries(apt-basicviewer apt-gui apt-core ${tools_LIBRARIES} ${QT_LIBS})
        set_property(TARGET apt-basicviewer PROPERTY POSITION_INDEPENDENT_CODE ON)
    endif()

endif()

################################################################################
## Options
################################################################################


message("")

option(CLUSTER_BUILD "Whether or not building on a cluster (i-e no gui and local linkage)" ON)
message("Cluster building mode is " ${CLUSTER_BUILD})
if(CLUSTER_BUILD)
    add_definitions(-DCLUSTER_BUILD)
endif()

option(WITH_DEBUG_INFO "Sets whether to maximize debug data collection" OFF)
message("With debug info " ${WITH_DEBUG_INFO})
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
message("Build tests " ${BUILD_TESTS})

option(NO_PRINTER "Sets whether to disable QPrinter related capabilities" OFF)
message("No printer " ${NO_PRINTER})
if (NO_PRINTER)
    add_definitions(-NO_PRINTER)
endif()

option(NO_DVG "Sets whether to disable SVG related capabilities" OFF)
message("No printer " ${NO_SVG})
if (NO_SVG)
    add_definitions(-NO_SVG)
endif()

################################################################################
## Make documentation
################################################################################

option(GENERATE_HELP "Whether or not to build the help through doxygen" OFF)
message("Help generation is " ${GENERATE_HELP})

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
    set(DOXYGEN_EXCLUDE_SYMBOLS "_details*"
                                "SAGFieldMetadata"
                                "CFILE")

    # Qt help
    set(DOXYGEN_GENERATE_QHP YES)
    set(DOXYGEN_QCH_FILE "../${PROJECT_NAME}.qch")
    set(DOXYGEN_QHP_NAMESPACE "kgd.apt")
    set(DOXYGEN_QHP_VIRTUAL_FOLDER "kgd_docs")
    set(DOXYGEN_QHG_LOCATION "qhelpgenerator")

    # Debugging
    if (CMAKE_BUILD_TYPE MATCHES Debug)
        message("Enabling full documentation")
        set(DOXYGEN_ENABLED_SECTIONS "internal")
        set(DOXYGEN_EXTRACT_PRIVATE YES)
    endif()

    doxygen_add_docs(docs ${SOURCES} COMMENT "Generating documentation")
endif()


################################################################################
## Package info generation / installation
################################################################################

install(TARGETS apt-core ARCHIVE DESTINATION lib/kgd/apt)
if (NOT CLUSTER_BUILD)
    install(TARGETS apt-gui ARCHIVE DESTINATION lib/kgd/apt)
endif()

if (BUILD_GENERIC_VIEWER AND NOT CLUSTER_BUILD)
    install(TARGETS apt-basicviewer
        RUNTIME DESTINATION bin/kgd/apt)
endif()

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

set(CONFIG "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake")
file(GENERATE
    OUTPUT ${CONFIG}
    CONTENT
"# CMake configuration settings for project ${PROJECT_NAME} (${CMAKE_BUILD_TYPE})

set(${PROJECT_NAME}_INCLUDE_DIRS
    ${tools_INCLUDE_DIRS}
    ${CMAKE_INSTALL_PREFIX}/include)

set(${PROJECT_NAME}_LIBRARIES
    $<JOIN:${tools_LIBRARIES};${LIB_BASE}/$<TARGET_FILE_NAME:apt-core>,
    >)

if(NOT CLUSTER_BUILD)
    set(${PROJECT_NAME}_GUI_LIBRARIES $<JOIN:${GUI_LIB},
        >)

    find_package(Qt5 REQUIRED COMPONENTS ${QT_COMPONENTS})
    message(\"Qt found at \" ${Qt5_DIR})
    set(QT_LIBS ${QT_LIBS})
endif()

set(${PROJECT_NAME}_BUILD_TYPE ${CMAKE_BUILD_TYPE})
if(NOT \${${PROJECT_NAME}_BUILD_TYPE} STREQUAL \${CMAKE_BUILD_TYPE})
    message(FATAL_ERROR \"Mismatching build version!
        Linking against '\${${PROJECT_NAME}_BUILD_TYPE}' while using '\${CMAKE_BUILD_TYPE}'\")
endif()\n"
)
install(FILES ${CONFIG} DESTINATION cmake)

set(DEPENDENCY_GENERATOR ${tools_DIR}/generate_dependency_data.sh)
set(DEPENDENCY_FILE "${CMAKE_CURRENT_BINARY_DIR}/apt.dependency")
add_custom_target(dependencies ALL DEPENDS ${DEPENDENCY_FILE})
add_custom_command(
    OUTPUT ${DEPENDENCY_FILE} always_check
    COMMAND ${DEPENDENCY_GENERATOR}
        ${DEPENDENCY_FILE} APOGeT ${CMAKE_BUILD_TYPE}

    COMMENT "Generating dependency data"
)
install(FILES ${DEPENDENCY_FILE} DESTINATION cmake)
