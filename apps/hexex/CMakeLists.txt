#
# hexex - the HexEx command line application.
#
# Reads a tetrahedral mesh with an integer grid map, extracts the hex mesh and
# writes the hex mesh, and - if geogram is available - the trajectories of all
# extracted elements as geogram meshes.
#
# It uses CLI11 (third_party/CLI11, a git submodule).
#

set (HEXEX_APP_DEFINITIONS "")

if (HEXEX_GEOGRAM_FOUND)
    # HexExGeogram links HexEx publicly, so HexEx comes in transitively
    set (HEXEX_APP_LIBRARIES CLI11::CLI11 HexExGeogram)
    list (APPEND HEXEX_APP_DEFINITIONS HEXEX_HAVE_GEOGRAM=1)
else ()
    set (HEXEX_APP_LIBRARIES CLI11::CLI11 HexEx)
endif ()

add_executable (hexex ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)

target_link_libraries (hexex PRIVATE ${HEXEX_APP_LIBRARIES})
target_compile_definitions (hexex PRIVATE ${HEXEX_APP_DEFINITIONS})

set_target_properties (hexex PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build/bin")

target_compile_options (hexex PRIVATE
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra>
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra>)
