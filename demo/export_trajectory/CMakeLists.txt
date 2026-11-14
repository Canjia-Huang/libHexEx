#
# export_trajectory: extract a hex mesh and write the trajectories of all its
# elements as geogram meshes.
#
# Only built if geogram was found (see the HEXEX_WITH_GEOGRAM option).
#

add_executable (export_trajectory ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)

target_link_libraries (export_trajectory
    PRIVATE
        HexExGeogram
    )

target_include_directories (export_trajectory
    PRIVATE
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/src
    )

set_target_properties (export_trajectory PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build/bin")
