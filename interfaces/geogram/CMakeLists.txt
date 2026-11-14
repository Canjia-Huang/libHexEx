#
# The geogram bridge: converts the trajectory layout into GEO::Mesh, with all
# information stored in geogram attributes, and writes it with geogram's mesh IO.
#
# Only built if geogram was found (see cmake/FindGeogram.cmake and the
# HEXEX_WITH_GEOGRAM option in the top level CMakeLists.txt).
#

add_library (HexExGeogram ${CMAKE_CURRENT_SOURCE_DIR}/hexex_geogram.cc)
add_library (HexEx::HexExGeogram ALIAS HexExGeogram)

target_include_directories (HexExGeogram
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${GEOGRAM_INCLUDE_DIRS}
    )

target_link_libraries (HexExGeogram
    PUBLIC
        HexEx
        ${GEOGRAM_LIBRARIES}
    )

if (APPLE AND GEOGRAM_LIBRARY_DIR)
    # geogram's libraries refer to each other through @rpath and are not
    # installed, so the runtime search path has to point at the build directory
    set_target_properties (HexExGeogram PROPERTIES
        BUILD_RPATH "${GEOGRAM_LIBRARY_DIR}"
        INSTALL_RPATH "${GEOGRAM_LIBRARY_DIR}")
endif ()

set_target_properties (HexExGeogram PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    AUTOMOC OFF
    AUTOUIC OFF
    )
