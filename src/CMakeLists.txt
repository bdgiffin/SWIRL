set(SWIRL_dependencies
    PARENT_SCOPE
    )

set(SWIRL_sources
    ${CMAKE_CURRENT_LIST_DIR}/SWIRL_API.cpp
    )

set(SWIRL_headers
    ${CMAKE_CURRENT_LIST_DIR}/Logger.h
    ${CMAKE_CURRENT_LIST_DIR}/Parameters.h
    ${CMAKE_CURRENT_LIST_DIR}/SWIRL.h
    ${CMAKE_CURRENT_LIST_DIR}/Particles.h
    ${CMAKE_CURRENT_LIST_DIR}/SpatialHash.h
    ${CMAKE_CURRENT_LIST_DIR}/Structure.h
    ${CMAKE_CURRENT_LIST_DIR}/WindField.h
    PARENT_SCOPE
    )
  
blt_add_library(
    NAME        SWIRL
    SOURCES     ${SWIRL_sources}
    HEADERS     ${SWIRL_headers}
    DEPENDS_ON  ${SWIRL_dependencies}
    SHARED      TRUE
    )
