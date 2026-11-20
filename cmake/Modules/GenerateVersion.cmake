#######################################################################
# generate git-version signature
#######################################################################
include_directories(${CMAKE_BINARY_DIR})
find_package(Git)

set(oofem_version_cmd
  ${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/src/oofem_version.h.in
                   -D DST=${CMAKE_BINARY_DIR}/oofem_version.h
                   -D GIT_EXECUTABLE=${GIT_EXECUTABLE}
                   -P ${CMAKE_SOURCE_DIR}/src/GenerateVersionHeader.cmake
  )

# Generate once at configure time. The object libraries (core, sm_elements, ...)
# all include oofem_version.h via oofemcfg.h, but only liboofem and oofem get an
# add_dependencies() edge to the target below, so a parallel build races it.
execute_process(COMMAND ${oofem_version_cmd})

# Refresh on every build so the recorded git hash stays current.
add_custom_target(version ${oofem_version_cmd})
