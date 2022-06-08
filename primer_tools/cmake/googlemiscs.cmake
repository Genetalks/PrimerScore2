if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # when using the makefile generator, use the special variable $(MAKE) to invoke make
    # this enables the jobserver to work correctly
    set(MAKE_COMMAND "$(MAKE)")
else()
    # invoke make explicitly
    # in this case, we assume the parent build system is running in parallel already so no -j flag is added
    find_program(MAKE_COMMAND NAMES make gmake)
endif()

# Build gflags as an external project.
set(THIRD_PARTY_DIR_NAME vendor)
set(NAME gflags)
set (LIB_NAME ${NAME}lib)
set(${NAME}_PREFIX ${PROJECT_BINARY_DIR}/${THIRD_PARTY_DIR_NAME}/${LIB_NAME})
set(${NAME}_SOURCE_DIR ${CMAKE_SOURCE_DIR}/${THIRD_PARTY_DIR_NAME}/${NAME})
set(${NAME}_INSTALL_DIR ${${NAME}_PREFIX}/${NAME})
set(${NAME}_I_DIR ${${NAME}_INSTALL_DIR}/include)
set(${NAME}_L_DIR ${${NAME}_INSTALL_DIR}/lib)

ExternalProject_Add(${LIB_NAME}
    PREFIX      ${${NAME}_PREFIX}
    SOURCE_DIR  ${${NAME}_SOURCE_DIR}
    INSTALL_DIR ${${NAME}_INSTALL_DIR}
    CMAKE_ARGS  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_CXX_FLAGS=${EXTERNAL_PROJECT_CMAKE_CXX_FLAGS}
	-DCMAKE_INSTALL_PREFIX:PATH=${${NAME}_INSTALL_DIR}
    BUILD_COMMAND ${MAKE_COMMAND} -j
)
set(GFLAGS_LIB ${LIB_NAME})
include_directories(BEFORE SYSTEM ${${NAME}_I_DIR})
link_directories(${${NAME}_L_DIR})

# Build glog as an external project.
set(NAME glog)
set (LIB_NAME ${NAME}lib)
set(${NAME}_PREFIX ${PROJECT_BINARY_DIR}/${THIRD_PARTY_DIR_NAME}/${LIB_NAME})
set(${NAME}_SOURCE_DIR ${CMAKE_SOURCE_DIR}/${THIRD_PARTY_DIR_NAME}/${NAME})
set(${NAME}_INSTALL_DIR ${${NAME}_PREFIX}/${NAME})
set(${NAME}_I_DIR ${${NAME}_INSTALL_DIR}/include)
set(${NAME}_L_DIR ${${NAME}_INSTALL_DIR}/lib)
ExternalProject_Add(${LIB_NAME}
    DEPENDS    ${GFLAGS_LIB}
    PREFIX      ${${NAME}_PREFIX}
    SOURCE_DIR  ${${NAME}_SOURCE_DIR}
    INSTALL_DIR ${${NAME}_INSTALL_DIR}
    CMAKE_ARGS  -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_CXX_FLAGS=${EXTERNAL_PROJECT_CMAKE_CXX_FLAGS}
        -Dcobalt_gflags_DIR=${gflags_INSTALL_DIR}
	-DCMAKE_INSTALL_PREFIX:PATH=${${NAME}_INSTALL_DIR}
        -DCMAKE_INSTALL_LIBDIR:PATH=lib                
    BUILD_COMMAND ${MAKE_COMMAND} -j
)

include_directories(BEFORE SYSTEM ${${NAME}_I_DIR})
link_directories(${${NAME}_L_DIR})

set(glog_LIBS ${glog_L_DIR}/libglog.a ${gflags_L_DIR}/libgflags.a)


