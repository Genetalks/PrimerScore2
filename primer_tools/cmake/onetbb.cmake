# set(TBB_DIR ../vendor/onetbb/cmake)
# include(${CMAKE_SOURCE_DIR}/vendor/oneTBB/cmake/TBBBuild.cmake)
# tbb_build(TBB_ROOT ${CMAKE_SOURCE_DIR}/vendor/oneTBB CONFIG_DIR TBB_DIR)

if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # when using the makefile generator, use the special variable $(MAKE) to invoke make
    # this enables the jobserver to work correctly
    set(MAKE_COMMAND "$(MAKE)")
else()
    # invoke make explicitly
    # in this case, we assume the parent build system is running in parallel already so no -j flag is added
    find_program(MAKE_COMMAND NAMES make gmake)
endif()


set(THIRD_PARTY_DIR_NAME vendor)
set(NAME onetbb)
set (LIB_NAME ${NAME}lib)
set(${NAME}_PREFIX ${PROJECT_BINARY_DIR}/${THIRD_PARTY_DIR_NAME}/${LIB_NAME})
set(${NAME}_SOURCE_DIR ${PROJECT_SOURCE_DIR}/${THIRD_PARTY_DIR_NAME}/oneTBB)
set(${NAME}_BINARY_DIR ${${NAME}_PREFIX})
set(${NAME}_I_DIR ${${NAME}_SOURCE_DIR}/include)
set(${NAME}_L_DIR ${${NAME}_BINARY_DIR})

# build oneTBB
ExternalProject_Add(
    ${LIB_NAME}
    PREFIX ${${NAME}_PREFIX}
    SOURCE_DIR  ${${NAME}_SOURCE_DIR}
    CONFIGURE_COMMAND ""
    #BUILD_COMMAND ${MAKE_COMMAND} -j arch=intel64 tbb_root=${${NAME}_SOURCE_DIR} extra_inc=linux.inc extra_inc=../../tbb_static.inc  tbb_build_prefix=tbb tbb_build_dir=${${NAME}_BINARY_DIR}
    BUILD_COMMAND ${MAKE_COMMAND} -j arch=intel64 tbb_root=${${NAME}_SOURCE_DIR} extra_inc=big_iron.inc tbb_build_prefix=tbb tbb_build_dir=${${NAME}_BINARY_DIR}
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND ""
    LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
)

include_directories(BEFORE_SYSTEM ${${NAME}_I_DIR})
set(tbb_LIBS ${${NAME}_L_DIR}/tbb_release/libtbb.a)
