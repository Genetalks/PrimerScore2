if (CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    # when using the makefile generator, use the special variable $(MAKE) to invoke make
    # this enables the jobserver to work correctly
    set(MAKE_COMMAND "$(MAKE)")
else()
    # invoke make explicitly
    # in this case, we assume the parent build system is running in parallel already so no -j flag is added
    find_program(MAKE_COMMAND NAMES make gmake)
endif()


# build zlib
set(THIRD_PARTY_DIR_NAME vendor)
set(NAME zlib)
set(LIB_NAME zlib)
set(${NAME}_DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/${THIRD_PARTY_DIR_NAME})
set(${NAME}_URL http://zlib.net/fossils/zlib-1.2.11.tar.gz)
set(${NAME}_MD5 1c9f62f0778697a09d36121ead88e08e)
set(${NAME}_PREFIX ${PROJECT_BINARY_DIR}/${THIRD_PARTY_DIR_NAME}/${LIB_NAME})
set(${NAME}_SOURCE_DIR ${${NAME}_PREFIX}/zlib-1.2.11)
set(${NAME}_BINARY_DIR ${${NAME}_SOURCE_DIR})
set(${NAME}_INSTALL_DIR ${${NAME}_PREFIX}/${NAME})
set(${NAME}_I_DIR ${${NAME}_INSTALL_DIR}/include)
set(${NAME}_L_DIR ${${NAME}_INSTALL_DIR}/lib)

ExternalProject_Add(${LIB_NAME}
    PREFIX ${${NAME}_PREFIX}
    DOWNLOAD_DIR ${${NAME}_DOWNLOAD_DIR}
    URL          ${${NAME}_URL}
    URL_HASH     MD5=${${NAME}_MD5}
    SOURCE_DIR  ${${NAME}_SOURCE_DIR}
    CONFIGURE_COMMAND ${${NAME}_SOURCE_DIR}/configure --prefix=${${NAME}_INSTALL_DIR}
    BUILD_IN_SOURCE 1
    BUILD_COMMAND COMMAND ${MAKE_COMMAND} -j
    INSTALL_COMMAND ${MAKE_COMMAND} install
    LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
)
include_directories(BEFORE SYSTEM ${${NAME}_I_DIR})
set(${LIB_NAME}_LIB ${${NAME}_L_DIR}/libz.a)

# build deflate
set(NAME deflate)
set(LIB_NAME ${NAME}lib)
set(${NAME}_PREFIX ${PROJECT_BINARY_DIR}/${THIRD_PARTY_DIR_NAME}/${LIB_NAME})
set(${NAME}_SOURCE_DIR ${CMAKE_SOURCE_DIR}/${THIRD_PARTY_DIR_NAME}/libdeflate)
set(${NAME}_BINARY_DIR ${${NAME}_PREFIX}/lib${NAME})

ExternalProject_Add(
    ${LIB_NAME}
    SOURCE_DIR  ${${NAME}_SOURCE_DIR}
    PREFIX      ${${NAME}_PREFIX}
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${${NAME}_SOURCE_DIR}  ${${NAME}_BINARY_DIR}
    BINARY_DIR ${${NAME}_BINARY_DIR}
    BUILD_IN_SOURCE 0
    BUILD_COMMAND ${MAKE_COMMAND} -j CFLAGS=-fPIC libdeflate.a
    INSTALL_COMMAND ""
    LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
)

set(DEFLATELIB_BINARY_DIR ${${NAME}_BINARY_DIR})
set(${LIB_NAME}_LIBS ${${NAME}_BINARY_DIR}/libdeflate.a)

# build htslib
set(NAME hts)
set (LIB_NAME ${NAME}lib)
set(${NAME}_PREFIX ${PROJECT_BINARY_DIR}/${THIRD_PARTY_DIR_NAME}/${LIB_NAME})
set(${NAME}_SOURCE_DIR ${CMAKE_SOURCE_DIR}/${THIRD_PARTY_DIR_NAME}/${LIB_NAME})
set(${NAME}_BINARY_DIR ${${NAME}_PREFIX}/${LIB_NAME})
set(${NAME}_I_DIR ${${NAME}_BINARY_DIR})
set(${NAME}_L_DIR ${${NAME}_BINARY_DIR})

ExternalProject_Add(
    ${LIB_NAME}
    DEPENDS deflatelib
    PREFIX ${${NAME}_PREFIX}
    SOURCE_DIR  ${${NAME}_SOURCE_DIR}               
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E copy_directory ${${NAME}_SOURCE_DIR}  ${${NAME}_BINARY_DIR}
        COMMAND autoreconf
	COMMAND ${${NAME}_BINARY_DIR}/configure "CFLAGS=-I${DEFLATELIB_BINARY_DIR} -O2 -g" --disable-bz2 --disable-lzma --with-libdeflate --disable-gcs --disable-s3 --disable-libcurl --disable-plugins LDFLAGS=-L${DEFLATELIB_BINARY_DIR} --prefix=${${NAME}_BINARY_DIR}
    BINARY_DIR ${${NAME}_BINARY_DIR}
    BUILD_IN_SOURCE 0
    BUILD_COMMAND ${MAKE_COMMAND} -j
    LOG_DOWNLOAD 0
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
)

include_directories(BEFORE SYSTEM ${${NAME}_I_DIR})
#link_directories(${${NAME}_L_DIR})

set(htslib_LIBS
    "${${NAME}_BINARY_DIR}/libhts.a"
    ${deflatelib_LIBS}
    ${zlib_LIB}
    ${CMAKE_DL_LIBS}
    )
