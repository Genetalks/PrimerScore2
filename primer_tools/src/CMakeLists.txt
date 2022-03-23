set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(${TBB_INCLUDE_DIR})

# version related
EXECUTE_PROCESS(
    COMMAND git describe --tags master
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE PKG_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

if("${PKG_VERSION} " STREQUAL " ")
    set(PKG_VERSION "ptools-UNKNOWN")
endif()

set(VERSION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/version.h)
message(STATUS "version: " ${PKG_VERSION})

file(GENERATE OUTPUT ${VERSION_FILE}
    CONTENT "#define PTOOLS_VERSION \"${PKG_VERSION}\"\n"
    )

file(GLOB pt_src
    "*.cpp"
	"*.c"
)
set(PROG ptools)
add_executable(${PROG} ${pt_src})

add_dependencies(${PROG} CLI11lib deflatelib zlib htslib onetbblib)
set_property(TARGET ${PROG} PROPERTY CXX_STANDARD 11)
set_target_properties(${PROG} PROPERTIES INSTALL_RPATH "./")

target_link_libraries(${PROG} LINK_PUBLIC ${tbb_LIBS} ${htslib_LIBS} pthread )
