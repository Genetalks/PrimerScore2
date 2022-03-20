cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
SET(CMAKE_CXX_FLAGS_PROFILING "-DNDEBUG -O3 -g -fno-omit-frame-pointer" CACHE STRING "")
SET(CMAKE_C_FLAGS_PROFILING "-DNDEBUG -O3 -g -fno-omit-frame-pointer" CACHE STRING "")
MARK_AS_ADVANCED(CMAKE_CXX_FLAGS_PROFILING CMAKE_C_FLAGS_PROFILING)

project(primer_tools)

include(ExternalProject) 

add_compile_options("-Wall")
add_compile_options("-std=c++11")
add_compile_options("-m64")
add_compile_options("-fPIC")
#add_compile_options("-mavx2")
#add_compile_options("-mavx512f")
#add_compile_options("-mavx512vl")
# add_compile_options("-fprefetch-loop-arrays")
# add_compile_options("-mprefetchwt1")
# add_compile_options("-static-libgcc")
add_compile_options("-static-libstdc++")


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()


if (${USE_VMP_PROTECT})
    set (VMPROTECT_LIB ${PROJECT_SOURCE_DIR}/vendor/vmprotect/libVMProtectSDK64.so)
    message("--- build stitch vmprotect version ---")
    ADD_DEFINITIONS(-DUSE_VMP)
endif()

include("cmake/onetbb.cmake")
include("cmake/htslib.cmake")
include("cmake/CLI11.cmake")

include_directories(src)

add_subdirectory(src)

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release -DBUILD_DOCS=ON ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(profile
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Profiling ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Profile"
  )


