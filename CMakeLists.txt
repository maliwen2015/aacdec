cmake_minimum_required(VERSION 3.10)
project(aacdec VERSION 1.7.0 LANGUAGES C)

option(AACDEC_BUILD_EXAMPLE "Build the command line decoder" ON)
option(AACDEC_BUILD_TESTS "Build the test programs" ON)
option(AACDEC_ENABLE_LTO "Enable link time optimisation" OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(AACDEC_SOURCES
    src/aacdec.c
    src/aac_core.c
    src/aac_usac.c
    src/aac_usac_ac.c
    src/aac_tables.c
    src/aac_huffman.c
    src/aac_ps.c
    src/aac_sbr.c)

add_library(aacdec ${AACDEC_SOURCES})
add_library(AACDEC::aacdec ALIAS aacdec)

target_include_directories(aacdec
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

if(BUILD_SHARED_LIBS)
    target_compile_definitions(aacdec PUBLIC AACDECLIB_SHARED)
    target_compile_definitions(aacdec PRIVATE AACDECLIB_BUILD)
    set_target_properties(aacdec PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION 1)
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(aacdec PRIVATE -Wall -Wextra)
endif()

if(UNIX)
    target_link_libraries(aacdec PRIVATE m)
endif()

if(AACDEC_ENABLE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT AACDEC_IPO_OK OUTPUT AACDEC_IPO_MSG)
    if(AACDEC_IPO_OK)
        set_property(TARGET aacdec PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(WARNING "LTO is not supported: ${AACDEC_IPO_MSG}")
    endif()
endif()

if(AACDEC_BUILD_EXAMPLE)
    add_executable(aacdec_example example_decoder.c)
    target_link_libraries(aacdec_example PRIVATE aacdec)
    if(UNIX)
        target_link_libraries(aacdec_example PRIVATE m)
    endif()
    set_target_properties(aacdec_example PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

if(AACDEC_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

install(TARGETS aacdec
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
install(FILES aacdec.h DESTINATION include)
