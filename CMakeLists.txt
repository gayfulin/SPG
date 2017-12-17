cmake_minimum_required(VERSION 2.8.5)
project(SPG Fortran)
set(CMAKE_Fortran_FLAGS "-fdefault-real-8")
add_library(rng STATIC
            Src/rng/ggnml.f
            Src/rng/ggubs.f
            Src/rng/mdnris.f
            Src/rng/merfi.f
            Src/rng/uertst.f
            Src/rng/ugetio.f)
add_library(generator STATIC
            Src/fft99f.f
            Src/generator.f90
            Src/set_gen_params.f90)
target_link_libraries(generator rng)

add_executable(SPG.exe Src/main.f90)
target_link_libraries(SPG.exe generator)

install(TARGETS SPG.exe
        RUNTIME DESTINATION bin)
include (CTest)
add_test(NAME SPG-testrun
  COMMAND $<TARGET_FILE:SPG.exe> ${CMAKE_SOURCE_DIR}/Config/testrun.cfg)

# A verbose test command with valgrind checker
add_custom_target(check ${CMAKE_COMMAND} -E env CTEST_OUTPUT_ON_FAILURE=1
                  ${CMAKE_CTEST_COMMAND} -T memcheck -C $<CONFIG> --verbose
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

