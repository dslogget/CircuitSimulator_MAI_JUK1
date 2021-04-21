cmake_minimum_required( VERSION 3.18 )

project( Matrix VERSION 0.1 )

configure_file( MatrixConfig.h.in MatrixConfig.h )

# file( GLOB SOURCES "*.cpp" )
add_executable( Matrix "test.cpp" )

target_include_directories( Matrix PUBLIC
                            "${PROJECT_SOURCE_DIR}/includes" )

add_executable( AD-Test "AutoDiffTest.cpp" )

target_include_directories( AD-Test PUBLIC
                            "${PROJECT_SOURCE_DIR}/includes" )

add_executable( FFT-Test "fft.cpp" )

target_include_directories( FFT-Test PUBLIC
                            "${PROJECT_SOURCE_DIR}/includes" )

file( GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/CircuitElements/*.cpp" )

add_executable( CircuitSimulator "CircuitSimulator.cpp" "${SOURCES}" )
#-----------------------------------------------------------------------------------------------
# Find Python
find_package(Python3 COMPONENTS Development NumPy )

if (Python3_FOUND)
   target_include_directories( CircuitSimulator PUBLIC "${Python3_INCLUDE_DIRS}"
                                                       "${Python3_NumPy_INCLUDE_DIRS}" )
   target_link_libraries( CircuitSimulator Python3::Python Python3::NumPy )
   target_compile_definitions( CircuitSimulator PUBLIC WITH_PYTHON )
   message(STATUS "Found Python3")
else()
   message(WARNING "Did not find Python3")
endif()

#-----------------------------------------------------------------------------------------------
# Find Matlab
find_package(Matlab)

if (Matlab_FOUND)
   target_include_directories( CircuitSimulator PUBLIC "${Matlab_ROOT_DIR}/extern/include" )
   target_compile_definitions( CircuitSimulator PUBLIC WITH_MATLAB )
   message(STATUS "Found Matlab")

   # Link C++ API
   if(WIN32)
      message(STATUS "Linking WINDOWS C++")
      target_link_libraries(CircuitSimulator ${Matlab_ROOT_DIR}/extern/lib/win64/microsoft/libMatlabEngine.lib)
      target_link_libraries(CircuitSimulator ${Matlab_ROOT_DIR}/extern/lib/win64/microsoft/libMatlabDataArray.lib)
   endif(WIN32)
   if(UNIX)
      message(STATUS "Linking UNIX C++")
      target_link_libraries(CircuitSimulator ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabEngine.so)
      target_link_libraries(CircuitSimulator ${Matlab_ROOT_DIR}/extern/bin/glnxa64/libMatlabDataArray.so)
      find_package(Threads)
      target_link_libraries(CircuitSimulator ${CMAKE_THREAD_LIBS_INIT})
   endif(UNIX)

   add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${PROJECT_SOURCE_DIR}/Matlab"
                "${PROJECT_BINARY_DIR}/Matlab")
   add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                "${PROJECT_SOURCE_DIR}/RunSimulatorMatlab.bat"
                "${PROJECT_BINARY_DIR}/RunSimulatorMatlab.bat")

else()
   message(WARNING "Did not find Matlab")
endif()


target_include_directories( CircuitSimulator PUBLIC
                            "${PROJECT_SOURCE_DIR}/includes" )


add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${PROJECT_SOURCE_DIR}/Netlists
                ${PROJECT_BINARY_DIR}/Netlists)
add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${PROJECT_SOURCE_DIR}/S-Parameters
                ${PROJECT_BINARY_DIR}/S-Parameters)
add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${PROJECT_SOURCE_DIR}/PRR
                ${PROJECT_BINARY_DIR}/PRR)
add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${PROJECT_SOURCE_DIR}/TimeSeries
                ${PROJECT_BINARY_DIR}/TimeSeries)
add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${PROJECT_SOURCE_DIR}/DataDumps
                ${PROJECT_BINARY_DIR}/DataDumps)
add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                "${PROJECT_SOURCE_DIR}/Circuit Simulator Run Transistors.bat"
                "${PROJECT_BINARY_DIR}/Circuit Simulator Run Transistors.bat")
add_custom_command(
        TARGET CircuitSimulator PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                "${PROJECT_SOURCE_DIR}/Run Benchmarks.bat"
                "${PROJECT_BINARY_DIR}/Run Benchmarks.bat")


# specify the C++ standard
set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_STANDARD_REQUIRED True )

set( CMAKE_CXX_FLAGS  "-g --std=c++20 -Wall -Wpedantic -Werror -march=native -O3" )
set( CMAKE_EXPORT_COMPILE_COMMANDS ON )

