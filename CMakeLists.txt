cmake_minimum_required(VERSION 3.2)

set(GSL_FIND_QUIETLY TRUE)
find_package(GSL)
if(GSL_FOUND)
  include_directories("${GSL_INCLUDE_DIRS}")
  add_definitions(-DGSL_SUPPORT)
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Enable C99
if(CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else()
  set(CMAKE_C_STANDARD 99)
endif()

set(PROGRAMS crop split_raw lowpass combi_ckr asymptotic_nc_filter concatenate_images ponomarenko_fit_raw_multiple channel_equalization affine0255 equalization create_burst add_noise rgb2raw compute difference_images compare_homography mean_and_std)
set(LIBS m jpeg png tiff fftw3f fftw3 ${GSL_LIBRARIES})

set(IMSCRIPT src/imscript/)
include_directories(${IMSCRIPT})

set(EXTERNAL external/)
include_directories(${EXTERNAL})

foreach(prg ${PROGRAMS})
    add_executable(${prg} src/${prg}.c ${IMSCRIPT}/iio.c)
    target_link_libraries(${prg} ${LIBS})
endforeach(prg ${PROGRAMS})
