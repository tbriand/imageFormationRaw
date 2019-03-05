cmake_minimum_required(VERSION 3.0)

set (CMAKE_C_STANDARD 99)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(PROGRAMS crop split_raw lowpass combi_ckr asymptotic_nc_filter concatenate_images ponomarenko_fit_raw_multiple channel_equalization affine0255 equalization)
set(LIBS m jpeg png tiff fftw3f fftw3)

set(IMSCRIPT src/imscript/)
include_directories(${IMSCRIPT})

foreach(prg ${PROGRAMS})
    add_executable(${prg} src/${prg}.c ${IMSCRIPT}/iio.c)
    target_link_libraries(${prg} ${LIBS})
endforeach(prg ${PROGRAMS})
