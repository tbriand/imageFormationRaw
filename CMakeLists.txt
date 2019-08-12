cmake_minimum_required(VERSION 3.2)

# Enable C99
if(CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
else()
  set(CMAKE_C_STANDARD 99)
endif()

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

# include source code directory
set(SRC src)
include_directories(${SRC})

# include external code directory
set(EXTERNAL external)
include_directories(${EXTERNAL})

# include imscript directory
set(IMSCRIPT ${EXTERNAL}/imscript)
include_directories(${IMSCRIPT})

# include bspline directory
set(BSPLINE ${EXTERNAL}/bspline)
include_directories(${BSPLINE})

# include NFFT
add_subdirectory(${EXTERNAL}/nfft-3.5.0)
include_directories(${NFFT3_INCPATH})
link_directories(   ${NFFT3_LIBPATH})

set(PROGRAMS crop split_raw lowpass combi_ckr asymptotic_nc_filter concatenate_images ponomarenko_fit_raw_multiple channel_equalization affine0255 equalization create_burst add_noise rgb2raw compute difference_images compare_homography mean_and_std periodic_plus_smooth)

set(LIBS m jpeg png tiff fftw3_threads fftw3  ${GSL_LIBRARIES})

foreach(prg ${PROGRAMS})
    add_executable(${prg} src/main_${prg}.c ${IMSCRIPT}/iio.c ${BSPLINE}/splinter.c ${BSPLINE}/bspline.c src/bicubic.c src/combi_ckr_core.c src/fft_core.c src/filter_core.c src/homography_core.c src/interpolation_core.c src/periodic_plus_smooth.c src/ponomarenko_core.c src/tpi.c)
	add_dependencies(${prg} nfft-3.5.0)
    target_link_libraries(${prg} ${LIBS} nfft3_threads)
endforeach(prg ${PROGRAMS})
