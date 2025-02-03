cmake_minimum_required(VERSION "3.14")
project(Clima LANGUAGES Fortran VERSION "0.5.6")

add_library(fractal SHARED
  src/adgaquad_mod.F90
  src/adgaquad_types_mod.F90
  src/fractal_meanfield_mod.F90
  src/lusolvec_mod.F90
  src/miess_mod.F90
  src/system_mod.F90
  src/fractal_meanfield_mod_c.F90
)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(fractal PRIVATE -ffree-line-length-none -Wimplicit-interface -fimplicit-none)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(fractal PRIVATE -fcheck=all,no-array-temps)
  endif()
endif()

install(TARGETS fractal DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/)