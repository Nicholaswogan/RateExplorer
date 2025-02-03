cmake_minimum_required(VERSION "3.14")
project(Clima LANGUAGES Fortran C VERSION "0.5.6")

# Download fractal sources
file(DOWNLOAD 
  "https://github.com/storyofthewolf/fractal_optics_coreshell/archive/276197ad01a5d29dcb8f1ddf0f6ed2b2cd6c1efd.tar.gz"
  "fractal_optics_coreshell.tar.gz"
)
file(ARCHIVE_EXTRACT 
  INPUT "fractal_optics_coreshell.tar.gz"
  DESTINATION ${CMAKE_BINARY_DIR}
)
set(FRACTAL_DIR "${CMAKE_BINARY_DIR}/fractal_optics_coreshell-276197ad01a5d29dcb8f1ddf0f6ed2b2cd6c1efd/src")

add_library(fractal SHARED
  "${FRACTAL_DIR}/adgaquad_mod.F90"
  "${FRACTAL_DIR}/adgaquad_types_mod.F90"
  "${FRACTAL_DIR}/fractal_meanfield_mod.F90"
  "${FRACTAL_DIR}/lusolvec_mod.F90"
  "${FRACTAL_DIR}/miess_mod.F90"
  "${FRACTAL_DIR}/system_mod.F90"
  src/fractal_meanfield_mod_c.F90
)

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(fractal PRIVATE -ffree-line-length-none -Wimplicit-interface -fimplicit-none)
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(fractal PRIVATE -fcheck=all,no-array-temps)
  endif()
endif()

install(TARGETS fractal DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/)