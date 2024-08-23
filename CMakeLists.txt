cmake_minimum_required(VERSION 3.16)

project(POLLUX
  LANGUAGES Fortran
)

add_executable(pollux Polx_Xu.f npy.f90)
target_include_directories(pollux PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
