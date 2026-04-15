cmake_minimum_required(VERSION 3.15...4.1)

project(
  radiative_transfer
  VERSION 1.0
  LANGUAGES CXX)

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)
set (CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3")


# Find packages go here.
find_package (Eigen3 REQUIRED NO_MODULE)
find_package(pybind11 REQUIRED)

pybind11_add_module(loki_python_binding
cpp_lib/constants.cpp
cpp_lib/CppInterfaceClass.cpp
cpp_lib/FixedHeightModelClass.cpp
cpp_lib/functions.cpp
cpp_lib/integrator.cpp
cpp_lib/interface_class_binding.cpp
cpp_lib/interpolator.cpp
cpp_lib/PSRClass.cpp
cpp_lib/read_write.cpp
cpp_lib/SolverClass.cpp)

target_link_libraries (loki_python_binding PRIVATE Eigen3::Eigen)