cmake_minimum_required(VERSION 3.20)
project(g2oSlamExample)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Set the C++ standard to C++11 or higher
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(G2O_DIR "/usr/local/lib/cmake/g2o")

# Find G2O package
find_package(g2o REQUIRED)
include_directories(${g2o_INCLUDE_DIRS})

# Find other dependencies like Eigen (G2O requires Eigen)
find_package(Eigen3 REQUIRED)

# Include directories for G2O and Eigen
include_directories(${EIGEN3_INCLUDE_DIR})

# Create the executable target
add_executable(${PROJECT_NAME} main.cpp se2.h vertex_se2.h vertex_se2.cpp edge_se2.h edge_se2.cpp)

# Link G2O and Eigen libraries to your project
target_link_libraries(${PROJECT_NAME} PUBLIC  g2o::core  g2o::solver_eigen g2o::stuff  ${EIGEN3_LIBRARIES}) 

