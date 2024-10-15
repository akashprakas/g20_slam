cmake_minimum_required(VERSION 3.22)
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

# Add source files (your C++ files)
set(SOURCES
    src/main.cpp

)

# Create the executable target
add_executable(${PROJECT_NAME} ${SOURCES})

# Link G2O and Eigen libraries to your project
# target_link_libraries(${PROJECT_NAME} PUBLIC  g2o_core g2o_stuff g2o_solver_eigen)

target_link_libraries(${PROJECT_NAME} 
    ${g2o_CORE_LIBRARY} 
    ${g2o_STUFF_LIBRARY}
    ${g2o_SOLVER_CSPARSE}  # Add any additional solvers you need here
    ${g2o_SOLVER_DENSE}    # For example, CSPARSE and DENSE
    ${g2o_TYPES_SLAM2D}    # Add types as per your application
    ${g2o_TYPES_SLAM3D}
    ${EIGEN3_LIBRARIES}
)
