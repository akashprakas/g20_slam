#include "Eigen/Eigen"
#include <vector>

#ifndef ICP_H
#define ICP_H

#define N_pt 30          // # of points in the datasets
#define N_tests 100      // # of test iterations
#define noise_sigma 0.01 // standard deviation error to be added
#define translation 0.1  // max translation of the test set
#define rotation 0.1     // max rotation (radians) of the test set

// ICP output structure for both 2D and 3D
template <int Dim> struct ICP_OUT {
  Eigen::Matrix<double, Dim + 1, Dim + 1>
      trans; // Transformation matrix (3x3 for 2D, 4x4 for 3D)
  std::vector<float> distances; // Distances between corresponding points
  int iter;                     // Number of iterations taken
};

// Nearest neighbor output structure for both 2D and 3D
struct NEIGHBOR {
  std::vector<float> distances; // Nearest distances
  std::vector<int> indices;     // Indices of nearest points
};

// Function to calculate best fit transformation (2D and 3D)
template <int Dim>
Eigen::Matrix<double, Dim + 1, Dim + 1>
best_fit_transform(const Eigen::Matrix<double, Dim, Eigen::Dynamic> &A,
                   const Eigen::Matrix<double, Dim, Eigen::Dynamic> &B);

// ICP algorithm for both 2D and 3D
template <int Dim>
ICP_OUT<Dim> icp(const Eigen::Matrix<double, Dim, Eigen::Dynamic> &A,
                 const Eigen::Matrix<double, Dim, Eigen::Dynamic> &B,
                 int max_iterations = 20, double tolerance = 0.001);

// Nearest neighbor search function (works for both 2D and 3D)
template <int Dim>
NEIGHBOR
nearest_neighbor(const Eigen::Matrix<double, Dim, Eigen::Dynamic> &src,
                 const Eigen::Matrix<double, Dim, Eigen::Dynamic> &dst);

// Distance calculation function for 2D
float dist(const Eigen::Vector2d &pta, const Eigen::Vector2d &ptb);

// Distance calculation function for 3D
float dist(const Eigen::Vector3d &pta, const Eigen::Vector3d &ptb);

#endif
