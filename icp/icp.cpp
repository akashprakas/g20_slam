#include "icp.h"
#include "Eigen/Eigen"
#include <iostream>
#include <numeric>

using namespace std;
using namespace Eigen;

// Best fit transformation for 2D points
Eigen::Matrix3d best_fit_transform(const Eigen::MatrixXd &A,
                                   const Eigen::MatrixXd &B) {
  /*
  Notice:
  1/ JacobiSVD return U,S,V, S as a vector, "use U*S*Vt" to get original Matrix;
  2/ matrix type 'MatrixXd' or 'MatrixXf' matters.
  */
  Eigen::Matrix3d T =
      Eigen::MatrixXd::Identity(3, 3); // Identity for 2D (3x3 matrix)
  Eigen::Vector2d centroid_A(0, 0);
  Eigen::Vector2d centroid_B(0, 0);
  Eigen::MatrixXd AA = A;
  Eigen::MatrixXd BB = B;
  int row = A.rows();

  // Calculate centroids
  for (int i = 0; i < row; i++) {
    centroid_A += A.block<1, 2>(i, 0).transpose();
    centroid_B += B.block<1, 2>(i, 0).transpose();
  }
  centroid_A /= row;
  centroid_B /= row;

  // Subtract centroids from point sets
  for (int i = 0; i < row; i++) {
    AA.block<1, 2>(i, 0) = A.block<1, 2>(i, 0) - centroid_A.transpose();
    BB.block<1, 2>(i, 0) = B.block<1, 2>(i, 0) - centroid_B.transpose();
  }

  // Compute the covariance matrix H
  Eigen::MatrixXd H = AA.transpose() * BB;
  Eigen::MatrixXd U, V, Vt;
  Eigen::VectorXd S;
  Eigen::Matrix2d R;
  Eigen::Vector2d t;

  // Perform Singular Value Decomposition (SVD)
  JacobiSVD<Eigen::MatrixXd> svd(H, ComputeFullU | ComputeFullV);
  U = svd.matrixU();
  S = svd.singularValues();
  V = svd.matrixV();
  Vt = V.transpose();

  // Compute rotation matrix R
  R = Vt.transpose() * U.transpose();

  // Ensure a proper rotation matrix (right-handed)
  if (R.determinant() < 0) {
    Vt.block<1, 2>(1, 0) *= -1;
    R = Vt.transpose() * U.transpose();
  }

  // Compute translation vector t
  t = centroid_B - R * centroid_A;

  // Combine into transformation matrix T
  T.block<2, 2>(0, 0) = R;
  T.block<2, 1>(0, 2) = t;
  return T;
}

// ICP for 2D points
ICP_OUT<2> icp(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B,
               int max_iterations, double tolerance) {
  int row = A.rows();
  Eigen::MatrixXd src =
      Eigen::MatrixXd::Ones(3, row); // 2D homogeneous coordinates (3xN)
  Eigen::MatrixXd src2d = Eigen::MatrixXd::Ones(2, row); // 2D points (2xN)
  Eigen::MatrixXd dst = Eigen::MatrixXd::Ones(3, row);
  NEIGHBOR neighbor;
  Eigen::Matrix3d T;
  Eigen::MatrixXd dst_chorder =
      Eigen::MatrixXd::Ones(2, row); // Matching points in B
  ICP_OUT<2> result;
  int iter = 0;

  // Initialize src and dst matrices with points
  for (int i = 0; i < row; i++) {
    src.block<2, 1>(0, i) = A.block<1, 2>(i, 0).transpose(); // 2D points
    src2d.block<2, 1>(0, i) = A.block<1, 2>(i, 0).transpose();
    dst.block<2, 1>(0, i) = B.block<1, 2>(i, 0).transpose();
  }

  double prev_error = 0;
  double mean_error = 0;
  for (int i = 0; i < max_iterations; i++) {
    // Find nearest neighbors between src and B
    neighbor = nearest_neighbor<2>(src2d.transpose(), B);

    // Match the nearest points from B to src
    for (int j = 0; j < row; j++) {
      dst_chorder.block<2, 1>(0, j) = dst.block<2, 1>(0, neighbor.indices[j]);
    }

    // Compute the transformation between src and the nearest points in B
    T = best_fit_transform(src2d.transpose(), dst_chorder.transpose());

    // Apply the transformation to src
    src = T * src;
    for (int j = 0; j < row; j++) {
      src2d.block<2, 1>(0, j) = src.block<2, 1>(0, j);
    }

    // Compute the mean error between the nearest points
    mean_error = std::accumulate(neighbor.distances.begin(),
                                 neighbor.distances.end(), 0.0) /
                 neighbor.distances.size();
    if (abs(prev_error - mean_error) < tolerance) {
      break;
    }
    prev_error = mean_error;
    iter = i + 1;
  }

  // Final transformation
  T = best_fit_transform(A, src2d.transpose());
  result.trans = T;
  result.distances = neighbor.distances;
  result.iter = iter;

  return result;
}

// Nearest neighbor search for 2D points
NEIGHBOR nearest_neighbor(const Eigen::MatrixXd &src,
                          const Eigen::MatrixXd &dst) {
  int row_src = src.rows();
  int row_dst = dst.rows();
  Eigen::Vector2d vec_src;
  Eigen::Vector2d vec_dst;
  NEIGHBOR neigh;
  float min = 100;
  int index = 0;
  float dist_temp = 0;

  // Find nearest neighbor for each point in src
  for (int ii = 0; ii < row_src; ii++) {
    vec_src = src.block<1, 2>(ii, 0).transpose();
    min = 100;
    index = 0;
    dist_temp = 0;
    for (int jj = 0; jj < row_dst; jj++) {
      vec_dst = dst.block<1, 2>(jj, 0).transpose();
      dist_temp = dist(vec_src, vec_dst);
      if (dist_temp < min) {
        min = dist_temp;
        index = jj;
      }
    }
    neigh.distances.push_back(min);
    neigh.indices.push_back(index);
  }

  return neigh;
}

// Distance function for 2D points
float dist(const Eigen::Vector2d &pta, const Eigen::Vector2d &ptb) {
  return sqrt((pta[0] - ptb[0]) * (pta[0] - ptb[0]) +
              (pta[1] - ptb[1]) * (pta[1] - ptb[1]));
}
