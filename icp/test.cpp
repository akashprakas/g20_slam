#include "icp2d.h"
#include <Eigen/Dense>
#include <cmath>
#include <fstream> // Include for file handling
#include <iostream>
// Function to convert degrees to radians
double degreesToRadians(double degrees) { return degrees * (M_PI / 180.0); }

int main() {
  double angle = 45; // this is in degrees.

  // Convert degrees to radians
  double rad = degreesToRadians(static_cast<double>(angle));
  Eigen::Rotation2D<double> R_true(rad);
  Eigen::Vector2d t_true;
  t_true << -2.0, 5.0;

  const int num_points = 30;
  Eigen::Matrix<double, 2, num_points> true_data;
  Eigen::Matrix<double, 2, num_points> transformed_data;

  // Generate data
  for (int i = 0; i < num_points; ++i) {
    true_data(0, i) = static_cast<double>(i); // x-coordinates
    true_data(1, i) =
        0.2f * true_data(0, i) * sin(0.5f * true_data(0, i)); // y-coordinates
  }

  // Apply rotation and translation
  for (int i = 0; i < num_points; ++i) {
    Eigen::Vector2d point = R_true * true_data.col(i) + t_true;
    transformed_data.col(i) = point;
  }

  std::ofstream outfile("matrix_log.txt");
  if (outfile.is_open()) {
    // Log true_data
    outfile << "True Data:\n" << true_data << "\n\n";

    // Log transformed_data
    outfile << "Transformed Data:\n" << transformed_data << "\n";

    // Close the file
    outfile.close();
    std::cout << "Data successfully logged to matrix_log.txt" << std::endl;
  } else {
    std::cerr << "Error opening file for logging!" << std::endl;
  }

  // Eigen::MatrixXd A(2, N_pt); // 2D points (2xN matrix)
  // Eigen::MatrixXd B(2, N_pt);
  auto result = icp(true_data.transpose(), transformed_data.transpose());

  std::cout << result.trans << std::endl;

  // Eigen::MatrixXd A(3, N_pt); // 3D points (3xN matrix)
  // Eigen::MatrixXd B(3, N_pt);
  // auto result = icp<3>(A, B);

  return 0;
}
