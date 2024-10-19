#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cassert>

#include <g2o/config.h>
#include <g2o/core/base_binary_edge.h>
#include <g2o/core/base_vertex.h>
#include <g2o/core/hyper_graph_action.h>
#include <g2o/stuff/macros.h>
#include <g2o/stuff/misc.h>

#include <g2o/core/block_solver.h>
#include <g2o/core/factory.h>
#include <g2o/core/optimization_algorithm_factory.h>
#include <g2o/core/optimization_algorithm_gauss_newton.h>
#include <g2o/core/sparse_optimizer.h>
#include <g2o/solvers/eigen/linear_solver_eigen.h>

#include "edge_se2.h"
#include "se2.h"
#include "vertex_se2.h"

#include <fstream>

void saveGraph(const g2o::SparseOptimizer &optimizer,
               const std::string &filename) {
  std::ofstream outFile(filename);
  if (!outFile.is_open()) {
    std::cerr << "Failed to open file for writing: " << filename << std::endl;
    return;
  }

  // Save vertices
  for (const auto &vertexPair : optimizer.vertices()) {
    const auto *vertex = dynamic_cast<const VertexSE2 *>(vertexPair.second);
    if (vertex) {
      outFile << "VERTEX_SE2 " << vertex->id() << " "
              << vertex->estimate().translation().x() << " "
              << vertex->estimate().translation().y() << " "
              << vertex->estimate().rotation().angle() << std::endl;
    }
  }

  // Save edges
  for (const auto &edge : optimizer.edges()) {
    const auto *edgeSE2 = dynamic_cast<const EdgeSE2 *>(edge);
    if (edgeSE2) {
      outFile << "EDGE_SE2 " << edgeSE2->vertices()[0]->id() << " "
              << edgeSE2->vertices()[1]->id() << " "
              << edgeSE2->measurement().translation().x() << " "
              << edgeSE2->measurement().translation().y() << " "
              << edgeSE2->measurement().rotation().angle() << " ";
      for (int i = 0; i < 3; ++i) {
        for (int j = i; j < 3; ++j) {
          outFile << edgeSE2->information()(i, j) << " ";
        }
      }
      outFile << std::endl;
    }
  }

  outFile.close();
  std::cout << "Graph saved to " << filename << std::endl;
}

struct VertexHolder {
  int id;
  double x;
  double y;
  double theta;
};

struct EdgeHolder {
  int from;
  int to;
  float dx;
  float dy;
  float dtheta;
  float i11, i12, i13, i22, i23, i33;
};

int main() {
  std::vector<VertexHolder> vertices;
  std::vector<EdgeHolder> edges;
  std::ifstream file("/home/akash/Documents/code/g2o_slam/input_INTEL_g2o.g2o");
  if (!file.is_open()) {
    std::cerr << "Failed to open the file." << std::endl;
  }

  std::string line;
  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string keyword;
    ss >> keyword;
    if (keyword == "VERTEX_SE2") {
      int id;
      float x, y, theta;
      ss >> id >> x >> y >> theta;
      VertexHolder temp{id, x, y, theta};
      vertices.push_back(temp);
    }
  }
  std::ifstream file2(
      "/home/akash/Documents/code/g2o_slam/input_INTEL_g2o.g2o");
  if (!file2.is_open()) {
    std::cerr << "Failed to open the file." << std::endl;
  }
  while (std::getline(file2, line)) {
    std::stringstream ss(line);
    std::string keyword;
    ss >> keyword;
    if (keyword == "EDGE_SE2") {
      int from, to;
      float dx, dy, dtheta;
      float i11, i12, i13, i22, i23, i33;
      ss >> from >> to >> dx >> dy >> dtheta >> i11 >> i12 >> i13 >> i22 >>
          i23 >> i33;
      EdgeHolder temp{from, to, dx, dy, dtheta, i11, i12, i13, i22, i23, i33};
      edges.push_back(temp);
    }
  }

  typedef g2o::BlockSolver<g2o::BlockSolverTraits<3, 3>> SlamBlockSolver;
  typedef g2o::LinearSolverEigen<SlamBlockSolver::PoseMatrixType>
      SlamLinearSolver;

  // // allocating the optimizer

  g2o::SparseOptimizer optimizer;
  auto linearSolver = std::make_unique<SlamLinearSolver>();
  linearSolver->setBlockOrdering(false);
  g2o::OptimizationAlgorithmGaussNewton *solver =
      new g2o::OptimizationAlgorithmGaussNewton(
          std::make_unique<SlamBlockSolver>(std::move(linearSolver)));

  optimizer.setAlgorithm(solver);

  std::cout << "Adding vertex" << std::endl;

  for (size_t i = 0; i < vertices.size(); ++i) {
    const VertexHolder &vertex = vertices[i];
    SE2 temp{vertex.x, vertex.y, vertex.theta};
    VertexSE2 *robot = new VertexSE2;
    robot->setId(vertex.id);
    robot->setEstimate(temp);
    optimizer.addVertex(robot);
  }
  std::cerr << "Optimization: Adding odometry measurements ... ";

  for (size_t i = 0; i < edges.size(); ++i) {
    const auto &edge = edges[i];
    SE2 tranformation{edge.dx, edge.dy, edge.dtheta};
    Eigen::Matrix3d information;

    information(0, 0) = edge.i11;
    information(0, 1) = 0;
    information(0, 2) = 0;
    information(1, 0) = 0;
    information(1, 1) = edge.i22;
    information(1, 2) = 0;
    information(2, 0) = 0;
    information(2, 1) = 0;
    information(2, 2) = edge.i33;

    EdgeSE2 *odometry = new EdgeSE2;
    odometry->vertices()[0] = optimizer.vertex(edge.from);
    odometry->vertices()[1] = optimizer.vertex(edge.to);
    odometry->setMeasurement(tranformation);
    odometry->setInformation(information);
    optimizer.addEdge(odometry);
  }
  // dump initial state to the disk
  // optimizer.save("tutorial_before.g2o");
  saveGraph(optimizer, "before_optimization.g2o");

  VertexSE2 *firstRobotPose = dynamic_cast<VertexSE2 *>(optimizer.vertex(0));
  firstRobotPose->setFixed(true);
  optimizer.setVerbose(true);

  std::cerr << "Optimizing" << std::endl;
  std::cout << "Initial chi2: " << optimizer.chi2() << std::endl;
  optimizer.initializeOptimization();
  optimizer.optimize(10);
  std::cout << "Final chi2: " << optimizer.chi2() << std::endl;
  std::cerr << "done." << std::endl;

  saveGraph(optimizer, "after_optimization.g2o");

  // freeing the graph memory
  optimizer.clear();

  return 0;
}
