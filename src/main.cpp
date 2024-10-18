// #include "edge_se2.h"
// #include "se2.h"
// #include "vertex_se2.h"
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

class SE2 {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  SE2() : _R(0), _t(0, 0) {}

  SE2(double x, double y, double theta) : _R(theta), _t(x, y) {}

  const Eigen::Vector2d &translation() const { return _t; }

  Eigen::Vector2d &translation() { return _t; }

  const Eigen::Rotation2Dd &rotation() const { return _R; }

  Eigen::Rotation2Dd &rotation() { return _R; }

  SE2 operator*(const SE2 &tr2) const {
    SE2 result(*this);
    result._t += _R * tr2._t;
    result._R.angle() += tr2._R.angle();
    result._R.angle() = g2o::normalize_theta(result._R.angle());
    return result;
  }

  SE2 &operator*=(const SE2 &tr2) {
    _t += _R * tr2._t;
    _R.angle() += tr2._R.angle();
    _R.angle() = g2o::normalize_theta(_R.angle());
    return *this;
  }

  Eigen::Vector2d operator*(const Eigen::Vector2d &v) const {
    return _t + _R * v;
  }

  SE2 inverse() const {
    SE2 ret;
    ret._R = _R.inverse();
    ret._R.angle() = g2o::normalize_theta(ret._R.angle());
    ret._t = ret._R * (Eigen::Vector2d(-1 * _t));
    return ret;
  }

  double operator[](int i) const {
    assert(i >= 0 && i < 3);
    if (i < 2)
      return _t(i);
    return _R.angle();
  }

  double &operator[](int i) {
    assert(i >= 0 && i < 3);
    if (i < 2)
      return _t(i);
    return _R.angle();
  }

  void fromVector(const Eigen::Vector3d &v) { *this = SE2(v[0], v[1], v[2]); }

  Eigen::Vector3d toVector() const {
    Eigen::Vector3d ret;
    for (int i = 0; i < 3; i++) {
      ret(i) = (*this)[i];
    }
    return ret;
  }

protected:
  Eigen::Rotation2Dd _R;
  Eigen::Vector2d _t;
};

class VertexSE2 : public g2o::BaseVertex<3, SE2> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  VertexSE2() : g2o::BaseVertex<3, SE2>(){};

  virtual void setToOriginImpl() { _estimate = SE2(); }

  virtual void oplusImpl(const double *update) {
    SE2 up(update[0], update[1], update[2]);
    _estimate *= up;
  }
  virtual bool read(std::istream &is) {
    Eigen::Vector3d p;
    is >> p[0] >> p[1] >> p[2];
    _estimate.fromVector(p);
    return true;
  };
  virtual bool write(std::ostream &os) const {
    Eigen::Vector3d p = estimate().toVector();
    os << p[0] << " " << p[1] << " " << p[2];
    return os.good();
  };
};

class EdgeSE2 : public g2o::BaseBinaryEdge<3, SE2, VertexSE2, VertexSE2> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  EdgeSE2() : BaseBinaryEdge<3, SE2, VertexSE2, VertexSE2>(){};

  void computeError() {
    const VertexSE2 *v1 = static_cast<const VertexSE2 *>(_vertices[0]);
    const VertexSE2 *v2 = static_cast<const VertexSE2 *>(_vertices[1]);
    SE2 delta =
        _inverseMeasurement * (v1->estimate().inverse() * v2->estimate());
    _error = delta.toVector();
  }

  void setMeasurement(const SE2 &m) {
    _measurement = m;
    _inverseMeasurement = m.inverse();
  }

  virtual bool read(std::istream &is) {
    Eigen::Vector3d p;
    is >> p[0] >> p[1] >> p[2];
    _measurement.fromVector(p);
    _inverseMeasurement = measurement().inverse();
    for (int i = 0; i < 3; ++i)
      for (int j = i; j < 3; ++j) {
        is >> information()(i, j);
        if (i != j)
          information()(j, i) = information()(i, j);
      }
    return true;
  };

  virtual bool write(std::ostream &os) const {
    Eigen::Vector3d p = measurement().toVector();
    os << p.x() << " " << p.y() << " " << p.z();
    for (int i = 0; i < 3; ++i)
      for (int j = i; j < 3; ++j)
        os << " " << information()(i, j);
    return os.good();
  };

protected:
  SE2 _inverseMeasurement;
};

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
  optimizer.save("tutorial_before.g2o");

  VertexSE2 *firstRobotPose = dynamic_cast<VertexSE2 *>(optimizer.vertex(0));
  firstRobotPose->setFixed(true);
  optimizer.setVerbose(true);

  std::cerr << "Optimizing" << std::endl;
  std::cout << "Initial chi2: " << optimizer.chi2() << std::endl;
  optimizer.initializeOptimization();
  optimizer.optimize(10);
  std::cout << "Final chi2: " << optimizer.chi2() << std::endl;
  std::cerr << "done." << std::endl;

  optimizer.save("tutorial_after.g2o");

  // freeing the graph memory
  optimizer.clear();

  return 0;
}

// struct EdgeHolder {
//   int from;
//   int to;
//   float dx;
//   float dy;
//   float dtheta;
//   float i11, i12, i13, i22, i23, i33;
// };