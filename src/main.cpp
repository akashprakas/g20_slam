// #include "edge_se2.h"
// #include "se2.h"
// #include "vertex_se2.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cassert>

#include "g2o/core/base_binary_edge.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/hyper_graph_action.h"
#include "g2o/stuff/macros.h"
#include "g2o/stuff/misc.h"

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
};

class EdgeSE2 : public g2o::BaseBinaryEdge<3, SE2, VertexSE2, VertexSE2> {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  EdgeSE2();

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

  std::cout << "demo" << std::endl;
}