/************************************************************************
 * Software License Agreement (BSD License)
 *
 * @file non_uniform_bspline.cpp
 * @author liyong2018@zju.edu.cn
 * @brief bspline curve fitting for raw path
 * @version  1.0.0
 * @date 2025-10-9
 * Copyright (c) 2025, CVTE.
 * All rights reserved.
 *     ______ ____  _____  _____        _______            ______
      |__   _/ __ \|  __ \|  __ \      |  __ \ \         /   __  \
        | | | |  | | |__) | |__) |_____| |  | \ \  /\  / /| |__) |
        | | | |  | |  ___/|  ___/______| |  | |\ \/  \/ / |  _  /
        | | | |__| | |    | |          | |__| | \  /\  /  | | \ \
        |_|  \____/|_|    |_|          |_____/   \/  \/   |_|  \_\

 ************************************************************************/
#include "non_uniform_bspline.hpp"
namespace TOPP_DWR {
NonUniformBspline::NonUniformBspline(const Trajectory &path_data,
                                     const int &order, const double &interval) {
  Eigen::MatrixXd ctrl_pts(path_data.size() + 2, 3);
  {
    Eigen::Vector3d first_control_point =
        Eigen::Vector3d(2 * path_data[0].x - path_data[1].x,
                        2 * path_data[0].y - path_data[1].y, 0);
    ctrl_pts(0, 0) = first_control_point(0);
    ctrl_pts(0, 1) = first_control_point(1);
    ctrl_pts(0, 2) = 0;
    for (size_t i = 0; i < path_data.size(); ++i) {
      ctrl_pts(i + 1, 0) = path_data[i].x;
      ctrl_pts(i + 1, 1) = path_data[i].y;
      ctrl_pts(i + 1, 2) = 0;
    }
    Eigen::Vector3d last_control_point =
        Eigen::Vector3d(2 * path_data[path_data.size() - 1].x -
                            path_data[path_data.size() - 2].x,
                        2 * path_data[path_data.size() - 1].y -
                            path_data[path_data.size() - 2].y,
                        0);
    ctrl_pts(path_data.size() + 1, 0) = last_control_point(0);
    ctrl_pts(path_data.size() + 1, 1) = last_control_point(1);
    ctrl_pts(path_data.size() + 1, 2) = 0;
  }
  control_points_ = ctrl_pts;
  p_ = order;
  interval_ = interval;
  n_ = ctrl_pts.rows() - 1;
  m_ = n_ + p_ + 1;
  u_ = Eigen::VectorXd::Zero(m_ + 1);

  for (int i = 0; i <= u_.size() - 1; ++i) {
    if (i <= p_) {
      u_(i) = static_cast<double>(-p_ + i) * interval_;
    } else if ((i > p_) && (i <= u_.size() - 1 - p_)) {
      auto first = path_data.at(i - p_ - 1);
      auto second = path_data.at(i - p_);
      double dis = first.distanceTo(second);
      double interal = dis / 0.2 * interval_;
      u_(i) = u_(i - 1) + interal;
    } else {
      u_(i) = u_(i - 1) + interval_;
    }
  }
}
void NonUniformBspline::getTimeSpan(double &um, double &um_p) {
  um = u_(p_);
  um_p = u_(m_ - p_);
}
void NonUniformBspline::setKnot(const Eigen::VectorXd &knot) { u_ = knot; }
Eigen::VectorXd NonUniformBspline::getKnot() { return u_; }
Eigen::VectorXd NonUniformBspline::evaluateDeBoor(const double &u, int &k,
                                                  int &k_upper) {
  double ub = std::min(std::max(u_(p_), u), u_(m_ - p_));
  k_upper = std::min(k_upper, m_ - p_);

  // determine which [ui,ui+1] lay in
  int end = m_ - p_;
  end = k_upper;
  if (k < p_)
    k = p_;
  // Determine the interval using binary search
  while (u_(k + 1) < ub) {
    int mid = (k + end) / 2;
    if (u_(mid) > ub) {
      end = mid;
    } else {
      k = mid;
    }
  }
  k_upper = k + 1;

  /* deBoor's alg */
  std::vector<Eigen::VectorXd> d;
  for (int i = 0; i <= p_; ++i) {
    d.emplace_back(control_points_.row(k - p_ + i));
  }

  for (int r = 1; r <= p_; ++r) {
    for (int i = p_; i >= r; --i) {
      double alpha =
          (ub - u_[i + k - p_]) / (u_[i + 1 + k - r] - u_[i + k - p_]);
      d[i] = (1 - alpha) * d[i - 1] + alpha * d[i];
    }
  }
  return d[p_];
}
} // namespace TOPP_DWR