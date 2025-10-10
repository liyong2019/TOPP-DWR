/************************************************************************
 * Software License Agreement (BSD License)
 *
 * @file non_uniform_bspline.hpp
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
#include "data_type.hpp"
#include "eigen3/Eigen/Core"
namespace TOPP_DWR {
class NonUniformBspline {
private:
  Eigen::MatrixXd control_points_;
  int p_, n_, m_;     // p degree, n+1 control points, m = n+p+1
  Eigen::VectorXd u_; // knots vector
  double interval_;   // knot span \delta t
public:
  void getTimeSpan(double &um, double &um_p);
  Eigen::VectorXd getKnot();
  Eigen::VectorXd evaluateDeBoor(const double &u, int &k, int &k_upper);
  void setKnot(const Eigen::VectorXd &knot);
  NonUniformBspline(const Trajectory &path_data, const int &order,
                    const double &interval);
  NonUniformBspline(){};
};
} // namespace TOPP_DWR