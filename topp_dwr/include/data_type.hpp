/************************************************************************
 * Software License Agreement (BSD License)
 *
 * @file data_type.hpp
 * @author liyong2018@zju.edu.cn
 * @brief basic data type for TOPP-DWR
 * @version  1.0.0
 * @date 2025-10.9
 * Copyright (c) 2025, CVTE.
 * All rights reserved.
 *     ______ ____  _____  _____        _______            ______
      |__   _/ __ \|  __ \|  __ \      |  __ \ \         /   __  \
        | | | |  | | |__) | |__) |_____| |  | \ \  /\  / /| |__) |
        | | | |  | |  ___/|  ___/______| |  | |\ \/  \/ / |  _  /
        | | | |__| | |    | |          | |__| | \  /\  /  | | \ \
        |_|  \____/|_|    |_|          |_____/   \/  \/   |_|  \_\

 ************************************************************************/
#pragma once
#include <cmath>
#include <vector>
namespace TOPP_DWR {
struct TrajectoryPoint {
  TrajectoryPoint() {}
  TrajectoryPoint(const TrajectoryPoint &obj) {
    this->x = obj.x;
    this->y = obj.y;
    this->yaw = obj.yaw;
    this->v = obj.v;
    this->w = obj.w;
    this->s = obj.s;
    this->time = obj.time;
    this->curve = obj.curve;
  }
  double x;
  double y;
  double yaw;
  double v;
  double w;
  double s;
  double time;
  double curve;
  double distanceTo(const TrajectoryPoint &other) const {
    double dx = other.x - x;
    double dy = other.y - y;
    return std::sqrt(dx * dx + dy * dy);
  }
  double angle_diff(const TrajectoryPoint &second) {
    double a = yaw;
    double b = second.yaw;
    double d1, d2;
    a = std::atan2(std::sin(a), std::cos(a));
    b = std::atan2(std::sin(b), std::cos(b));
    d1 = a - b;
    d2 = 2 * M_PI - std::fabs(d1);
    if (d1 > 0) {
      d2 *= -1.0;
    }
    if (fabs(d1) < std::fabs(d2)) {
      return (d1);
    } else {
      return (d2);
    }
  }
};
typedef std::vector<TrajectoryPoint> Trajectory;
struct TOPP_DWR_Params {
  double v_max;
  double a_max_l;
  double a_max_n;
  double vr_max;
  double vl_max;
  double d;
  double w_max;
  double start_vel;
  double end_vel;
  TOPP_DWR_Params()
      : v_max(0.6), a_max_l(1.0), a_max_n(0.6), vr_max(0.75), vl_max(0.75),
        d(0.35), w_max(2.0), start_vel(0.0), end_vel(2.0) {}
};

typedef struct {
  double x;
  double y;
  double z;
  double w;
} Quaternion;
static Quaternion yawToQuaternion(const double roll, const double pitch,
                                  const double yaw) {
  double halfYaw = yaw * 0.5;
  double halfPitch = pitch * 0.5;
  double halfRoll = roll * 0.5;
  double cosYaw = std::cos(halfYaw);
  double sinYaw = std::sin(halfYaw);
  double cosPitch = std::cos(halfPitch);
  double sinPitch = std::sin(halfPitch);
  double cosRoll = std::cos(halfRoll);
  double sinRoll = std::sin(halfRoll);
  Quaternion q;
  q.x = sinRoll * cosPitch * cosYaw - cosRoll * sinPitch * sinYaw;
  q.y = cosRoll * sinPitch * cosYaw + sinRoll * cosPitch * sinYaw;
  q.z = cosRoll * cosPitch * sinYaw - sinRoll * sinPitch * cosYaw;
  q.w = cosRoll * cosPitch * cosYaw + sinRoll * sinPitch * sinYaw;
  return q;
}
} // namespace TOPP_DWR