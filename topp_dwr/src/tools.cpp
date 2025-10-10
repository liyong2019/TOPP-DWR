/************************************************************************
 * Software License Agreement (BSD License)
 *
 * @file topps.cpp
 * @author liyong2018@zju.edu.cn
 * @brief tools for topp_dwr
 * @version  1.0
 * @date 2025-09-30
 * Copyright (c) 2025, CVTE.
 * All rights reserved.
 *     ______ ____  _____  _____        _______            ______
      |__   _/ __ \|  __ \|  __ \      |  __ \ \         /   __  \
        | | | |  | | |__) | |__) |_____| |  | \ \  /\  / /| |__) |
        | | | |  | |  ___/|  ___/______| |  | |\ \/  \/ / |  _  /
        | | | |__| | |    | |          | |__| | \  /\  /  | | \ \
        |_|  \____/|_|    |_|          |_____/   \/  \/   |_|  \_\

 ************************************************************************/
#include "tools.hpp"
#include "ament_index_cpp/get_package_prefix.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
namespace TOPP_DWR {
static bool isCommonTrajType(const std::string &traj_type) {
  return (traj_type == "c1" || traj_type == "c2" || traj_type == "c3");
}

// 2. Implement basic scenario-based parameter reading function
bool Tools::readTOPP_DWR_Params(const std::string &param_file_path,
                                const std::string &traj_type,
                                TOPP_DWR_Params &out_params) {
  // Step 1: Check if the parameter file exists
  if (!std::filesystem::exists(param_file_path)) {
    std::cerr << "[Error] Param file not found: " << param_file_path
              << std::endl;
    return false;
  }
  try {
    // Step 2: Load the YAML file
    YAML::Node config = YAML::LoadFile(param_file_path);
    // Step 3: Select parameter group based on trajectory type
    YAML::Node target_param_group;
    if (isCommonTrajType(traj_type)) {
      // For c1/c2/c3: Read common parameter group
      if (!config["common_params"].IsDefined()) {
        std::cerr << "[Error] Common param group 'common_params' not found in "
                  << param_file_path << std::endl;
        return false;
      }
      target_param_group = config["common_params"];
      std::cout << "[Info] Reading common params (for " << traj_type << ")"
                << std::endl;
    } else if (traj_type == "c4") {
      // For c4: Read RViz exclusive parameter group
      if (!config["rviz_click_params"].IsDefined()) {
        std::cerr
            << "[Error] RViz param group 'rviz_click_params' not found in "
            << param_file_path << std::endl;
        return false;
      }
      target_param_group = config["rviz_click_params"];
      std::cout << "[Info] Reading RViz exclusive params (for c4)" << std::endl;
    } else {
      std::cerr << "[Error] Invalid trajectory type: " << traj_type
                << " (only support c1/c2/c3/c4)" << std::endl;
      return false;
    }

    // Step 4: Read specific fields in the parameter group (use default values
    // if fields are missing)
    if (target_param_group["v_max"].IsDefined())
      out_params.v_max = target_param_group["v_max"].as<double>();
    if (target_param_group["a_max_l"].IsDefined())
      out_params.a_max_l = target_param_group["a_max_l"].as<double>();
    if (target_param_group["a_max_n"].IsDefined())
      out_params.a_max_n = target_param_group["a_max_n"].as<double>();
    if (target_param_group["vr_max"].IsDefined())
      out_params.vr_max = target_param_group["vr_max"].as<double>();
    if (target_param_group["vl_max"].IsDefined())
      out_params.vl_max = target_param_group["vl_max"].as<double>();
    if (target_param_group["d"].IsDefined())
      out_params.d = target_param_group["d"].as<double>();
    if (target_param_group["w_max"].IsDefined())
      out_params.w_max = target_param_group["w_max"].as<double>();
    if (target_param_group["start_vel"].IsDefined())
      out_params.start_vel = target_param_group["start_vel"].as<double>();
    if (target_param_group["end_vel"].IsDefined())
      out_params.end_vel = target_param_group["end_vel"].as<double>();

    // Step 5: Print the read parameters (for debugging)
    std::cout << "[Info] Successfully read params:" << std::endl;
    std::cout << "v_max: " << out_params.v_max << std::endl;
    std::cout << "a_max_l: " << out_params.a_max_l << std::endl;
    std::cout << "a_max_n: " << out_params.a_max_n << std::endl;
    std::cout << "vr_max: " << out_params.vr_max << std::endl;
    std::cout << "vl_max: " << out_params.vl_max << std::endl;
    std::cout << "d: " << out_params.d << std::endl;
    std::cout << "w_max: " << out_params.w_max << std::endl;
    std::cout << "start_vel: " << out_params.start_vel << std::endl;
    std::cout << "end_vel: " << out_params.end_vel << std::endl;
    return true;
  } catch (const YAML::Exception &e) {
    std::cerr << "[Error] YAML parse failed: " << e.what() << std::endl;
    return false;
  } catch (const std::exception &e) {
    std::cerr << "[Error] Param read failed: " << e.what() << std::endl;
    return false;
  }
}
// Implement scenario-based function that reads parameters from package
// directory
bool Tools::readTOPP_DWR_ParamsFromPackage(const std::string &traj_type,
                                           TOPP_DWR_Params &out_params) {
  try {
    std::string pkg_prefix = ament_index_cpp::get_package_prefix("topp_dwr");
    std::string param_file_path = pkg_prefix + "/params/params.yaml";
    return readTOPP_DWR_Params(param_file_path, traj_type, out_params);
  } catch (const std::exception &e) {
    std::cerr << "[Error] Get package prefix failed: " << e.what() << std::endl;
    return false;
  }
}
bool Tools::parsePoint(const YAML::Node &point, Trajectory &reference_path) {
  TrajectoryPoint trajectory_point;
  // Parse x and y coordinates
  if (!point["x"] || !point["y"]) {
    std::cout << "Point missing x or y coordinate" << std::endl;
    return false;
  }
  trajectory_point.x = point["x"].as<double>();
  trajectory_point.y = point["y"].as<double>();
  trajectory_point.yaw = point["yaw"].as<double>(0.0);

  reference_path.emplace_back(trajectory_point);
  return true;
}
bool Tools::getPathFromFile(const std::string &filepath,
                            Trajectory &reference_path) {
  try {
    YAML::Node doc = YAML::LoadFile(filepath);
    const YAML::Node &point_node = doc["point"];

    if (!point_node || point_node.size() == 0) {
      std::cout << "Path file error: empty or missing 'point' node"
                << std::endl;
      return false;
    }
    size_t path_size = point_node.size();
    reference_path.reserve(path_size);

    for (const auto &point : point_node) {
      if (!parsePoint(point, reference_path)) {
        return false;
      }
    }
    std::cout << " reference_path size: " << reference_path.size() << std::endl;
    {
      auto tmp_path = reference_path;
      tmp_path.clear();
      tmp_path.emplace_back(reference_path.front());
      for (size_t i = 1; i < reference_path.size(); ++i) {
        double dis = reference_path.at(i).distanceTo(tmp_path.back());
        if (dis > 0.05) {
          tmp_path.emplace_back(reference_path.at(i));
        } else {
          std::cout << " reducted path point: " << i << std::endl;
        }
      }
      std::cout << " tmp_path size: " << tmp_path.size() << std::endl;
      reference_path = tmp_path;
    }

    return true;
  } catch (const YAML::Exception &e) {
    std::cout << "YAML parsing error: " << e.what() << " " << filepath
              << std::endl;
  } catch (const std::exception &e) {
    std::cout << "Error reading path file: " << std::endl;
  }
  return false;
}
bool Tools::loadCurve(const std::string file, Trajectory &trajectory) {
  std::string yaml_path;
  try {
    std::string pkg_path = ament_index_cpp::get_package_prefix("topp_dwr");
    yaml_path = pkg_path + "/files/" + file + ".yaml";
  } catch (const std::exception &e) {
    std::cerr << "[error] get_package_prefix failed" << e.what() << std::endl;
    return false;
  }
  bool succ = getPathFromFile(yaml_path, trajectory);
  if (!succ)
    return false;

  NonUniformBspline bspline;
  double interval = 1.0 / 3.0;
  bspline = NonUniformBspline(trajectory, 3, interval);
  trajectory = getRawPathWithBspline(bspline, 0.05);
  return succ;
}
bool Tools::generateRVIZPath(Trajectory &trajectory) {
  Trajectory interpolated_traj; // Store interpolated trajectory points
  const double STEP =
      1.0; // Interpolation resolution: 0.2m (core modification point)

  // Iterate through each pair of adjacent original clicked points for linear
  // interpolation
  for (size_t i = 0; i < trajectory.size() - 1; ++i) {
    const auto &p_prev = trajectory[i];     // Previous original point
    const auto &p_curr = trajectory[i + 1]; // Next original point

    // Calculate distance and increment between the two points
    double dx = p_curr.x - p_prev.x; // X-direction increment
    double dy = p_curr.y - p_prev.y; // Y-direction increment
    double dist_between_points =
        std::sqrt(dx * dx + dy * dy); // Distance between the two points

    // Calculate the number of interpolation steps needed for the current
    // segment (avoid division by 0, minimum 1 step)
    int num_steps = static_cast<int>(dist_between_points / STEP);
    num_steps = (num_steps < 1)
                    ? 1
                    : num_steps; // Ensure at least 1 step (prevent no
                                 // interpolation when points are too close)

    // Generate interpolated points step by step
    for (int j = 0; j < num_steps; ++j) {
      double ratio =
          static_cast<double>(j) / num_steps; // Interpolation ratio (0~1)
      TrajectoryPoint new_point;

      // 1. Position interpolation (linear interpolation for x/y)
      new_point.x = p_prev.x + ratio * dx;
      new_point.y = p_prev.y + ratio * dy;

      // Add to interpolated trajectory
      interpolated_traj.push_back(new_point);
    }
  }
  // Replace original trajectory with interpolated trajectory (free memory of
  // original points)
  trajectory.swap(interpolated_traj);
  for (size_t i = 0; i < trajectory.size(); ++i) {
    std::cout << " interpolated_traj: " << i << " x " << trajectory.at(i).x
              << " y " << trajectory.at(i).y << std::endl;
  }
  NonUniformBspline bspline;
  double interval = 1.0 / 3.0;
  bspline = NonUniformBspline(trajectory, 3, interval);
  trajectory = getRawPathWithBspline(bspline, 0.05);
  return true;
}

Trajectory Tools::getRawPathWithBspline(NonUniformBspline &bspline,
                                        double delta_t) {
  Trajectory trajectory;
  double tm = 0.0, tmp = 0.0;
  bspline.getTimeSpan(tm, tmp);
  double path_time = tmp - tm;
  int k_pos = 3; // 3rd-order B-spline
  int k_pos_upper = k_pos + 1;
  Eigen::Vector3d last_pt = bspline.evaluateDeBoor(0, k_pos, k_pos_upper);
  double last_length = 0;
  trajectory.reserve(static_cast<int>(path_time / delta_t));
  std::cout << " bsplineee start_u: " << tm << " end_u: " << tmp;
  bool end = false;
  Eigen::VectorXd pos_knot = bspline.getKnot();
  Eigen::Vector3d vel;
  Eigen::Vector3d pt;
  for (double t = tm; t <= tmp + delta_t; t += delta_t) {
    if (end)
      break;
    if (t > tmp) {
      t = tmp;
      end = true;
    }
    for (int i = k_pos; i < pos_knot.rows(); ++i) {
      if (pos_knot(i, 0) >= t) {
        k_pos_upper = i;
        break;
      }
    }
    pt = bspline.evaluateDeBoor(t, k_pos, k_pos_upper);
    double dis = (pt.block<2, 1>(0, 0) - last_pt.block<2, 1>(0, 0)).norm();
    if (dis < 1e-3)
      continue;
    TrajectoryPoint trajectory_point;
    trajectory_point.x = pt(0);
    trajectory_point.y = pt(1);
    trajectory_point.s = last_length + dis;
    trajectory_point.time = t;
    last_pt = pt;
    last_length = trajectory_point.s;
    trajectory.emplace_back(trajectory_point);
  }
  getPathKappa(trajectory);
  {
    for (size_t index = 1; index < trajectory.size(); ++index) {
      auto first = trajectory.at(index - 1);
      auto second = trajectory.at(index);
      double theta = std::atan2(second.y - first.y, second.x - first.x);
      trajectory.at(index - 1).yaw = theta;
    }
    if (trajectory.size() >= 2) {
      trajectory.back().yaw = trajectory.at(trajectory.size() - 2).yaw;
    }
  }
  return trajectory;
}

void Tools::getPathKappa(Trajectory &trajectory) {

  if (trajectory.size() <= 3) {
    std::cout << " can not get kappa with the size of the path points: "
              << trajectory.size() << std::endl;
    return;
  }
  std::vector<double> delt_x_over_delt_s;
  std::vector<double> delt_y_over_delt_s;
  std::vector<double> delt_x_over_delt_s2;
  std::vector<double> delt_y_over_delt_s2;
  delt_x_over_delt_s.reserve(trajectory.size());
  delt_y_over_delt_s.reserve(trajectory.size());
  delt_x_over_delt_s2.reserve(trajectory.size());
  delt_y_over_delt_s2.reserve(trajectory.size());

  for (size_t i = 0; i < trajectory.size(); ++i) {
    double x_ds = 0.0;
    double y_ds = 0.0;
    if (i == 0) {
      x_ds = (trajectory[i + 1].x - trajectory[i].x) /
             (trajectory[i + 1].s - trajectory[i].s);
      y_ds = (trajectory[i + 1].y - trajectory[i].y) /
             (trajectory[i + 1].s - trajectory[i].s);
    } else if (i == trajectory.size() - 1) {
      x_ds = (trajectory[i].x - trajectory[i - 1].x) /
             (trajectory[i].s - trajectory[i - 1].s);
      y_ds = (trajectory[i].y - trajectory[i - 1].y) /
             (trajectory[i].s - trajectory[i - 1].s);
    } else {
      x_ds = (trajectory[i + 1].x - trajectory[i - 1].x) /
             (trajectory[i + 1].s - trajectory[i - 1].s);
      y_ds = (trajectory[i + 1].y - trajectory[i - 1].y) /
             (trajectory[i + 1].s - trajectory[i - 1].s);
    }
    delt_x_over_delt_s.emplace_back(x_ds);
    delt_y_over_delt_s.emplace_back(y_ds);
  }
  for (size_t i = 0; i < trajectory.size(); ++i) {
    double x_dds = 0.0;
    double y_dds = 0.0;
    if (i == 0) {
      x_dds = (delt_x_over_delt_s[i + 1] - delt_x_over_delt_s[i]) /
              (trajectory[i + 1].s - trajectory[i].s);
      y_dds = (delt_y_over_delt_s[i + 1] - delt_y_over_delt_s[i]) /
              (trajectory[i + 1].s - trajectory[i].s);
    } else if (i == trajectory.size() - 1) {
      x_dds = (delt_x_over_delt_s[i] - delt_x_over_delt_s[i - 1]) /
              (trajectory[i].s - trajectory[i - 1].s + 1e-3);
      y_dds = (delt_y_over_delt_s[i] - delt_y_over_delt_s[i - 1]) /
              (trajectory[i].s - trajectory[i - 1].s + 1e-3);
    } else {
      x_dds = (delt_x_over_delt_s[i + 1] - delt_x_over_delt_s[i - 1]) /
              (trajectory[i + 1].s - trajectory[i - 1].s);
      y_dds = (delt_y_over_delt_s[i + 1] - delt_y_over_delt_s[i - 1]) /
              (trajectory[i + 1].s - trajectory[i - 1].s);
    }
    delt_x_over_delt_s2.emplace_back(x_dds);
    delt_y_over_delt_s2.emplace_back(y_dds);
  }

  for (size_t i = 0; i < trajectory.size(); ++i) {
    double xds = delt_x_over_delt_s[i];
    double yds = delt_y_over_delt_s[i];
    double xdds = delt_x_over_delt_s2[i];
    double ydds = delt_y_over_delt_s2[i];
    double num = std::sqrt((xds * xds + yds * yds)) * (xds * xds + yds * yds);
    double kappa = (xds * ydds - yds * xdds) / (num + 1e-9);
    trajectory.at(i).curve = kappa;
    if (i == 1)
      trajectory.at(i - 1).curve = kappa;
    if (i == trajectory.size() - 1)
      trajectory.at(i).curve = trajectory.at(i - 1).curve;
  }
}
void Tools::generateLissajousCurvePath(Trajectory &liss_trajectory) {

  Trajectory trajectory;
  double last_x, last_y;
  size_t N = 300 * 50;
  for (size_t i = 0; i < N; ++i) {
    double t = 300 * i / (N - 1.0);
    double x, y, vx, vy, ax, ay;
    double A = 10;
    double B = 2;
    double C = 50;
    x = A * (std::cos(M_PI / 4) - std::cos(6 * M_PI / C * t + M_PI / 4));
    y = B * (1 - std::cos(4 * M_PI / C * t));
    vx = A * 6 * M_PI / C * std::sin(6 * M_PI / C * t + M_PI / 4);
    vy = B * 4 * M_PI / C * sin(4 * M_PI / C * t);
    ax =
        A * 6 * M_PI / C * 6 * M_PI / C * std::cos(6 * M_PI / C * t + M_PI / 4);
    ay = B * 4 * M_PI / C * 4 * M_PI / C * sin(4 * M_PI / C * t);
    TrajectoryPoint trajectory_point;
    trajectory_point.x = x;
    trajectory_point.y = y;
    trajectory_point.yaw = std::atan2(vy, vx);
    trajectory_point.v = std::sqrt(vx * vx + vy * vy);
    trajectory_point.w = 0;
    trajectory_point.s = 0;
    if (i > 0)
      trajectory_point.s =
          trajectory.back().s +
          std::sqrt((x - last_x) * (x - last_x) + (y - last_y) * (y - last_y));
    trajectory_point.time = 0;
    trajectory_point.curve =
        std::fabs((vx * ay - vy * ax) / std::pow(vx * vx + vy * vy, 1.5));

    trajectory.emplace_back(trajectory_point);
    last_x = x;
    last_y = y;
  }
  liss_trajectory = trajectory;
}
} // namespace TOPP_DWR