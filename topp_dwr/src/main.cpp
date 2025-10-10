/************************************************************************
 * Software License Agreement (BSD License)
 *
 * @file main.cpp
 * @author liyong2018@zju.edu.cn
 * @brief test code for TOPP-DWR
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
#include "data_type.hpp"
#include "geometry_msgs/msg/pose_stamped.hpp"
#include "geometry_msgs/msg/pose_with_covariance_stamped.hpp" // subscribe /initialpose
#include "rclcpp/rclcpp.hpp"
#include "tools.hpp"
#include "topp_dwr.hpp"
#include "visualization_msgs/msg/marker_array.hpp"
#include <cmath>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

using namespace TOPP_DWR;

// --- Added: Global Variables and Callback Function ---// Global buffer: Stores
// points clicked by the user via 2D Pose Estimate
std::vector<TOPP_DWR::TrajectoryPoint> g_rviz_click_points;
// Mutex: Protects g_rviz_click_points to avoid race conditions between
// subscription callback and main thread
std::mutex g_click_points_mutex;
// Trigger flag: Set to true after the user clicks 2D Goal Pose to trigger
// trajectory generation
bool g_goal_triggered = false;

void publishTrajectory(
    const rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr
        publisher,
    const Trajectory &trajectory) {
  if (trajectory.empty() || publisher == nullptr) {
    std::cout << "Empty trajectory or empty publisher" << std::endl;
    return;
  }
  std::cout << "global_trajectory_size: " << trajectory.size() << std::endl;
  visualization_msgs::msg::MarkerArray trajectory_color_list;
  {
    // clear old trajectory in rviz
    visualization_msgs::msg::Marker trajectory_marker;
    trajectory_marker.id = 0;
    trajectory_marker.ns = "trajectory";
    trajectory_marker.type = visualization_msgs::msg::Marker::CUBE;
    trajectory_marker.header.frame_id = "map";
    trajectory_marker.scale.x = 0.05;
    trajectory_marker.scale.y = 0.05;
    trajectory_marker.scale.z = 0.05;
    auto duration = builtin_interfaces::msg::Duration();
    duration.set__sec(1e6);
    trajectory_marker.lifetime = duration;
    trajectory_marker.action = visualization_msgs::msg::Marker::DELETEALL;
    geometry_msgs::msg::PoseStamped temp_point;
    temp_point.pose.position.x = 0;
    temp_point.pose.position.y = 0;
    temp_point.pose.position.z = 0.15;
    trajectory_marker.pose.position = temp_point.pose.position;
    trajectory_marker.pose.orientation.w = 1.0;
    trajectory_color_list.markers.emplace_back(trajectory_marker);
  }

  double max_vel = 0;
  double min_vel = 2.0;
  for (unsigned int i = 0; i < trajectory.size(); ++i) {
    auto vel = trajectory.at(i).v;
    if (vel > max_vel)
      max_vel = vel;
    if (vel < min_vel)
      min_vel = vel;
  }
  std::cout << "min_vel: " << min_vel << " max_vel: " << max_vel << std::endl;
  if (max_vel - min_vel < 1e-3) {
    min_vel = 0;
    max_vel = 0.6;
  }
  int jump_num = 1;
  for (unsigned int index = 0; index < trajectory.size();
       index = index + jump_num) {
    auto vel = trajectory.at(index).v;
    std_msgs::msg::ColorRGBA color;
    color.r = 139.0 / 255.0;
    color.g = 1.0 - (std::fabs(max_vel - std::fabs(vel))) / (max_vel - min_vel);
    color.b = (std::fabs(max_vel - std::fabs(vel))) / (max_vel - min_vel) /
              255.0 * 50.0;
    bool visualize = true;
    if (visualize)
      color.a = 1.0;
    else
      color.a = 0;
    Quaternion q_yaw = yawToQuaternion(0, 0, trajectory.at(index).yaw);
    visualization_msgs::msg::Marker trajectory_marker;
    trajectory_marker.id = index;
    trajectory_marker.ns = "trajectory";
    trajectory_marker.type = visualization_msgs::msg::Marker::CUBE;
    trajectory_marker.header.frame_id = "map";
    trajectory_marker.scale.x = 0.05;
    trajectory_marker.scale.y = 0.05;
    trajectory_marker.scale.z = 0.05;
    trajectory_marker.color = color;
    auto duration = builtin_interfaces::msg::Duration();
    duration.set__sec(1e6);
    trajectory_marker.lifetime = duration;
    trajectory_marker.action = visualization_msgs::msg::Marker::ADD;
    geometry_msgs::msg::PoseStamped temp_point;
    temp_point.pose.position.x = trajectory.at(index).x;
    temp_point.pose.position.y = trajectory.at(index).y;
    temp_point.pose.position.z = 0.15;
    trajectory_marker.pose.position = temp_point.pose.position;
    trajectory_marker.pose.orientation.x = q_yaw.x;
    trajectory_marker.pose.orientation.y = q_yaw.y;
    trajectory_marker.pose.orientation.z = q_yaw.z;
    trajectory_marker.pose.orientation.w = q_yaw.w;

    trajectory_color_list.markers.emplace_back(trajectory_marker);
  }
  try {
    publisher->publish(trajectory_color_list);
  } catch (std::exception &e) {
    std::cout << "visualization_msgs::msg::MarkerArray pub exception";
  }
}
/**
    @brief Handles RViz 2D Pose Estimate clicks (/initialpose topic)
    Buffers the clicked points and performs distance checks with real-time
   prompts.
*/
void initialpose_callback(
    const geometry_msgs::msg::PoseWithCovarianceStamped::SharedPtr msg) {
  std::lock_guard<std::mutex> lock(g_click_points_mutex);

  // Check if enough points have been collected to prevent duplicate clicks
  if (g_rviz_click_points.size() >= 5) {
    std::cout << "[Warn] Already collected 5 points. Please click '2D Goal "
                 "Pose' to generate trajectory."
              << std::endl;
  }
  TOPP_DWR::TrajectoryPoint new_point;
  new_point.x = msg->pose.pose.position.x;
  new_point.y = msg->pose.pose.position.y;
  double qx = msg->pose.pose.orientation.x;
  double qy = msg->pose.pose.orientation.y;
  double qz = msg->pose.pose.orientation.z;
  double qw = msg->pose.pose.orientation.w;
  new_point.yaw = atan2(2 * (qw * qz + qx * qy), 1 - 2 * (qy * qy + qz * qz));
  new_point.v = 0.0;

  // --- distance check ---
  const double MIN_DISTANCE = 0.2;
  if (!g_rviz_click_points.empty()) {
    const auto &last_point = g_rviz_click_points.back();
    double dx = new_point.x - last_point.x;
    double dy = new_point.y - last_point.y;
    double distance = std::sqrt(dx * dx + dy * dy);

    if (distance < MIN_DISTANCE) {
      std::cerr
          << "[Error] New point is too close to the previous one! (Distance: "
          << distance << "m, Min required: " << MIN_DISTANCE << "m)"
          << std::endl;
      return;
    }
  }
  g_rviz_click_points.push_back(new_point);
  int collected = g_rviz_click_points.size();
  int remaining = 5 - collected;
  std::cout << "[Info] Collected point " << collected << "/5. "
            << (remaining > 0
                    ? "Need " + std::to_string(remaining) + " more point(s)."
                    : "Ready to generate!")
            << std::endl;
}

/**
 * @brief handle RViz 2D Goal Pose Click（/goal_pose）
 */
void goalpose_callback(const geometry_msgs::msg::PoseStamped::SharedPtr msg) {
  std::lock_guard<std::mutex> lock(g_click_points_mutex);
  if (g_rviz_click_points.size() < 5) {
    std::cerr << "[Error] Not enough 2D Pose Estimate points! Need 5, current "
              << g_rviz_click_points.size() << std::endl;
    return;
  }
  g_goal_triggered = true;
  std::cout
      << "[Info] Received 2D Goal Pose! Triggering trajectory generation..."
      << std::endl;
}

bool g_should_exit = false;
void signal_handler(int signal_num) {
  if (signal_num == SIGINT || signal_num == SIGTERM) {
    std::cout << "\nReceived Ctrl+C . Exiting gracefully..." << std::endl;
    g_should_exit = true;

    if (rclcpp::ok()) {
      rclcpp::shutdown();
      std::cout << "ROS 2 node shutdown successfully." << std::endl;
    }
  }
}

void print_help() {
  std::cout << "Usage: ./topp_dwr_node [OPTIONS]" << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "  --algo [type]    Select algorithm (required):" << std::endl;
  std::cout << "                    mosek_c    : Use mosek 10.2 C API"
            << std::endl;
  std::cout << "                    mosek_cpp  : Use mosek 10.2 C++ API"
            << std::endl;
  std::cout << "                    ecos       : Use ecos API" << std::endl;
  std::cout << "  --traj [type]     Select trajectory type (required):"
            << std::endl;
  std::cout << "                    c1         : Standard Lissajous curve"
            << std::endl;
  std::cout << "                    c2         : Complex CVTE curve"
            << std::endl;
  std::cout << "                    c3         : Full-coverage path"
            << std::endl;
  std::cout << "                    c4         : Generate path from RViz "
               "clicks (2D Pose Estimate x5 + 2D Goal Pose)"
            << std::endl;
  std::cout << "  --help            Show this help message" << std::endl;
  std::cout << "Example:" << std::endl;
  std::cout << "  ./install/topp_dwr/lib/topp_dwr/topp_dwr_node --algo mosek_c "
               "--traj c1"
            << std::endl;
  std::cout << "  ./install/topp_dwr/lib/topp_dwr/topp_dwr_node --algo ecos "
               "--traj c2"
            << std::endl;
  std::cout << "  ./install/topp_dwr/lib/topp_dwr/topp_dwr_node --algo "
               "mosek_cpp --traj c4"
            << std::endl;
}

int main(int argc, char **argv) {
  signal(SIGINT, signal_handler);

  std::string selected_algo = "";
  std::string selected_traj = "";
  const std::vector<std::string> VALID_ALGOS = {"mosek_c", "mosek_cpp", "ecos"};
  const std::vector<std::string> VALID_TRAJS = {"c1", "c2", "c3", "c4"};

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--help") {
      print_help();
      return 0;
    } else if (arg == "--algo" && i + 1 < argc) {
      selected_algo = argv[++i];
      if (std::find(VALID_ALGOS.begin(), VALID_ALGOS.end(), selected_algo) ==
          VALID_ALGOS.end()) {
        std::cerr << "[Error] Invalid algorithm: " << selected_algo
                  << std::endl;
        print_help();
        return -1;
      }
    } else if (arg == "--traj" && i + 1 < argc) {
      selected_traj = argv[++i];
      if (std::find(VALID_TRAJS.begin(), VALID_TRAJS.end(), selected_traj) ==
          VALID_TRAJS.end()) {
        std::cerr << "[Error] Invalid trajectory type: " << selected_traj
                  << std::endl;
        print_help();
        return -1;
      }
    } else {
      std::cerr << "[Error] Unknown argument: " << arg << std::endl;
      print_help();
      return -1;
    }
  }

  if (selected_algo.empty() || selected_traj.empty()) {
    std::cerr << "[Error] Missing required arguments (--algo and --traj must "
                 "be specified)!"
              << std::endl;
    print_help();
    return -1;
  }

  rclcpp::init(argc, argv);
  auto vis_node = rclcpp::Node::make_shared("topp_dwr");
  auto marker_publisher =
      vis_node->create_publisher<visualization_msgs::msg::MarkerArray>(
          "trajectory", 10);

  std::cout << "[Info] TOPP-DWR node started. Selected: Algo=" << selected_algo
            << ", Traj=" << selected_traj << std::endl;

  TOPP_DWR::Tools tools;
  std::vector<TOPP_DWR::TrajectoryPoint> trajectory;
  trajectory.clear();

  try {
    if (selected_traj == "c1") {
      tools.generateLissajousCurvePath(trajectory);
    } else if (selected_traj == "c2") {
      std::cout << "[Info] Loading CVTE complex curve..." << std::endl;
      std::string file = "cvte";
      tools.loadCurve(file, trajectory);
    } else if (selected_traj == "c3") {
      std::cout << "[Info] Loading full-coverage path..." << std::endl;
      std::string file = "ccpp";
      tools.loadCurve(file, trajectory);
    } else if (selected_traj == "c4") {
      std::cout << "[Info] Entering RViz interactive mode:" << std::endl;
      std::cout << "  1. Click '2D Pose Estimate' in RViz at least 5 times."
                << std::endl;
      std::cout << "     - Minimum distance between consecutive points: 0.2m."
                << std::endl;
      std::cout << "  2. Click '2D Goal Pose' in RViz to trigger trajectory "
                   "generation."
                << std::endl;
      auto initialpose_sub = vis_node->create_subscription<
          geometry_msgs::msg::PoseWithCovarianceStamped>("/initialpose", 10,
                                                         initialpose_callback);
      auto goalpose_sub =
          vis_node->create_subscription<geometry_msgs::msg::PoseStamped>(
              "/goal_pose", 10, goalpose_callback);

      // Loop and wait until the user clicks the Goal or exits with Ctrl+C
      while (!g_goal_triggered && !g_should_exit) {
        rclcpp::spin_some(vis_node); // Process subscription callbacks
        std::this_thread::sleep_for(
            std::chrono::milliseconds(100)); // Reduce CPU usage
      }

      // Check whether it was triggered normally or interrupted
      if (g_should_exit) {
        throw std::runtime_error("Program interrupted by user.");
      }

      // Copy the buffered points to the final trajectory
      {
        std::lock_guard<std::mutex> lock(g_click_points_mutex);
        trajectory = g_rviz_click_points;
      }

      // Call the tool function to smooth/interpolate the original clicked
      // points
      std::cout << "[Info] Generating trajectory from " << trajectory.size()
                << " RViz clicks..." << std::endl;
      tools.generateRVIZPath(trajectory);
      std::cout << "[Info] Trajectory generated with " << trajectory.size()
                << " points." << std::endl;
    }

  } catch (const std::exception &e) {
    std::cerr << "[Error] Failed to generate trajectory: " << e.what()
              << std::endl;
    rclcpp::shutdown();
    return -1;
  }

  if (trajectory.empty()) {
    std::cerr << "[Error] Generated trajectory is empty! Check trajectory "
                 "function implementation."
              << std::endl;
    rclcpp::shutdown();
    return -1;
  }
  std::cout << "[Info] Trajectory generated successfully, total points: "
            << trajectory.size() << std::endl;

  TOPP_DWR::TOPP_DWR_Params params;
  if (!TOPP_DWR::Tools::readTOPP_DWR_ParamsFromPackage(selected_traj, params)) {
    std::cerr << "[Warn] Use default params (failed to read config file)!"
              << std::endl;
  }

  TOPP_DWR::TOPPDWR topp_dwr;

  try {
    if (selected_algo == "mosek_c") {
      std::cout << "[Info] Running solveWithMosekC..." << std::endl;
      topp_dwr.solveWithMosekC(trajectory, params);
    } else if (selected_algo == "mosek_cpp") {
      std::cout << "[Info] Running solveWithMosekCPP..." << std::endl;
      topp_dwr.solveWithMosekCPP(trajectory, params);
    } else if (selected_algo == "ecos") {
      std::cout << "[Info] Running solveWithEcos..." << std::endl;
      topp_dwr.solveWithEcos(trajectory, params);
    }
  } catch (const std::exception &e) {
    std::cerr << "[Error] Algorithm execution failed: " << e.what()
              << std::endl;
    rclcpp::shutdown();
    return -1;
  }
  std::cout << "[Info] Algorithm finished successfully." << std::endl;
  std::cout << "[Info] Starting trajectory publication..." << std::endl;
  for (size_t i = 0; i < 5 && !g_should_exit; ++i) {
    publishTrajectory(marker_publisher, trajectory);
    std::this_thread::sleep_for(std::chrono::seconds(1));
  }

  std::cout << "[Info] Program exiting. Cleaning up ROS 2 resources..."
            << std::endl;
  if (rclcpp::ok()) {
    rclcpp::shutdown();
  }

  return 0;
}
