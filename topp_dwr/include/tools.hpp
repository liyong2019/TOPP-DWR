/************************************************************************
 * Software License Agreement (BSD License)
 *
 * @file tools.hpp
 * @author liyong2018@zju.edu.cn
 * @brief tools for TOPP-WDR
 * @version  1.0.0
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
#pragma once
#include "data_type.hpp"
#include "non_uniform_bspline.hpp"
#include "yaml-cpp/yaml.h"
#include <vector>
namespace TOPP_DWR {
class Tools {
private:
  static bool parsePoint(const YAML::Node &point, Trajectory &raw_trajectory);
  static bool getPathFromFile(const std::string &filepath,
                              Trajectory &reference_path);

public:
  Tools(){};
  ~Tools(){};
  void generateLissajousCurvePath(Trajectory &path);
  bool generateRVIZPath(Trajectory &path);
  bool loadCurve(const std::string file, Trajectory &trajectory);

  static bool readTOPP_DWR_Params(const std::string &param_file_path,
                                  const std::string &traj_type,
                                  TOPP_DWR_Params &out_params);

  static bool readTOPP_DWR_ParamsFromPackage(const std::string &traj_type,
                                             TOPP_DWR_Params &out_params);
  Trajectory getRawPathWithBspline(NonUniformBspline &bspline, double delta_t);
  void getPathKappa(Trajectory &trajectory);
};
} // namespace TOPP_DWR
