/************************************************************************
 * Software License Agreement (BSD License)
 *
 * @file TOPPDWR.hpp
 * @author liyong2018@zju.edu.cn
 * @brief open-source code for IROS2025 paper: TOPP-DWR：Time-Optimal Path
 *        Parameterization of Differential-Driven Wheeled Robots Considering
 *        Piecewise-Constant Angular Velocity Constraints
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
#include "tools.hpp"
namespace TOPP_DWR {
class TOPPDWR {
private:
public:
  TOPPDWR(){};
  ~TOPPDWR(){};
  void solveWithMosekC(Trajectory &trajectory, const TOPP_DWR_Params &params);
  void solveWithMosekCPP(Trajectory &trajectory, const TOPP_DWR_Params &params);
  void solveWithEcos(Trajectory &trajectory, const TOPP_DWR_Params &params);
};
} // namespace TOPP_DWR