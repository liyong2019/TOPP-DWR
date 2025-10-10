/************************************************************************
 * Software License Agreement (BSD License)
 *
 * @file topp_dwr.cpp
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
#include "topp_dwr.hpp"
#include "chrono"
#include "ecos.h"
#include "fusion.h"
#include <cmath>
#include <iostream>
using namespace mosek::fusion;
using namespace monty;
bool PRINT_ECOS = false;
namespace TOPP_DWR {
static void MSKAPI printstr(void *handle, const char str[]) {
  (void)handle;
  std::cout << str;
}
void TOPPDWR::solveWithMosekC(Trajectory &trajectory,
                              const TOPP_DWR_Params &params) {
  auto time_before_init = std::chrono::system_clock::now();
  std::cout << "path length: " << trajectory.back().s << " size "
            << trajectory.size() << std::endl;
  double v_max = params.v_max;
  double a_max_l = params.a_max_l;
  double a_max_n = params.a_max_n;
  double vr_max = params.vr_max;
  double vl_max = params.vl_max;
  double d = params.d;
  double w_max = params.w_max;
  double v_0 = std::min(params.start_vel, v_max);
  double v_s = std::min(params.end_vel, v_max);

  std::vector<double> v_bar;
  std::vector<double> dis;
  std::vector<double> theta;
  std::vector<double> g;
  std::vector<double> h_theta;
  size_t N = trajectory.size();
  v_bar.reserve(N);
  dis.reserve(N - 1);
  theta.reserve(N - 1);
  g.reserve(N - 1);
  h_theta.reserve(N - 1);

  double kappa = std::fabs(trajectory.front().curve);
  double v_bar_0 = std::min(v_max, std::sqrt(a_max_n / kappa));
  v_bar.emplace_back(v_bar_0);
  for (size_t i = 1; i < N; ++i) {
    auto first = trajectory.at(i - 1);
    auto second = trajectory.at(i);
    double dis_i_j = first.distanceTo(second);
    double theta_i_j = first.angle_diff(second);
    double kappa = std::fabs(trajectory.at(i).curve);
    double v_bar_i = std::min(v_max, std::sqrt(a_max_n / kappa));
    v_bar.emplace_back(v_bar_i);
    dis.emplace_back(dis_i_j);
    theta.emplace_back(theta_i_j);
    g.emplace_back(d * theta_i_j / (4.0 * dis_i_j));
    h_theta.emplace_back(
        std::min(2.0 * dis_i_j / std::fabs(theta_i_j) * w_max, 1e6));
  }
  if (v_0 > v_bar.front() + 1e-3)
    std::cout << " start vel error!, v_0: " << v_0 << " v_bar " << v_bar.front()
              << std::endl;
  if (v_s > v_bar.back() + 1e-3)
    std::cout << " end vel error!, v_s: " << v_s << " v_bar " << v_bar.back()
              << std::endl;
  auto time_after_init = std::chrono::system_clock::now();
  std::chrono::duration<double> time_init = time_after_init - time_before_init;
  std::cout << "TIIME: trajectory plan initializtion: " << time_init.count()
            << std::endl;
  {
    auto time_before_construct = std::chrono::system_clock::now();
    MSKrescodee r;
    const MSKint32t numv = N;
    const MSKint32t numy = N - 1;
    // const MSKint32t numt = 0;
    const MSKint32t numvar = numv + numy; // v_i->N, y_i->N-1,std::sqrt(2)->N-1
    const MSKint32t numacc = N - 1;
    const MSKint32t numdec = N - 1;
    const MSKint32t numw = N - 1;
    const MSKint32t numvr = N - 1;
    const MSKint32t numvl = N - 1;

    const MSKint32t numcon = numacc + numdec + numw + numvr + numvl;
    const MSKint32t numf =
        3 * numacc + 3 * numdec + 2 * numw + 2 * numvr + 2 * numvl;
    const MSKint64t numafe = 3 * numy;
    const MSKint64t numafcc = numy;
    const MSKint64t f_nnz = 3 * numy;

    // 1.optimization objuective
    std::vector<double> c;
    c.reserve((size_t)numvar);
    c.resize((size_t)numvar);
    // double c[numvar];
    for (size_t i = 0; i < (size_t)numvar; ++i) {
      if (i < N)
        c[i] = 0.0;
      else if (i < 2 * N - 1) {
        c[i] = 2.0 * dis.at(i - N);
      } else
        c[i] = 0.0;
    }
    // 2. variable constraints(v_i,y_i,std::sqrt(2))
    std::vector<MSKboundkeye> bkx;
    bkx.reserve((size_t)numvar);
    bkx.resize((size_t)numvar);
    std::vector<MSKrealt> blx;
    blx.reserve((size_t)numvar);
    blx.resize((size_t)numvar);
    std::vector<MSKrealt> bux;
    bux.reserve((size_t)numvar);
    bux.resize((size_t)numvar);

    for (size_t i = 0; i < N; ++i) {
      // v_i
      if (i == 0) {
        bkx[i] = MSK_BK_FX;
        blx[i] = v_0;
        bux[i] = v_0;
      } else if (i == N - 1) {
        bkx[i] = MSK_BK_FX;
        blx[i] = v_s;
        bux[i] = v_s;
      } else {
        bkx[i] = MSK_BK_RA;
        blx[i] = 0.0;
        bux[i] = v_bar.at(i);
      }
    }
    for (size_t i = N; i < 2 * N - 1; ++i) {
      // y_i
      bkx[i] = MSK_BK_LO;
      blx[i] = 1.0 / (2.0 * v_max); //
      // blx[i] = 0.0;
      bux[i] = +MSK_INFINITY;
    }

    // 3.1 left part for linear constraints
    std::vector<MSKint32t> aptrb;
    std::vector<MSKint32t> aptre;
    std::shared_ptr<MSKint32t[]> asub =
        std::shared_ptr<MSKint32t[]>(new MSKint32t[numf]);
    std::shared_ptr<double[]> aval =
        std::shared_ptr<double[]>(new double[numf]);
    aptrb.reserve((size_t)numcon);
    aptrb.resize((size_t)numcon);
    aptre.reserve((size_t)numcon);
    aptre.resize((size_t)numcon);

    // ub*y_i>=v_(i+1)-v_i
    for (MSKint32t i = 0; i < numacc; ++i) {
      aval[3 * i] = -1.0;
      aval[3 * i + 1] = 1.0;
      aval[3 * i + 2] = -2.0 * dis.at(i) * a_max_l;

      asub[3 * i] = i;
      asub[3 * i + 1] = i + 1;
      asub[3 * i + 2] = N + i;

      aptrb[i] = 3 * i;
      aptre[i] = 3 * i + 3;
    }
    std::cout << " nnnum numacc: " << numacc << std::endl;
    // ub*y_i>=v_i-v_(i+1)
    for (MSKint32t i = 0; i < numdec; ++i) {
      aval[3 * numacc + 3 * i] = 1.0;
      aval[3 * numacc + 3 * i + 1] = -1.0;
      aval[3 * numacc + 3 * i + 2] = -2.0 * dis.at(i) * a_max_l;

      asub[3 * numacc + 3 * i] = i;
      asub[3 * numacc + 3 * i + 1] = i + 1;
      asub[3 * numacc + 3 * i + 2] = N + i;

      aptrb[numacc + i] = 3 * numacc + 3 * i;
      aptre[numacc + i] = 3 * numacc + 3 * i + 3;
    }
    std::cout << " nnnum numacc+numdec: " << numacc + numdec << std::endl;
    // 0=<v_i+v_(i+1)<=w_max
    for (MSKint32t i = 0; i < numw; ++i) {
      aval[3 * numacc + 3 * numdec + 2 * i] = 1.0;
      aval[3 * numacc + 3 * numdec + 2 * i + 1] = 1.0;

      asub[3 * numacc + 3 * numdec + 2 * i] = i;
      asub[3 * numacc + 3 * numdec + 2 * i + 1] = i + 1;

      aptrb[numacc + numdec + i] = 3 * numacc + 3 * numdec + 2 * i;
      aptre[numacc + numdec + i] = 3 * numacc + 3 * numdec + 2 * i + 2;
    }
    std::cout << " nnnum numacc + numdec + numw: " << numacc + numdec + numw
              << std::endl;

    // -v_max=<(1+g)*v_i+g*v_(i+1)<=v_max
    for (MSKint32t i = 0; i < numvr; ++i) {
      aval[3 * numacc + 3 * numdec + 2 * numw + 2 * i] = 1.0 + g.at(i);
      aval[3 * numacc + 3 * numdec + 2 * numw + 2 * i + 1] = g.at(i);

      asub[3 * numacc + 3 * numdec + 2 * numw + 2 * i] = i;
      asub[3 * numacc + 3 * numdec + 2 * numw + 2 * i + 1] = i + 1;

      aptrb[numacc + numdec + numw + i] =
          3 * numacc + 3 * numdec + 2 * numw + 2 * i;
      aptre[numacc + numdec + numw + i] =
          3 * numacc + 3 * numdec + 2 * numw + 2 * i + 2;
    }
    std::cout << " nnnum numacc + numdec + numw+numvr: "
              << numacc + numdec + numw + numvr << std::endl;
    // -v_max=<(1-g)*v_i-g*v_(i+1)<=v_max
    for (MSKint32t i = 0; i < numvl; ++i) {
      aval[3 * numacc + 3 * numdec + 2 * numw + 2 * numvr + 2 * i] =
          1.0 - g.at(i);
      aval[3 * numacc + 3 * numdec + 2 * numw + 2 * numvr + 2 * i + 1] =
          -g.at(i);

      asub[3 * numacc + 3 * numdec + 2 * numw + 2 * numvr + 2 * i] = i;
      asub[3 * numacc + 3 * numdec + 2 * numw + 2 * numvr + 2 * i + 1] = i + 1;

      aptrb[numacc + numdec + numw + numvr + i] =
          3 * numacc + 3 * numdec + 2 * numw + 2 * numvr + 2 * i;
      aptre[numacc + numdec + numw + numvr + i] =
          3 * numacc + 3 * numdec + 2 * numw + 2 * numvr + 2 * i + 2;
    }
    std::cout << " nnnum numacc + numdec + numw + numvr + numvl: "
              << numacc + numdec + numw + numvr + numvl << std::endl;

    // 3.2 right part for linear constraints
    std::vector<MSKboundkeye> bkc;
    std::vector<double> blc;
    std::vector<double> buc;
    bkc.reserve((size_t)numcon);
    bkc.resize((size_t)numcon);
    blc.reserve((size_t)numcon);
    blc.resize((size_t)numcon);
    buc.reserve((size_t)numcon);
    buc.resize((size_t)numcon);

    for (MSKint32t i = 0; i < numcon; ++i) {
      if (i < numacc + numdec) {
        bkc[i] = MSK_BK_UP;
        blc[i] = -MSK_INFINITY;
        buc[i] = 0.0;
      } else if (i < numacc + numdec + numw) {
        bkc[i] = MSK_BK_RA;
        blc[i] = 0.0;
        buc[i] = h_theta.at(i - (numacc + numdec));
      } else {
        bkc[i] = MSK_BK_UP;
        blc[i] = -MSK_INFINITY;
        buc[i] = vr_max;
      }
    }

    std::shared_ptr<MSKrealt[]> f_val =
        std::shared_ptr<MSKrealt[]>(new MSKrealt[3 * (N - 1)]); // value
    std::shared_ptr<MSKint64t[]> afeidx =
        std::shared_ptr<MSKint64t[]>(new MSKint64t[3 * (N - 1)]); // index i
    std::shared_ptr<MSKint32t[]> varidx =
        std::shared_ptr<MSKint32t[]>(new MSKint32t[3 * (N - 1)]); // index j

    MSKrealt g = std::sqrt(2.0); // g value
    for (size_t i = 0; i < N - 1; ++i) {
      f_val[3 * i] = 1.0;
      afeidx[3 * i] = 3 * i;
      varidx[3 * i] = N + i;

      f_val[3 * i + 1] = 1.0;
      afeidx[3 * i + 1] = 3 * i + 1;
      varidx[3 * i + 1] = i;

      f_val[3 * i + 2] = 1.0;
      afeidx[3 * i + 2] = 3 * i + 1;
      varidx[3 * i + 2] = i + 1;
    }

    std::shared_ptr<MSKint64t[]> domidx =
        std::shared_ptr<MSKint64t[]>(new MSKint64t[numafcc]); // value
    for (int i = 0; i < numafcc; ++i) {
      domidx[i] = i;
    }
    MSKint32t i, j;
    MSKenv_t env = NULL;
    MSKtask_t task = NULL;

    /* Create the mosek environment. */
    r = MSK_makeenv(&env, NULL);

    if (r == MSK_RES_OK) {
      /* Create the optimization task. */
      r = MSK_maketask(env, numcon, numvar, &task);
      // MSK_putintparam(task, MSK_IPAR_NUM_THREADS, 32);

      if (r == MSK_RES_OK) {
        // MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);
        MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr);

        /* Append 'numvar' variables.The variables will initially be fixed
         * at zero (x=0). */
        if (r == MSK_RES_OK)
          r = MSK_appendvars(task, numvar);
        if (r != MSK_RES_OK)
          std::cout << " append MSK_appendvars failed! " << std::endl;
        /* Append 'numafe' affine expressions.The affine expressions will
         * initially be empty. */
        if (r == MSK_RES_OK)
          r = MSK_appendafes(task, numafe);

        if (r != MSK_RES_OK)
          std::cout << " append MSK_appendafes failed! " << std::endl;

        for (j = 0; j < numvar && r == MSK_RES_OK; ++j) {
          /* Set the linear term c_j in the objective.*/
          if (r == MSK_RES_OK)
            r = MSK_putcj(task, j, c[j]);

          /* Set the bounds on variable j.
          blx[j] <= x_j <= bux[j] */
          if (r == MSK_RES_OK)
            r = MSK_putvarbound(task, j, /* Index of variable.*/
                                bkx[j],  /* Bound key.*/
                                blx[j],  /* Numerical value of lower bound.*/
                                bux[j]); /* Numericalvalue of upper bound.*/
        }
        if (r != MSK_RES_OK)
          std::cout << "MSK_putvarbound failed! " << std::endl;
        {
          // /* Append 'numcon' empty constraints.
          // The constraints will initially have no bounds. */
          if (r == MSK_RES_OK)
            r = MSK_appendcons(task, numcon);
          for (i = 0; i < numcon && r == MSK_RES_OK; ++i) {
            /* Input column j of A */
            if (r == MSK_RES_OK) {
              r = MSK_putarow(
                  task, i,               /* Variable (column) index.*/
                  aptre[i] - aptrb[i],   /* Number of non-zeros in column j.*/
                  asub.get() + aptrb[i], /* Pointer to row indexes of columnj.*/
                  aval.get() + aptrb[i]); /* Pointer to Values of column j.*/
            }
          }
          if (r != MSK_RES_OK)
            std::cout << " MSK_putarow failed! " << std::endl;
          // /* Set the bounds on constraints.
          //  for i=1, ...,numcon : blc[i] <= constraint i <= buc[i] */
          for (i = 0; i < numcon && r == MSK_RES_OK; ++i)
            r = MSK_putconbound(task, i, /* Index of constraint.*/
                                bkc[i],  /* Bound key.*/
                                blc[i],  /* Numerical value of lower bound.*/
                                buc[i]); /* Numerical value of  upper bound.*/
        }
        if (r != MSK_RES_OK)
          std::cout << " MSK_putconbound failed! " << std::endl;

        if (r == MSK_RES_OK) {
          /* Set the non-zero entries of the F matrix */
          r = MSK_putafefentrylist(task, f_nnz, afeidx.get(), varidx.get(),
                                   f_val.get());
          if (r != MSK_RES_OK)
            std::cout << "MSK_putafefentrylist failed! " << std::endl;
          for (MSKint32t jj = 0; jj < numafcc && r == MSK_RES_OK; ++jj)
            r = MSK_putafeg(task, 3 * jj + 2, g);

          /* Append rotated quadratic cone domain */
          if (r == MSK_RES_OK) {
            for (size_t j = 0; j < (size_t)numafcc && r == MSK_RES_OK; ++j) {
              r = MSK_appendrquadraticconedomain(task, 3, &domidx[j]);
            }
            if (r != MSK_RES_OK)
              std::cout << " append rq failed! " << std::endl;
          }

          /* Append two ACCs made up of the AFEs and the domains defined
          above.
           */
          if (r == MSK_RES_OK) {
            r = MSK_appendaccsseq(task, numafcc, domidx.get(), numafe,
                                  afeidx[0], NULL);
            if (r != MSK_RES_OK)
              std::cout << " append acc failed! " << std::endl;
          }
        }
        std::cout << " MOSEK CONSRTUCT success" << std::endl;

        auto time_before_solve = std::chrono::system_clock::now();

        std::chrono::duration<double> time_construct =
            time_before_solve - time_before_construct;
        std::cout << "TIIME: trajectory plan construct: "
                  << time_construct.count() << std::endl;
        if (r == MSK_RES_OK) {
          MSKrescodee trmcode;
          /* Run optimizer */
          r = MSK_optimizetrm(task, &trmcode);
          auto time_after_solve = std::chrono::system_clock::now();
          std::chrono::duration<double> time_solve =
              time_after_solve - time_before_solve;
          std::cout << "TIIME: trajectory plan mosek solve time: "
                    << time_solve.count() << std::endl;

          if (r == MSK_RES_OK) {
            MSKsolstae solsta;
            MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
            switch (solsta) {
            case MSK_SOL_STA_OPTIMAL: {
              std::cout << "MOSEK C SOLVE SUCCESS!" << std::endl;
              double *xx = NULL;

              xx = reinterpret_cast<double *>(calloc(numvar, sizeof(double)));
              if (xx) {
                MSK_getxx(task, MSK_SOL_ITR, xx);

                std::cout << "Optimal primal solution:" << std::endl;
                bool apply = true;
                if (apply) {
                  //计算时间戳
                  trajectory.front().time = 0;
                  trajectory.front().v = xx[0];
                  for (size_t i = 1; i < trajectory.size(); ++i) {
                    double theta_i = theta.at(i - 1);
                    double dis_i = dis.at(i - 1);
                    double time = 2.0 * dis_i / (xx[i] + xx[i - 1]);
                    trajectory.at(i).time = trajectory.at(i - 1).time + time;
                    trajectory.at(i).v = xx[i];
                    trajectory.at(i).w = theta_i / time;
                  }
                  std::cout
                      << "trajectory duration time: " << trajectory.back().time
                      << std::endl;
                }
              } else {
                r = MSK_RES_ERR_SPACE;
              }
              free(xx);
            } break;
            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
              std::cout << "Primal or dual infeasibility certificate found."
                        << std::endl;
              break;
            case MSK_SOL_STA_UNKNOWN:
              std::cout << "The status of the solution could not be "
                           "determined. "
                           "Termination code: "
                        << trmcode << std::endl;
              break;
            default:
              std::cout << "Other solution status." << std::endl;
              break;
            }
          } else {
            std::cout << "Error while optimizing." << std::endl;
          }
        }

        if (r != MSK_RES_OK) {
          /* In case of an error print error code and description. */
          char symname[MSK_MAX_STR_LEN];
          char desc[MSK_MAX_STR_LEN];

          std::cout << "An error occurred while optimizing." << std::endl;
          MSK_getcodedesc(r, symname, desc);
          std::cout << "Error symname: " << symname << " desc " << desc
                    << std::endl;
        }
      }
      /* Delete the task and the associated data. */
      MSK_deletetask(&task);
      MSK_deleteenv(&env);
    } else {
      /* Delete the environment and the associated data. */
      MSK_deletetask(&task);
      MSK_deleteenv(&env);
    }
  }
  std::cout << "mosek C solver finished!" << std::endl;
}
void TOPPDWR::solveWithMosekCPP(Trajectory &trajectory,
                                const TOPP_DWR_Params &params) {
  auto time_before_init = std::chrono::system_clock::now();
  std::cout << "path length: " << trajectory.back().s << " size "
            << trajectory.size() << std::endl;
  double v_max = params.v_max;
  double a_max_l = params.a_max_l;
  double a_max_n = params.a_max_n;
  double vr_max = params.vr_max;
  double vl_max = params.vl_max;
  double d = params.d;
  double w_max = params.w_max;
  double v_0 = std::min(params.start_vel, v_max);
  double v_s = std::min(params.end_vel, v_max);

  std::vector<double> v_bar;
  std::vector<double> dis;
  std::vector<double> theta;
  std::vector<double> g;
  std::vector<double> h_theta;
  size_t N = trajectory.size();
  v_bar.reserve(N);
  dis.reserve(N - 1);
  theta.reserve(N - 1);
  g.reserve(N - 1);
  h_theta.reserve(N - 1);

  double kappa = std::fabs(trajectory.front().curve);
  double v_bar_0 = std::min(v_max, std::sqrt(a_max_n / kappa));
  v_bar.emplace_back(v_bar_0);
  for (size_t i = 1; i < N; ++i) {
    auto first = trajectory.at(i - 1);
    auto second = trajectory.at(i);
    double dis_i_j = first.distanceTo(second);
    double theta_i_j = first.angle_diff(second);
    double kappa = std::fabs(trajectory.at(i).curve);
    double v_bar_i = std::min(v_max, std::sqrt(a_max_n / kappa));
    v_bar.emplace_back(v_bar_i);
    dis.emplace_back(dis_i_j);
    theta.emplace_back(theta_i_j);
    g.emplace_back(d * theta_i_j / (4.0 * dis_i_j));
    h_theta.emplace_back(
        std::min(2.0 * dis_i_j / std::fabs(theta_i_j) * w_max, 1e6));
  }
  if (v_0 > v_bar.front() + 1e-3)
    std::cout << " start vel error!, v_0: " << v_0 << " v_bar " << v_bar.front()
              << std::endl;
  if (v_s > v_bar.back() + 1e-3)
    std::cout << " end vel error!, v_s: " << v_s << " v_bar " << v_bar.back()
              << std::endl;
  auto time_after_init = std::chrono::system_clock::now();
  std::chrono::duration<double> time_init = time_after_init - time_before_init;
  std::cout << "TIIME: trajectory plan initializtion: " << time_init.count()
            << std::endl;
  {
    auto time_before_construct = std::chrono::system_clock::now();
    mosek::fusion::Model::t M = new mosek::fusion::Model("cqo1");

    auto _M = monty::finally([&]() { M->dispose(); });

    mosek::fusion::Variable::t v =
        M->variable("v", N, mosek::fusion::Domain::greaterThan(0.0));
    mosek::fusion::Variable::t y =
        M->variable("y", N - 1, mosek::fusion::Domain::greaterThan(0));
    M->constraint("cx_0", v->index(0), mosek::fusion::Domain::equalsTo(v_0));
    for (size_t i = 1; i < N - 1; ++i) {
      std::string cx_name = "cx_" + std::to_string(i);
      M->constraint(cx_name, v->index(i),
                    mosek::fusion::Domain::lessThan(v_bar.at(i)));
    }
    M->constraint("cx_N", v->index(N - 1),
                  mosek::fusion::Domain::equalsTo(v_s));

    for (size_t i = 0; i < N - 1; ++i) {
      std::string c1_name = "lc1" + std::to_string(i);
      M->constraint(
          c1_name,
          mosek::fusion::Expr::sub(
              mosek::fusion::Expr::mul(2.0 * dis.at(i) * a_max_l, y->index(i)),
              mosek::fusion::Expr::sub(v->index(i + 1), v->index(i))),
          mosek::fusion::Domain::greaterThan(0.0));

      std::string c2_name = "lc2" + std::to_string(i);
      M->constraint(
          c2_name,
          mosek::fusion::Expr::sub(
              mosek::fusion::Expr::mul(2.0 * dis.at(i) * a_max_l, y->index(i)),
              mosek::fusion::Expr::sub(v->index(i), v->index(i + 1))),
          mosek::fusion::Domain::greaterThan(0.0));

      std::string c3_name = "lc3" + std::to_string(i);
      M->constraint(c3_name,
                    mosek::fusion::Expr::add(v->index(i), v->index(i + 1)),
                    mosek::fusion::Domain::lessThan(h_theta.at(i)));

      std::string c4_name = "lc4" + std::to_string(i);
      M->constraint(c4_name,
                    mosek::fusion::Expr::add(v->index(i), v->index(i + 1)),
                    mosek::fusion::Domain::greaterThan(-h_theta.at(i)));

      {
        std::string c5_name = "lc5" + std::to_string(i);
        M->constraint(c5_name,
                      mosek::fusion::Expr::add(
                          mosek::fusion::Expr::mul(1 + g.at(i), v->index(i)),
                          mosek::fusion::Expr::mul(g.at(i), v->index(i + 1))),
                      mosek::fusion::Domain::lessThan(vr_max));

        std::string c7_name = "lc7" + std::to_string(i);
        M->constraint(c7_name,
                      mosek::fusion::Expr::sub(
                          mosek::fusion::Expr::mul(1 - g.at(i), v->index(i)),
                          mosek::fusion::Expr::mul(g.at(i), v->index(i + 1))),
                      mosek::fusion::Domain::lessThan(vr_max));
      }
    }
    for (size_t i = 0; i < N - 1; ++i) {
      std::string c3_name = "rc" + std::to_string(i);
      mosek::fusion::Constraint::t qc2 = M->constraint(
          c3_name,
          mosek::fusion::Expr::vstack(
              y->index(i),
              mosek::fusion::Expr::add(v->index(i + 1), v->index(i)),
              std::sqrt(2)),
          mosek::fusion::Domain::inRotatedQCone());
    }
    std::shared_ptr<ndarray<double, 1>> dis_expr = new_array_ptr(dis);
    M->objective("obj", mosek::fusion::ObjectiveSense::Minimize,
                 Expr::dot(dis_expr, y));
    auto time_before_solve = std::chrono::system_clock::now();

    std::chrono::duration<double> time_construct =
        time_before_solve - time_before_construct;
    std::cout << "TIIME: trajectory plan construct: " << time_construct.count()
              << std::endl;

    // Solve the problem
    M->solve();

    mosek::fusion::AccSolutionStatus status = M->getAcceptedSolutionStatus();
    if (status == mosek::fusion::AccSolutionStatus::Optimal) {
      std::cout << "The problem was solved to optimality." << std::endl;
    } else if (status == mosek::fusion::AccSolutionStatus::Feasible) {
      std::cout << "A primal feasible solution was found." << std::endl;
    } else {
      std::cout << "The solver did not find a satisfactory solution. "
                   "Status: "
                << status << std::endl;
    }
    auto time_after_solve = std::chrono::system_clock::now();
    std::chrono::duration<double> time_solve =
        time_after_solve - time_before_solve;
    std::cout << "TIIME: trajectory plan mosek solve time: "
              << time_solve.count() << std::endl;
    // Get the linear solution values
    auto xlvl = *(v->level());
    auto ylvl = *(y->level());
    // update time
    trajectory.front().time = 0;
    trajectory.front().v = xlvl[0];
    for (size_t i = 1; i < trajectory.size(); ++i) {
      double theta_i = theta.at(i - 1);
      double dis_i = dis.at(i - 1);
      double time = 2.0 * dis_i / (xlvl[i] + xlvl[i - 1]);
      trajectory.at(i).time = trajectory.at(i - 1).time + time;
      trajectory.at(i).v = xlvl[i];
      trajectory.at(i).w = theta_i / time;
    }
    std::cout << "trajectory duration time: " << trajectory.back().time
              << std::endl;
  }
  std::cout << "mosek CPP solver finished!" << std::endl;
}

void TOPPDWR::solveWithEcos(Trajectory &trajectory,
                            const TOPP_DWR_Params &params) {
  auto time_before_init = std::chrono::system_clock::now();
  std::cout << "path length: " << trajectory.back().s << " size "
            << trajectory.size() << std::endl;
  double v_max = params.v_max;
  double a_max_l = params.a_max_l;
  double a_max_n = params.a_max_n;
  double vr_max = params.vr_max;
  double vl_max = params.vl_max;
  double d = params.d;
  double w_max = params.w_max;
  double v_0 = std::min(params.start_vel, v_max);
  double v_s = std::min(params.end_vel, v_max);

  std::vector<double> v_bar;
  std::vector<double> dis;
  std::vector<double> theta;
  std::vector<double> g;
  std::vector<double> h_theta;
  size_t N = trajectory.size();
  v_bar.reserve(N);
  dis.reserve(N - 1);
  theta.reserve(N - 1);
  g.reserve(N - 1);
  h_theta.reserve(N - 1);

  double kappa = std::fabs(trajectory.front().curve);
  double v_bar_0 = std::min(v_max, std::sqrt(a_max_n / kappa));
  v_bar.emplace_back(v_bar_0);
  for (size_t i = 1; i < N; ++i) {
    auto first = trajectory.at(i - 1);
    auto second = trajectory.at(i);
    double dis_i_j = first.distanceTo(second);
    double theta_i_j = first.angle_diff(second);
    double kappa = std::fabs(trajectory.at(i).curve);
    double v_bar_i = std::min(v_max, std::sqrt(a_max_n / kappa));
    v_bar.emplace_back(v_bar_i);
    dis.emplace_back(dis_i_j);
    theta.emplace_back(theta_i_j);
    g.emplace_back(d * theta_i_j / (4.0 * dis_i_j));
    h_theta.emplace_back(
        std::min(2.0 * dis_i_j / std::fabs(theta_i_j) * w_max, 1e6));
  }
  if (v_0 > v_bar.front() + 1e-3)
    std::cout << " start vel error!, v_0: " << v_0 << " v_bar " << v_bar.front()
              << std::endl;
  if (v_s > v_bar.back() + 1e-3)
    std::cout << " end vel error!, v_s: " << v_s << " v_bar " << v_bar.back()
              << std::endl;
  auto time_after_init = std::chrono::system_clock::now();
  std::chrono::duration<double> time_init = time_after_init - time_before_init;
  std::cout << "TIIME: trajectory plan initializtion: " << time_init.count()
            << std::endl;
  int numv = N;
  int numc = N - 1;
  int num_vars = numv + numc;
  int numacc = N - 1;
  int numdec = N - 1;
  int numw = N - 1;
  int numvr = N - 1;
  int numvl = N - 1;
  int num_cones = N - 1;
  int num_equa = 2;

  int num_constri_per_v = 2;
  int num_constri_per_c = 1;
  int num_constri_per_w = 1;
  int num_constri_per_vr = 1;
  int num_constri_per_vl = 1;
  int num_constri_per_acc = 1;
  int num_constri_per_dec = 1;

  int num_constri_per_socs = 3;
  int num_constri_per_socs_v = 1;
  int num_constri_per_socs_c = 1;

  int num_constri_final_v = 2;

  int num_var_per_v = 1;
  int num_var_per_c = 1;
  int num_var_per_w = 2;
  int num_var_per_vr = 2;
  int num_var_per_vl = 2;
  int num_var_per_acc =
      2; //-v(i) + v(i+1) - 2*dis.at(i)*c(i)<=0,c(i) is considered later
  int num_var_c_per_acc =
      1; //-v(i) + v(i+1) - 2*dis.at(i)*c(i)<=0,c(i)is considered later
  int num_var_per_dec =
      2; // v(i) - v(i+1) - 2*dis.at(i)*c(i)<=0,c(i)is considered later
  int num_var_c_per_dec =
      1; // v(i) - v(i+1) - 2*dis.at(i)*c(i)<=0,c(i)is considered later
  int num_var_per_socs =
      2; // (c(i),v(i)+v(i+1),sqrt(2)),c(i)is considered later
  int num_var_c_per_socs =
      1; // (c(i),v(i)+v(i+1),sqrt(2)),c(i)is considered later

  int num_ineq_linear =
      num_constri_per_c * numc + num_constri_per_v * numv +
      num_constri_per_w * numw + num_constri_per_vr * numvr +
      num_constri_per_vl * numvl + num_constri_per_acc * numacc +
      num_constri_per_dec *
          numdec; // Total number of linear inequality constraints
  int num_ineq =
      num_ineq_linear + num_constri_per_socs *
                            num_cones; // Total number of inequality constraints
                                       // (linear + conic constraints)
  int num_ineq_per = num_constri_per_c + num_constri_per_v + num_constri_per_w +
                     num_constri_per_vr + num_constri_per_vl +
                     num_constri_per_acc +
                     num_constri_per_dec; // Number of linear inequality
                                          // constraints per variable 10
  int num_ineq_G_linear =
      num_var_per_c * num_constri_per_c * numc +
      num_var_per_v * num_constri_per_v * numv +
      num_var_per_w * num_constri_per_w * numw +
      num_var_per_vr * num_constri_per_vr * numvr +
      num_var_per_vl * num_constri_per_vl * numvl +
      (num_var_per_acc + num_var_c_per_acc) * num_constri_per_acc * numacc +
      (num_var_per_dec + num_var_c_per_dec) * num_constri_per_dec *
          numdec; // Number of linear constraints in matrix G

  int num_ineq_G = num_ineq_G_linear +
                   2 * num_constri_per_socs *
                       num_cones; // Total number of constraints in matrix G
  int num_ineq_G_c_per =
      num_var_per_c + num_constri_per_acc + num_constri_per_dec +
      2 * num_var_c_per_socs; // Number of elements in column i+N of matrix G
  int num_ineq_G_linear_per =
      0 * num_constri_per_c * num_var_per_c +
      num_var_per_v * num_constri_per_v + num_var_per_w * num_constri_per_w +
      num_var_per_vr * num_constri_per_vr +
      num_var_per_vl * num_constri_per_vl +
      num_var_per_acc * num_constri_per_acc +
      num_var_per_dec *
          num_constri_per_dec; // 16  Number of linear elements for each vi in G
  int num_ineq_G_v = num_ineq_G_linear_per * (numv - 1) + num_constri_final_v +
                     2 * num_var_per_socs * num_cones;
  int num_ineq_G_per_var =
      num_ineq_G_linear_per +
      4 * num_constri_per_socs_v; // Number of linear constraint elements in
                                  // column i (0<i<N) of matrix G
  int num_ineq_G_per_var_first =
      num_ineq_G_linear_per +
      2 * num_constri_per_socs_v; // Number of linear constraint elements in
                                  // column 0 of matrix G
  int *q = (int *)malloc(num_cones * sizeof(int));
  for (int i = 0; i < num_cones; i++) {
    q[i] =
        num_constri_per_socs; // Second-order cone dimension is 3, corresponding
                              // to (c_k, v_k + v_{k + 1}, \sqrt{2})
  }

  // Objective function coefficient vector c
  // std::shared_ptr<double[]> c = std::shared_ptr<double[]>(new double[N]);
  double *c = (double *)malloc(num_vars * sizeof(double));
  for (size_t i = 0; i < N; i++) {
    c[i] = 0.0; // Coefficient for v_k is initially set to 0
  }
  for (int i = N; i < num_vars; i++) {
    c[i] = 1.0 * 2 * dis[i - N]; // Coefficient for c_k
  }

  // Equality constraints: v_0, v_s
  pfloat *Apr = (double *)malloc(num_equa * sizeof(double));
  idxint *Ajc = (int *)malloc((num_vars + 1) * sizeof(int));
  idxint *Air = (int *)malloc(num_equa * sizeof(int));
  double *b = (double *)malloc(num_equa * sizeof(double));
  for (int i = 0; i < num_vars + 1; ++i) {
    if (i == 0) {
      Apr[0] = 1.0;
      Air[0] = 0;
      Ajc[0] = 0;
      b[0] = v_0;
    } else if (i < numv - 1) {
      Ajc[i] = Ajc[i - 1] + (i == 1 ? 1 : 0);
    } else if (i == numv - 1) {
      Apr[1] = 1.0;
      Air[1] = 1;
      b[1] = v_s;
      Ajc[i] = Ajc[i - 1];
    } else {
      Ajc[i] = Ajc[i - 1] + (i == numv ? 1 : 0);
    }
  }

  // Construction of matrices and vectors for inequality constraints
  pfloat *Gpr = (double *)malloc(num_ineq_G * sizeof(double));
  idxint *Gjc = (int *)malloc((num_vars + 1) * sizeof(int));
  idxint *Gir = (int *)malloc(num_ineq_G * sizeof(int));
  pfloat *h = (double *)malloc(num_ineq * sizeof(double));

  int index_cones_left = 0;
  int index_h_cone = 0;
  double cof = -std::sqrt(2.0) / 2.0;
  for (int i = 0; i < num_ineq; i++) {
    const int index_bias = (i > 0 && num_constri_per_socs > 0) ? -2 : 0;
    if (i <= numv - 2) {
      int index_Gpr = 0;
      int index_Gir = 0;
      int index_h = 0;
      int index_i = num_constri_per_c;
      int index_j = num_constri_per_v + num_constri_per_c;
      int index_c = num_constri_per_v + num_constri_per_c;
      int index_Gpr_c = 0;
      int index_Gir_c = 0;
      int num_ineq_G_per_col = num_ineq_G_per_var;
      // -v(i) <=0
      Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -1.0;
      Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
          num_ineq_per * i + index_i++;
      // v(i)<=v_bar(i)
      Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = 1.0;
      Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
          num_ineq_per * i + index_i++;
      //--------------Column i--------------------//
      if (num_constri_per_acc > 0) {
        // -v(i) + v(i+1) - 2*dis.at(i)*c(i)<=0: v(i)
        Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -1.0;
        Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
            num_ineq_per * i + index_i++;
      } /////
      if (num_constri_per_dec > 0) {
        // v(i) - v(i+1) - 2*dis.at(i)*c(i)<=0:v(i)
        Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = 1.0;
        Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
            num_ineq_per * i + index_i++;
      }
      // v(i) + v(i+1) <= w_bar(i): v(i)
      Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = 1.0;
      Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
          num_ineq_per * i + index_i++;
      // -(1+g(i))*v(i) - g(i)*v(i+1) <= vr_max: v(i)
      if (num_constri_per_vr > 1) {
        Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -(1 + g.at(i));
        Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
            num_ineq_per * i + index_i++;
      }
      // (1+g(i))*v(i) + g(i)*v(i+1) <= vr_max: v(i)
      Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = (1 + g.at(i));
      Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
          num_ineq_per * i + index_i++;
      // -(1-g(i))*v(i) + g(i)*v(i+1) <= vr_max: v(i)
      if (num_constri_per_vr > 1) {
        Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -(1 - g.at(i));
        Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
            num_ineq_per * i + index_i++;
      }
      // (1-g(i))*v(i) - g(i)*v(i+1) <= vr_max: v(i)
      Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = (1 - g.at(i));
      Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
          num_ineq_per * i + index_i++;

      if (num_constri_per_socs > 0) {
        if (i == 0) {
          // (sqrt(2), -(v(i)+v(i+1)),-c(i)): v(i)
          Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = 1.0 * cof;
          Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
              num_ineq_linear + index_cones_left;
          index_cones_left++;

          Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -1.0 * cof;
          Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
              num_ineq_linear + index_cones_left;
        }
        if (i > 0) {
          index_cones_left--;
          // (sqrt(2), -(v(i)+v(i+1)),-c(i)): v(i+1)
          Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = 1.0 * cof;
          Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
              num_ineq_linear + index_cones_left;
          index_cones_left++;
          Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -1.0 * cof;
          Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
              num_ineq_linear + index_cones_left;
          index_cones_left++;
          index_cones_left++;
          Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = 1.0 * cof;
          Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
              num_ineq_linear + index_cones_left;
          index_cones_left++;
          Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -1.0 * cof;
          Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
              num_ineq_linear + index_cones_left;
        }
      }
      //--------------Column i+1--------------------//
      if (num_constri_per_acc > 0) {
        // -v(i) + v(i+1) - 2*dis.at(i)*c(i)<=0: v(i+1)
        Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = 1.0;
        Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
            num_ineq_per * i + index_j++;
      }
      if (num_constri_per_dec > 0) {
        // v(i) - v(i+1) - 2*dis.at(i)*c(i)<=0: v(i+1)
        Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -1.0;
        Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
            num_ineq_per * i + index_j++;
      }
      // v(i) + v(i+1) <= w_bar(i): v(i+1)
      Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = 1.0;
      Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
          num_ineq_per * i + index_j++;
      // -(1+g(i))*v(i) - g(i)*v(i+1) <= vr_max: v(i+1)
      if (num_constri_per_vr > 1) {
        Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -g.at(i);
        Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
            num_ineq_per * i + index_j++;
      }
      // (1+g(i))*v(i) + g(i)*v(i+1) <= vr_max: v(i+1)
      Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = g.at(i);
      Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
          num_ineq_per * i + index_j++;
      // -(1-g(i))*v(i) + g(i)*v(i+1) <= vr_max: v(i+1)
      if (num_constri_per_vr > 1) {
        Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = g.at(i);
        Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
            num_ineq_per * i + index_j++;
      }
      // (1-g(i))*v(i) - g(i)*v(i+1) <= vr_max: v(i+1)
      Gpr[num_ineq_G_per_col * i + index_bias + index_Gpr++] = -g.at(i);
      Gir[num_ineq_G_per_col * i + index_bias + index_Gir++] =
          num_ineq_per * i + index_j++;

      //--------------Column i+N--------------------//
      if (num_constri_per_c > 0) {
        // -c(i)<=0
        Gpr[num_ineq_G_v + i * num_ineq_G_c_per + index_Gpr_c++] = -1.0;
        Gir[num_ineq_G_v + i * num_ineq_G_c_per + index_Gir_c++] =
            num_ineq_per * i + 0;
      }
      if (num_constri_per_acc > 0) {
        // -v(i) + v(i+1) - 2*dis.at(i)*c(i)<=0: c(i)
        Gpr[num_ineq_G_v + i * num_ineq_G_c_per + index_Gpr_c++] =
            -2 * dis.at(i) * a_max_l;
        Gir[num_ineq_G_v + i * num_ineq_G_c_per + index_Gir_c++] =
            num_ineq_per * i + index_c++;
      }
      if (num_constri_per_dec > 0) {
        // v(i) - v(i+1) - 2*dis.at(i)*c(i)<=0: c(i)
        Gpr[num_ineq_G_v + i * num_ineq_G_c_per + index_Gpr_c++] =
            -2 * dis.at(i) * a_max_l;
        Gir[num_ineq_G_v + i * num_ineq_G_c_per + index_Gir_c++] =
            num_ineq_per * i + index_c++;
      }
      if (num_constri_per_socs > 0) {
        // (sqrt(2), -(v(i)+v(i+1)),-c(i)): c(i)
        Gpr[num_ineq_G_v + i * num_ineq_G_c_per + index_Gpr_c++] = 1.0 * cof;
        Gir[num_ineq_G_v + i * num_ineq_G_c_per + index_Gir_c++] =
            num_ineq_linear + num_constri_per_socs * i;
        // (sqrt(2), -(v(i)+v(i+1)),-c(i)): c(i)
        Gpr[num_ineq_G_v + i * num_ineq_G_c_per + index_Gpr_c++] = 1.0 * cof;
        Gir[num_ineq_G_v + i * num_ineq_G_c_per + index_Gir_c++] =
            num_ineq_linear + num_constri_per_socs * i + 1;
      }
      if (num_constri_per_c > 0) {
        // -c(i) <=0
        h[num_ineq_per * i + index_h++] =
            0 * -1.0 / (v_bar.at(i) + v_bar.at(i + 1));
      }
      // -v(i) <=0
      h[num_ineq_per * i + index_h++] = 0.0;
      // v(i)<=v_bar(i)
      h[num_ineq_per * i + index_h++] = v_bar.at(i);
      // -c(i) <=0
      // h[num_ineq_per * i + index_h++] = 0.0;
      if (num_constri_per_acc > 0) {
        // -v(i) + v(i+1) - 2*dis.at(i)*c(i)<=0
        h[num_ineq_per * i + index_h++] = 0;
      }
      if (num_constri_per_dec > 0) {
        // v(i) - v(i+1) - 2*dis.at(i)*c(i)<=0
        h[num_ineq_per * i + index_h++] = 0;
      }
      // v(i) + v(i+1) <= w_bar(i)
      h[num_ineq_per * i + index_h++] = h_theta.at(i);
      // -(1+g(i))*v(i) - g(i)*v(i+1) <= -(vr_min)
      if (num_constri_per_vr > 1) {
        h[num_ineq_per * i + index_h++] = vr_max;
      }
      // (1+g(i))*v(i) + g(i)*v(i+1) <= vr_max
      h[num_ineq_per * i + index_h++] = vr_max;
      // -(1-g(i))*v(i) + g(i)*v(i+1) <= -(-vl_min)
      if (num_constri_per_vr > 1) {
        h[num_ineq_per * i + index_h++] = vr_max;
      }
      // (1-g(i))*v(i) - g(i)*v(i+1) <= vr_max
      h[num_ineq_per * i + index_h++] = vr_max;
    } else if (i == numv - 1) {
      int index_h = 0;
      // -v(n-1) <=0
      Gpr[num_ineq_G_per_var * i + index_bias] = -1.0;
      Gir[num_ineq_G_per_var * i + index_bias] = num_ineq_per * i;
      h[num_ineq_per * i + index_h++] = 0.0;
      // v(n-1)<=v_bar(n-1)
      Gpr[num_ineq_G_per_var * i + 1] = 1.0;
      Gir[num_ineq_G_per_var * i + index_bias + 1] = num_ineq_per * i + 1;
      h[num_ineq_per * i + index_h++] = v_bar.at(i);
      if (num_constri_per_socs > 0) {
        Gpr[num_ineq_G_per_var * i + index_bias + 2] = 1.0 * cof;
        Gir[num_ineq_G_per_var * i + index_bias + 2] =
            num_ineq - num_constri_per_socs;
        Gpr[num_ineq_G_per_var * i + index_bias + 3] = -1.0 * cof;
        Gir[num_ineq_G_per_var * i + index_bias + 3] =
            num_ineq - num_constri_per_socs + 1;
      }
    } else if (i <= num_vars - 1) {
      if (num_constri_per_socs) {
        h[num_ineq_linear + index_h_cone++] = 0.0;
        h[num_ineq_linear + index_h_cone++] = 0.0;
        h[num_ineq_linear + index_h_cone++] = -std::sqrt(2.0);
      }
    }
  }
  for (int i = 0; i < num_vars + 1; i++) {
    if (i == 0) {
      Gjc[i] = i;
    } else if (i <= numv - 1) {
      Gjc[i] =
          Gjc[i - 1] +
          (i == 1
               ? (num_ineq_per - num_constri_per_c + 2 * num_constri_per_socs_v)
               : (num_ineq_G_linear_per + 4 * num_constri_per_socs_v));
    } else {
      Gjc[i] =
          Gjc[i - 1] +
          (i == numv
               ? (num_ineq_per - num_constri_per_c + 2 * num_constri_per_socs_v)
               : (num_constri_per_c + num_constri_per_acc +
                  num_constri_per_dec + 2 * num_constri_per_socs_c));
    }
  }

  // initialize ECOS solver
  pwork *prob;
  std::cout << " EOCS NUMBER num_ineq_G_per_var: " << num_ineq_G_per_var
            << " num_ineq_per: " << num_ineq_per << " num_ineq: " << num_ineq
            << " num_ineq_G: " << num_ineq_G << " num_ineq_G_per_var_first "
            << num_ineq_G_per_var_first << " num_ineq_linear "
            << num_ineq_linear << std::endl;
  if (num_constri_per_socs > 0)
    prob = ECOS_setup(num_vars, num_ineq, num_equa, num_ineq_linear, num_cones,
                      q, 0, Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
  else
    prob = ECOS_setup(num_vars, num_ineq, num_equa, num_ineq_linear, 0, NULL, 0,
                      Gpr, Gjc, Gir, Apr, Ajc, Air, c, h, b);
  if (prob == NULL) {
    std::cout << " ECOS setup failed! " << std::endl;
    free(q);
    free(c);
    free(Ajc);
    free(Air);
    free(Apr);
    free(b);
    free(Gpr);
    free(Gjc);
    free(Gir);
    free(h);
  } else {
    // solve with ecos
    auto time_after_build = std::chrono::system_clock::now();
    std::chrono::duration<double> time_build =
        time_after_build - time_before_init;
    std::cout << "TIIME: build problem [ecos]: " << time_build.count()
              << " path time: " << trajectory.back().time << std::endl;
    idxint exitflag = ECOS_solve(prob);
    if (exitflag == ECOS_OPTIMAL) {
      std::cout << " Ecos: Optimal solution found! " << std::endl;
    } else {
      std::cout << " ECOS failed with exit flag: " << exitflag << std::endl;
    }
    bool apply = true;
    if (apply) {
      trajectory.front().time = 0;
      trajectory.front().v = prob->x[0];
      for (size_t i = 1; i < trajectory.size(); ++i) {
        double theta_i = theta.at(i - 1);
        double dis_i = dis.at(i - 1);
        double time = 2.0 * dis_i / (prob->x[i] + prob->x[i - 1]);
        trajectory.at(i).time = trajectory.at(i - 1).time + time;
        trajectory.at(i).v = prob->x[i];
        trajectory.at(i).w = theta_i / time;
      }
      auto time_after_solve = std::chrono::system_clock::now();
      std::chrono::duration<double> time_solve =
          time_after_solve - time_before_init;
      std::cout << "TIIME: trajectory plan [ecos]: " << time_solve.count()
                << ", trajectory duration time: " << trajectory.back().time
                << std::endl;
    }
    // clear
    ECOS_cleanup(prob, 0);
    free(q);
    free(c);
    free(Ajc);
    free(Air);
    free(Apr);
    free(b);
    free(Gpr);
    free(Gjc);
    free(Gir);
    free(h);
    std::cout << " ecos exitflag: " << exitflag << std::endl;
  }
}

} // namespace TOPP_DWR