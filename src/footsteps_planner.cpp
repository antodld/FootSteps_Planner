#include "footsteps_planner.h"

namespace mc_plugin
{
namespace footsteps_planner
{

void FootStepGen::init(std::string supportFootName,
                       Footstep P_f0,
                       const std::vector<sva::MotionVecd> V,
                       const std::vector<double> Tstep,
                       std::vector<Footstep> ref_pose)
{

  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();

  plan_.support_foot_name(supportFootName);
  plan_.support_foot(P_f0);
  N_steps = -1;
  supportFoot = supportFootName;
  v_inputs_ = V;
  P_traj_.clear();
  P_ = static_cast<int>(Tp_ / delta_);
  for(size_t k = 0; k < v_inputs_.size(); k++)
  {
    P_traj_.push_back(IntegrateVelProfile(k));
  }
  t_steps_inputs_ = Tstep;

  if(t_steps_inputs_.size() != 0)
  {
    t_steps_inputs_[0] = std::max(Ts_min_, t_steps_inputs_[0]);
  }
  // while(t_steps_inputs_.back() < Tp_)
  // {
  //   t_steps_inputs_.push_back(t_steps_inputs_.back() + Ts_);
  // }
  // while(t_steps_inputs_.back() > Tp_)
  // {
  //   t_steps_inputs_.pop_back();
  // }


  for(size_t k = 0; k < ref_pose.size(); k++)
  {
    if(std::abs(ref_pose[k].ori()) > M_PI)
    {
      ref_pose[k].ori(ref_pose[k].ori() - (ref_pose[k].ori() / std::abs(ref_pose[k].ori())) * 2 * M_PI);
    }
  }
  pose_reference_ = ref_pose;

  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  // double ProcessTime = time_span.count();
  // std::cout << "Footstep plan init done" << std::endl;
  // mc_rtc::log::success("FootSteps Initialized in : " + std::to_string(ProcessTime));
}

void FootStepGen::GetStepsTimings()
{
  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();

  StepsTimings_.clear();
  StepsTimings_ = t_steps_inputs_;
  StepsTimings_indx_.clear();
  FootSteps_indx_.clear();
  // 1 Generate reference trajectory and velocity according to order V -> Pf -> Ts

  if(pose_reference_.size() != 0)
  {
    
      // mc_rtc::log::info(
      //     "[Steps Timings Generation] generating ...");
      ref_traj_point P_f_im1 = IntegrateVelProfile(0);
      ref_traj_point P_f_i(pose_reference_[0].pose(), pose_reference_[0].ori());
      const double d_step = Eigen::Vector2d{d_h_x,d_h_y}.norm()/2;
      const size_t nsteps = static_cast<size_t>( (P_f_im1.pose() - P_f_i.pose()).norm() / d_step);
      if( nsteps < StepsTimings_.size())
      {
        P_traj_ = GetRefTrajectory(P_f_im1, P_f_i,StepsTimings_[nsteps]);
        while (P_traj_.size() < static_cast<size_t>(P_))
        {
          P_traj_.push_back(P_traj_.back());
        }
        
      }
      else
      {
        P_traj_ = GetRefTrajectory(P_f_im1, P_f_i,Tp_);
      }
      // mc_rtc::log::info("Traj L gen {}",P_traj_.size());

  }        
  

  for(size_t i = 0; i < StepsTimings_.size(); i++)
  {
    double t_i = StepsTimings_[i];
    // StepsTimings_indx_[i] = (int)std::round(t_i / delta_);
    StepsTimings_indx_.push_back((int)std::round(t_i / delta_) );
    FootSteps_indx_.push_back(-1);
  }
  while(StepsTimings_[StepsTimings_.size() - 1] > Tp_ - delta_)
  {
    StepsTimings_.pop_back();
    FootSteps_indx_.pop_back();
    StepsTimings_indx_.pop_back();
  }

  // mc_rtc::log::success("[Steps Timings Generation] 3 OK");
  F_ = static_cast<int>(StepsTimings_.size());
  // mc_rtc::log::info("Params //");
  // for(size_t k = 0; k < static_cast<size_t>(F_); k++)
  // {
  //   mc_rtc::log::info("Ts : {} ", StepsTimings_[k]);
  //   mc_rtc::log::info("L traj: {} ", P_traj_.size());
  //   mc_rtc::log::info(StepsTimings_indx_[k]);
  //   // mc_rtc::log::info("Foot: {} ", FootSteps_indx_[k]);
  //   mc_rtc::log::info("//");
  // }

  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  // double ProcessTime = time_span.count();
  // mc_rtc::log::success("StepTiming gened in: " + std::to_string(ProcessTime) + " ms");
}

Footsteps_plan FootStepGen::compute_plan()
{
  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();
  plan_.clear();
  GetStepsTimings();

  Theta_f_.resize(F_, 1);

  if(F_ <= 0)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("[footstepsplanner::compute_plan] No step to compute");
    return plan_;
  }

  // Solving QP 1 for Orientation

  Eigen::VectorXd Dtheta_upper_lim(Eigen::VectorXd::Ones(F_) * max_theta);
  Eigen::VectorXd Dtheta_lower_lim(Eigen::VectorXd::Ones(F_) * -max_theta);
  Eigen::MatrixXd Delta =
      Eigen::MatrixXd::Identity(F_,   F_); // Differenciation matrix between two steps orientation (for constraints)
  Aeq = Eigen::MatrixXd::Zero(F_, F_);
  beq = Eigen::VectorXd::Zero(Aeq.rows());

  for(int k = 0; k < F_; ++k)
  {
    if(k == 0)
    {
      Dtheta_upper_lim(k) += plan_.support_foot().ori();
      Dtheta_lower_lim(k) += plan_.support_foot().ori();
    }

    else
    {

      Delta(k, k - 1) = -1;
    }

    if(FootSteps_indx_[k] != -1)
    {
      double theta_cstr = pose_reference_[FootSteps_indx_[k]].ori();
      if(P_traj_[StepsTimings_indx_[k]].ori() - pose_reference_[FootSteps_indx_[k]].ori() > M_PI)
      {
        theta_cstr -=
            pose_reference_[FootSteps_indx_[k]].ori() / std::abs(pose_reference_[FootSteps_indx_[k]].ori()) * 2 * M_PI;
      }
      Aeq(k, k) = 1;
      beq(k) = theta_cstr;
    }
  }


  Aineq = Eigen::MatrixXd::Zero(2 * Delta.rows(), F_);
  bineq = Eigen::VectorXd::Zero(Aineq.rows());
  Aineq << Delta, -Delta;
  bineq << Dtheta_upper_lim, -Dtheta_lower_lim;

  Eigen::VectorXd b = Eigen::VectorXd::Zero(F_);
  for(int i = 0 ; i < F_ ; i++)
  {

    b(i) = P_traj_[StepsTimings_indx_[i]].ori();

  }
  Q_ = Eigen::MatrixXd::Identity(F_, F_);
  p_ = -b;
  Eigen::VectorXd theta = solveQP();
  if(!QPsuccess)
  {
    mc_rtc::log::error("[Footsteps planner] Step theta failed");
  }

  for(int k = 0; k < theta.size(); k++)
  {

    if(std::abs(theta(k)) > M_PI)
    {
      Theta_f_(k) = theta(k) - (theta(k) / std::abs(theta(k))) * 2 * M_PI;
    }
    else
    {
      Theta_f_(k) = theta(k);
    }
  }


  // Solving QP 2 For placement

  Delta =
      Eigen::MatrixXd::Identity(2 * F_, 2 * F_); // Differenciation matrix between two steps location (for constraints)
  Aeq = Eigen::MatrixXd::Zero(2 * F_, 2 * F_);
  beq = Eigen::VectorXd::Zero(2 * F_);

  double sgn = -1.0;
  if(supportFoot == "RightFoot")
  {
    sgn = 1.0;
  }
  double l = l_;
  if(d_h_y/2 > l_ - 0.5 * d_min)
  {
    l = 0.5 * (d_h_y + d_min);
  }

  std::vector<Eigen::VectorXd> cstr_vec;
  std::vector<Eigen::MatrixX2d> Normal_Vec;

  Eigen::Matrix2d R_Theta_0;
  double theta_0 = P_traj_[0].ori();
  R_Theta_0 = sva::RotZ(-theta_0).block(0, 0, 2, 2);

  Admissible_Region Kinematic_Rectangle(Eigen::Vector3d{0, 0, plan_.support_foot().ori()},
                                        Eigen::Vector3d{d_h_x, d_h_y, 0});

  Eigen::VectorXd bcstr = Kinematic_Rectangle.Offset
                          + Kinematic_Rectangle.Polygone_Normals.transpose()
                                * (plan_.support_foot().pose() + sgn * R_Theta_0 * Eigen::Vector2d{0, l});

  Normal_Vec.push_back(Kinematic_Rectangle.Polygone_Normals.transpose());
  cstr_vec.push_back(bcstr);

  for(int k = 1; k < F_; k++)
  {

    double theta_k = P_traj_[StepsTimings_indx_[k]].ori();
    double theta_f_km1 = Theta_f_(k - 1);

    Eigen::Matrix2d R_Theta_k;
    R_Theta_k = sva::RotZ(-theta_k).block(0, 0, 2, 2);


    Kinematic_Rectangle = Admissible_Region(Eigen::Vector3d{0, 0, theta_f_km1}, Eigen::Vector3d{d_h_x, d_h_y, 0});
    Eigen::Matrix2d R_Theta_f_km1 = sva::RotZ(-theta_f_km1).block(0, 0, 2, 2);

    bcstr = Kinematic_Rectangle.Offset
            + sgn * (1 - 2 * (((k) % 2))) * Kinematic_Rectangle.Polygone_Normals.transpose() * R_Theta_f_km1
                  * Eigen::Vector2d{0, l};

    Normal_Vec.push_back(Kinematic_Rectangle.Polygone_Normals.transpose());
    cstr_vec.push_back(bcstr);
  }

  if(FootSteps_indx_[0] != -1)
  {
    Aeq(0, 0) = 1;
    Aeq(1, 1) = 1;
    beq(0) = steps_inputs_[FootSteps_indx_[0]].pose().x();
    beq(1) = steps_inputs_[FootSteps_indx_[0]].pose().y();
  }
  for(int k = 1; k < F_; k++)
  {

    Delta(2 * k, 2 * (k - 1)) = -1;
    Delta(2 * k + 1, 2 * (k - 1) + 1) = -1;

    if(FootSteps_indx_[k] != -1)
    {
      Aeq(2 * k, 2 * k) = 1;
      Aeq(2 * k + 1, 2 * k + 1) = 1;
      beq(2 * k) = steps_inputs_[FootSteps_indx_[k]].pose().x();
      beq(2 * k + 1) = steps_inputs_[FootSteps_indx_[k]].pose().y();
    }
  }

  Eigen::Index N_cstr = 0;
  for(size_t k = 0; k < static_cast<size_t>(Normal_Vec.size()); k++)
  {
    N_cstr += Normal_Vec[k].rows();
  }

  Aineq = Eigen::MatrixXd::Zero(N_cstr, 2 * F_);
  bineq = Eigen::VectorXd::Zero(N_cstr);

  Eigen::Index step_indx = 0;
  Eigen::Index cstr_index = 0;
  for(size_t i_ineq = 0; i_ineq < Normal_Vec.size(); i_ineq++)
  {

    Eigen::MatrixXd n_vec = Normal_Vec[i_ineq];
    Eigen::VectorXd ineq = cstr_vec[i_ineq];

    for(Eigen::Index cstr = 0; cstr < n_vec.rows(); cstr++)
    {

      Aineq(cstr_index + cstr, step_indx) = n_vec(cstr, 0);
      Aineq(cstr_index + cstr, step_indx + 1) = n_vec(cstr, 1);

      bineq(cstr_index + cstr) = ineq(cstr);
    }

    step_indx += 2;
    cstr_index += n_vec.rows();
  }

  Aineq = Aineq * Delta;

  b = Eigen::VectorXd::Zero(2 * F_);
  for(int i = 0; i < F_; i++)
  {
    double theta_i = P_traj_[i].ori();
    Eigen::Matrix2d R_Theta_i;
    R_Theta_i = sva::RotZ(-theta_i).block(0, 0, 2, 2);
    Eigen::Vector2d dl = (sgn * (1 - 2 * (((i) % 2))) * R_Theta_i * Eigen::Vector2d{0, l_/2});
    b(2 * i) = P_traj_[StepsTimings_indx_[i]].pose().x() + dl.x();
    b(2 * i + 1) = P_traj_[StepsTimings_indx_[i]].pose().y() + dl.y();
  }

  Q_ = Eigen::MatrixXd::Identity(2 * F_, 2 * F_);
  p_ = -b;
  // Aeq = Q_; beq = b;
  // mc_rtc::log::info("Step Aineq {}",Aineq);
  // mc_rtc::log::info("Step bineq {}",bineq);
  // mc_rtc::log::info("Step Aeq {}",Aeq);
  // mc_rtc::log::info("Step beq {}",beq);
  // mc_rtc::log::info("Step Q {}",Q_);
  // mc_rtc::log::info("Step p {}",p_);
  Eigen::VectorXd XY(solveQP());
  if(!QPsuccess)
  {
    mc_rtc::log::error("[Footsteps planner] Step QP failed");
  }

  // mc_rtc::log::info("Step out F {}\n{}", F_, XY);

  plan_.ori_offset(theta_offset_);

  for(Eigen::Index k = 0; k < F_; k++)
  {
    double xf = XY(2 * k);
    double yf = XY(2 * k + 1);
    plan_.add(Footstep(Eigen::Vector2d{xf, yf}, Theta_f_(k), StepsTimings_[static_cast<size_t>(k)],
                       Eigen::Vector2d{0.1, 0.1}));
  }

  const double theta_legs_0 = plan_.support_foot().ori();
  const double theta_legs_1 = Theta_f_(0);

  Eigen::Matrix2d A;
  A.block(0,0,2,1) << cos(theta_legs_0),sin(theta_legs_0);
  A.block(0,1,2,1) << -cos(theta_legs_1),-sin(theta_legs_1);
  if(A.determinant() != 0)
  {
    const Eigen::Matrix2d R_sup_0 =  plan_.support_foot().PT_pose().rotation().block(0,0,2,2).transpose();
    // const Eigen::Matrix2d R_sup_1 =  plan_.steps_PTpose()[0].rotation().block(0,0,2,2).transpose();
    const Eigen::Vector2d coeff = A.inverse() * (plan_.steps_pose()[0].segment(0,2) - plan_.support_foot().pose());
    intersec = plan_.support_foot().pose() + R_sup_0 * Eigen::Vector2d{1,0} * coeff(0);
    const double proj = (plan_.steps_pose()[0].segment(0,2) - plan_.support_foot().pose()).normalized().transpose() * (intersec - plan_.support_foot().pose());

    r = (intersec - plan_.support_foot().pose()) - (plan_.steps_pose()[0].segment(0,2) - plan_.support_foot().pose()).normalized() * proj;
    if((R_sup_0.transpose() * r).x() > 0 && r.norm() < 5)
    {
      mc_rtc::log::warning("[Footsteps planner] Potential Legs collision");
      // Eigen::Vector2d step_pose = plan_.steps_PTpose()[0].translation().segment(0,2) + R_sup_1.transpose() * Eigen::Vector2d{0,sgn * 0.1};
      // plan_.edit(Footstep(step_pose, (plan_.steps()[0]).ori_, StepsTimings_[0],Eigen::Vector2d{0.1, 0.1}),0);
      plan_.edit(Footstep(plan_.steps()[0].pose_, plan_.support_foot().ori(), StepsTimings_[0],Eigen::Vector2d{0.1, 0.1}),0);

    }
  }


  // mc_rtc::log::success("Position OK");
  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  // double ProcessTime = time_span.count();
  // mc_rtc::log::success("FootSteps computed in : {} ms",ProcessTime);
  return plan_;
}

std::vector<ref_traj_point> FootStepGen::GetRefTrajectory(ref_traj_point & P_s_0, ref_traj_point & P_s_1,const double duration)
{
// std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();
  const Eigen::Matrix2d R_s0_0 = P_s_0.PT_pose().rotation().block(0,0,2,2).transpose();
  // const Eigen::Matrix2d R_s1_0 = P_s_1.PT_pose().rotation().block(0,0,2,2).transpose();

  bool shuffle = false;
  bool backward = false;

  if( std::abs( (R_s0_0.transpose()* (P_s_1.pose() - P_s_0.pose())).x()) < 3e-1)
  {
    shuffle = true;
  }
  if( (R_s0_0.transpose()* (P_s_1.pose() - P_s_0.pose())).x() < 0 && !shuffle)
  {
    backward = true;
  }

  const Eigen::Vector2d init_pose = P_s_0.pose();
  const Eigen::Vector2d target_pose = P_s_1.pose();

  if(!shuffle)
  {
    Eigen::Vector2d init_ori = {cos(P_s_0.ori()), sin(P_s_0.ori())};
    Eigen::Vector2d target_ori = {cos(P_s_1.ori()), sin(P_s_1.ori())};
    if(backward)
    {
      init_ori = {-cos(P_s_0.ori()), -sin(P_s_0.ori())};
      target_ori = {-cos(P_s_1.ori()), -sin(P_s_1.ori())};
    }
    path.reset(init_pose, init_ori, target_pose, target_ori);
  }

  const int N = static_cast<int>(duration/delta_);

  std::vector<ref_traj_point> Output;
  for(int k = 0; k < N + 1; k++)
  {
    double t = ((double)k) / ((double)N);
    if(!shuffle)
    {
      Eigen::Vector2d pos_t = path.pos(t);
      Eigen::Vector2d ori_t = path.tangent(t);
      double theta = atan2(ori_t.y(), ori_t.x());
      // std::cout << "theta :" << theta  << std::endl;
      if((init_pose - target_pose).norm() < 5e-2)
      {
        theta = P_s_1.ori();
      }
      if(backward)
      {
        Output.push_back(ref_traj_point(pos_t, theta - M_PI));
      }
      else
      {
        Output.push_back(ref_traj_point(pos_t, theta));
      }
      // std::cout << "in traj :" <<Output.back().ori() << std::endl;
    }
    else
    {
      Output.push_back(ref_traj_point(P_s_0.pose() + t * (P_s_1.pose() - P_s_0.pose()) , P_s_1.ori()));
    }

    
  }

  // std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  // double ProcessTime = time_span.count();
  // mc_rtc::log::success("Ref trajectory computed in : {}", ProcessTime);

  return Output;
}

Eigen::VectorXd FootStepGen::solveQP()
{

  Eigen::QuadProgDense QP;

  int Nvar = static_cast<int>(Q_.rows());
  int NIneqConstr = static_cast<int>(Aineq.rows());
  int NEqConstr = static_cast<int>(Aeq.rows());
  // QP.tolerance(1e-4);
  QP.problem(Nvar, NEqConstr, NIneqConstr);
  QPsuccess = QP.solve(Q_, p_, Aeq, beq, Aineq, bineq);

  return QP.result();
}

ref_traj_point FootStepGen::IntegrateVelProfile(size_t k_end)
{

  size_t k_start = 0;
  if(k_end < P_traj_.size())
  {
    return P_traj_[k_end];
  }

  ref_traj_point Output(plan_.support_foot().pose(), plan_.support_foot().ori());
  if(P_traj_.size() == 0)
  {

    if(supportFoot == "RightFoot")
    {
      Output.pose_ += sva::RotZ(-Output.ori_).block(0, 0, 2, 2) * Eigen::Vector2d{0, l_ / 2};
    }
    else
    {
      Output.pose_ -= sva::RotZ(-Output.ori_).block(0, 0, 2, 2) * Eigen::Vector2d{0, l_ / 2};
    }
  }
  else
  {
    Output.pose_ = P_traj_.back().pose();
    Output.ori_ = P_traj_.back().ori();
    k_start = P_traj_.size() - 1;
  }

  for(size_t i = k_start; i < std::min(k_end, v_inputs_.size()); i++)
  {
    Output.ori_ += v_inputs_[i].angular().z() * delta_;
    Output.pose_ += sva::RotZ(-Output.ori_).block(0, 0, 2, 2) * v_inputs_[i].linear().segment(0, 2) * delta_;
  }

  return Output;
}

} // namespace footsteps_planner
} // namespace mc_plugin
