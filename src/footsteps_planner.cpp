#include "footsteps_planner.h"

namespace mc_plugin
{
namespace footsteps_planner
{

void FootStepGen::Init(std::string supportFootName,
                       Footstep P_f0,
                       const std::vector<sva::MotionVecd> & V,
                       const std::vector<double> & Tstep,
                       std::vector<Footstep> & Pf)
{

  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();

  plan_.clear();
  plan_.support_foot(P_f0);
  N_steps = -1;
  supportFoot = supportFootName;
  v_inputs_ = V;

  P_traj_.clear();
  P_traj_.push_back(IntegrateVelProfile(0));
  for(int i = 0; i < 3; i++)
  {
    if(std::abs(P_traj_.back().pose()(i)) < 1e-5)
    {
      P_traj_.back().pose_(i) = 0;
    }
  }
  for(int k = 1; k < (int)v_inputs_.size(); k++)
  {
    P_traj_.push_back(IntegrateVelProfile(k));
  }
  t_steps_inputs_ = Tstep;

  if(t_steps_inputs_.size() != 0)
  {
    t_steps_inputs_[0] = std::max(Ts_min_, t_steps_inputs_[0]);
  }

  for(int k = 0; k < Pf.size(); k++)
  {
    if(std::abs(Pf[k].ori()) > M_PI)
    {
      Pf[k].ori(  Pf[k].ori() - (Pf[k].ori() / std::abs(Pf[k].ori())) * 2 * M_PI );
    }

  }

  // mc_rtc::log::info("Sizes at start");
  // mc_rtc::log::info(P_traj_.size());
  // mc_rtc::log::info(Vx_.size());

  // mc_rtc::log::info(P_f_0.z());
  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  double ProcessTime = time_span.count();
  // mc_rtc::log::success("FootSteps Initialized in : " + std::to_string(ProcessTime));
}

void FootStepGen::GetStepsTimings()
{
  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();

  StepsTimings_.clear();
  StepsTimings_indx_.clear();
  FootSteps_indx_.clear();
  // 1 Generate reference trajectory and velocity according to order V -> Pf -> Ts

  if(v_inputs_.size() != P_)
  {
    // 1-1 Check if the reference trajectory and vel can be drawn from Pf_in;
   
    if(steps_inputs_.size() != 0)
    {
      mc_rtc::log::info("[Steps Timings Generation] Generation by Pf_in");
      int ki_start = 0;
      for(int kfoot = 0; kfoot < steps_inputs_.size(); kfoot++)
      {

        int ki = ki_start;
        while(ki <= v_inputs_.size())
        {
          if(v_inputs_.size() == P_)
          {
            break;
          }

          ref_traj_point P_ki = P_traj_[std::min(ki, (int)P_traj_.size() - 1)];

     
          double dist_i_x = std::abs((sva::RotZ(-P_ki.ori()).block(0,0,2,2).transpose() * ( P_ki.pose() - steps_inputs_[kfoot].pose())).x());
          double dist_i_y = std::abs((sva::RotZ(-P_ki.ori()).block(0,0,2,2).transpose() * ( P_ki.pose() - steps_inputs_[kfoot].pose())).x());
          double dist_i_theta = P_ki.ori() - steps_inputs_[kfoot].ori();

          if(std::abs(dist_i_theta) > M_PI)
          {
            dist_i_theta -= (dist_i_theta / std::abs(dist_i_theta)) * 2 * M_PI;
          }
          if(dist_i_x < d_h_x / 4 && dist_i_y < l_ / 2 + d_h_y
             && std::abs(dist_i_theta) < max_theta / 2) // Input step near reference trajectory
          {
            // mc_rtc::log::info("[Steps Timings Generation] Input step near reference trajectory");
            ki = ki_start;
            break;
          }
          if(ki == v_inputs_.size()) // Input step too far from reference trajectory, generating ...
          {
            // mc_rtc::log::info(
            //     "[Steps Timings Generation] Input step too far from reference trajectory, generating ...");
            ref_traj_point P_f_im1 = P_traj_[P_traj_.size() - 1];
            ref_traj_point P_f_i(steps_inputs_[kfoot].pose(),steps_inputs_[kfoot].ori());

            ref_traj_point P_f_i_pl(P_f_i.pose() + sva::RotZ(-P_f_i.ori()).block(0,0,2,2) * Eigen::Vector2d{0, l_ / 2} , P_f_i.ori());
            ref_traj_point P_f_i_ml(P_f_i.pose() + sva::RotZ(-P_f_i.ori()).block(0,0,2,2) * Eigen::Vector2d{0, - l_ / 2} , P_f_i.ori());

            if((P_f_i.pose() - P_f_im1.pose()).norm() < 1)
            {

              if((sva::RotZ(-P_f_im1.ori()).transpose().block(0,0,2,2) 
                  * (P_f_i.pose() - P_f_im1.pose())).y() < 0)
              {
                P_f_i = P_f_i_pl;
              }
              else
              {
                P_f_i = P_f_i_ml;
              }
            }

            std::vector<ref_traj_point> Traj = GetRefTrajectory(P_f_im1, P_f_i);
            std::vector<sva::MotionVecd> V;
            if(Traj.size() > 2)
            {
              for(int k = 0; k < Traj.size() - 1; k++)
              {
                sva::MotionVecd v_k((Traj[k+1].rpy_ori() - Traj[k].rpy_ori())/delta_  , (Traj[k+1].vec3_pose() - Traj[k].vec3_pose())/delta_  );
                V.push_back(v_k);
              }
            }
            else
            {
              V.push_back(sva::MotionVecd::Zero());
              V.push_back(sva::MotionVecd::Zero());
            }

            for(int k = 0; k < std::min((int)V.size(), (int) P_); k++)
            {

              v_inputs_.push_back(V[k]);

              ref_traj_point next_P_traj = P_traj_.back();
              next_P_traj.ori_ += v_inputs_.back().angular().z() * delta_;
              next_P_traj.pose_ += sva::RotZ(-next_P_traj.ori()).block(0,0,2,2) * v_inputs_[k].linear() * delta_;
              for(int i = 0; i < 3; i++)
              {
                if(std::abs(next_P_traj.pose()(i)) < 1e-5)
                {
                  next_P_traj.pose_(i) = 0;
                }
              }
              P_traj_.push_back(next_P_traj);
            }

            ki_start = v_inputs_.size() - 1;
          }

          ki += 1;
        }

        // if (kfoot + 1 != Xs_in_.size()) //delete Pf_s outside Tp
        // {
        //   for (int k = 0 ; k < Xs_in_.size() - kfoot - 1 ; k++)
        //   {
        //     Xs_in_.pop_back();
        //     Ys_in_.pop_back();
        //     Theta_s_in_.pop_back();
        //   }
        // }
      }
    }
   
    // 1-2 Use timesteps to generate forward velocity
   
    if(t_steps_inputs_.size() != 0 && v_inputs_.size() != P_)
    {

      double ti = v_inputs_.size() * delta_;
      int kTsteps = 0;

      while(t_steps_inputs_[kTsteps] - ti < Ts_min_ && kTsteps < t_steps_inputs_.size())
      {
        kTsteps += 1;
      }
      // mc_rtc::log::info("Generation by Ts_in , ksteps start = {}", kTsteps);
      while(kTsteps < t_steps_inputs_.size())
      {
        double nextTs = t_steps_inputs_[kTsteps];
        // mc_rtc::log::info("next {}",nextTs);
        if(nextTs - ti > Ts_min_)
        {
          double f = 1 / (nextTs - ti);
          double v = std::pow(f * robot_height_ / (0.157 * 172), 2) * 1e-2;
          // mc_rtc::log::info("speed v {}", v);
          for(int k = 0; k < (int)((nextTs - ti) / delta_); k++)
          {

            v_inputs_.push_back(sva::MotionVecd(Eigen::Vector3d::Zero() , Eigen::Vector3d{v,0,0} ));
            ref_traj_point next_P_traj = P_traj_.back();
            next_P_traj.ori_ += v_inputs_.back().angular().z() * delta_;
            next_P_traj.pose_ += sva::RotZ(-next_P_traj.ori()).block(0,0,2,2) * v_inputs_.back().linear() * delta_;
            for(int i = 0; i < 3; i++)
            {
              if(std::abs(next_P_traj.pose_(i)) < 1e-5)
              {
                next_P_traj.pose_(i) = 0;
              }
            }
            P_traj_.push_back(next_P_traj);
          }
          // mc_rtc::log::info("V size after t {}",Vx_.size());
        }
        ti = nextTs;
        kTsteps += 1;
      }
    }
   
    // 1-3 Generate a 0 velocity
   
    if(v_inputs_.size() != P_)
    {
      // mc_rtc::log::info("Generation by 0  vel");
      while(v_inputs_.size() < P_)
      {
        v_inputs_.push_back(sva::MotionVecd::Zero());
        P_traj_.push_back(P_traj_.back());
      }
    }
  }

  // mc_rtc::log::success("[Steps Timings Generation] 1 OK , V size {} , P_ {}",Vx_.size(),P_);

  // mc_rtc::log::info("Vx : {}", Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Vx_.data(), Vx_.size()));
  // mc_rtc::log::info("Vy : {}", Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Vy_.data(), Vy_.size()));
  // mc_rtc::log::info("omega : {}", Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Omega_.data(), Omega_.size()));

  // 2 Using V and Pf_in, generate Candidate Ts and link them to step inputs (if a ts is link to a pf then create its
  // boundaries)

  double ti = 0;
  int ki = 0;
  int kfoot = 0;
  std::vector<Eigen::Vector2d> StepsTimings_Upper_Lower_cstr;

  while(ti < Tp_)
  {
    // mc_rtc::log::info("speed at ki {} , {}",ki,Eigen::Vector2d{Vx_[ki], Vy_[ki]});
    if(v_inputs_[ki].linear().norm() > 1e-4)
    {
      Eigen::Vector2d v_i = v_inputs_[ki].linear().segment(0,2) * 1e2;
      double next_ts = (robot_height_ / (sqrt(v_i.norm()) * 0.157 * 172));
      mc_filter::utils::clampInPlace(next_ts, Ts_min_, Ts_max_);
      ti += next_ts;

      // mc_rtc::log::info("next t by deans {} at ki {}",ti,ki);
    }
    else
    {
      ti += Ts_;
      // mc_rtc::log::info("next t by mean {} at ki {}",ti,ki);
    }
    if(kfoot < steps_inputs_.size())
    {

      ref_traj_point P_ki = P_traj_[std::min(ki, (int)P_traj_.size() - 1)];

      double dist_i_x = std::abs((sva::RotZ(-P_ki.ori()).block(0,0,2,2).transpose() * ( P_ki.pose() - steps_inputs_[kfoot].pose())).x());
      double dist_i_y = std::abs((sva::RotZ(-P_ki.ori()).block(0,0,2,2).transpose() * ( P_ki.pose() - steps_inputs_[kfoot].pose())).x());
      double dist_i_theta = P_ki.ori() - steps_inputs_[kfoot].ori();
      if(std::abs(dist_i_theta) > M_PI)
      {
        dist_i_theta -= (dist_i_theta / std::abs(dist_i_theta)) * 2 * M_PI;
      }
      // mc_rtc::log::info("Foot dist x {} , y {} , theta {}" , dist_i_x , dist_i_y, dist_i_theta );
      if(dist_i_x < d_h_x / 2 && dist_i_y < l_ / 2 + d_h_y && std::abs(dist_i_theta) < max_theta)
      {
        int n = FootSteps_indx_.size();
        // mc_rtc::log::info("Step No: {} ",n+1);
        // mc_rtc::log::info("DIST Y : {} & Dtheta {} at ki : {}",dist_i_y,dist_i_theta,ki);
        // mc_rtc::log::info("Pki.z() : {}",P_ki.z());
        ki = Get_ki(ki, kfoot);
        // mc_rtc::log::info("ki : {}",ki);

        P_ki = P_traj_[ki];
        double dist_i = (P_ki.pose() - steps_inputs_[kfoot].pose()).norm(); 
        dist_i_y = (sva::RotZ(-P_ki.ori()).transpose().block(0,0,2,2) * (P_ki.pose() - steps_inputs_[kfoot].pose())).y();

        if(supportFoot == "RightFoot")
        {
          dist_i_y *= (2 * (n % 2) - 1);
        }
        else
        {
          dist_i_y *= (2 * ((n + 1) % 2) - 1);
        }
        // mc_rtc::log::info("diy {} ",dist_i_y  );
        if(dist_i_y > 0)
        {
          // mc_rtc::log::info("foot link kfoot {} at ki : {} ",kfoot,ki);
          // mc_rtc::log::info("DIST F : {}",dist_i);
          ti = ((double)ki) * delta_;
          FootSteps_indx_.push_back(kfoot);
          StepsTimings_Upper_Lower_cstr.push_back(Eigen::Vector2d{ti + 0.5 * Ts_min_, ti - 0.5 * Ts_min_});
          N_steps = FootSteps_indx_.size();
        }
        else
        {
          FootSteps_indx_.push_back(-1);
          StepsTimings_Upper_Lower_cstr.push_back(Eigen::Vector2d{-1, -1});
          FootSteps_indx_.push_back(kfoot);
          StepsTimings_Upper_Lower_cstr.push_back(Eigen::Vector2d{ti + 0.5 * Ts_min_, ti - 0.5 * Ts_min_});
          N_steps = FootSteps_indx_.size();
          StepsTimings_.push_back(ti);
          StepsTimings_indx_.push_back(ki);

          ti += Ts_;

          ki = (int)(ti / delta_);
          // mc_rtc::log::info("foot link kfoot for next step {} at ki : {} ",kfoot,ki);
        }
        kfoot += 1;
      }
      else
      {
        // mc_rtc::log::info("No foot close enough");
        FootSteps_indx_.push_back(-1);
        StepsTimings_Upper_Lower_cstr.push_back(Eigen::Vector2d{-1, -1});
      }
    }
    else
    {
      // mc_rtc::log::info("No foot to link");
      FootSteps_indx_.push_back(-1);
      StepsTimings_Upper_Lower_cstr.push_back(Eigen::Vector2d{-1, -1});
    }

    ki = (int)std::round(ti / delta_);
    StepsTimings_.push_back(ti);
    StepsTimings_indx_.push_back(ki);
  }
  // mc_rtc::log::success("[Steps Timings Generation] 2 OK, Ts size {} ; Ps cstr size {} ", StepsTimings_.size()
  // ,StepsTimings_Upper_Lower_cstr.size());
  while(StepsTimings_[StepsTimings_.size() - 1] > Tp_)
  {
    StepsTimings_.pop_back();
    FootSteps_indx_.pop_back();
    StepsTimings_indx_.pop_back();
  }
  // 3 Find the ts index to fit the input Ts

  Eigen::VectorXd Ts_Candidate =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(StepsTimings_.data(), StepsTimings_.size());

  Steps_timings_output cstr_Ts = Get_constrained_Ts(Ts_Candidate, StepsTimings_Upper_Lower_cstr);
  // mc_rtc::log::info("Qp has : {}", cstr_Ts.QPsuccess);

  for(int i = 0; i < cstr_Ts.Ts.size(); i++)
  {
    StepsTimings_[i] = cstr_Ts.Ts(i);
  }

  for(int i = 0; i < StepsTimings_.size(); i++)
  {
    double t_i = StepsTimings_[i];
    StepsTimings_indx_[i] = (int)std::round(t_i / delta_);
  }
  while(StepsTimings_[StepsTimings_.size() - 1] > Tp_ - delta_)
  {
    StepsTimings_.pop_back();
    FootSteps_indx_.pop_back();
    StepsTimings_indx_.pop_back();
  }

  // mc_rtc::log::success("[Steps Timings Generation] 3 OK");
  F_ = StepsTimings_.size();
  // mc_rtc::log::info("Params //");
  //   for (int k = 0; k<F_; k++){
  //     mc_rtc::log::info("Ts : {} ",StepsTimings_[k]);
  //     mc_rtc::log::info(StepsTimings_indx_[k]);
  //     mc_rtc::log::info("Foot: {} ",FootSteps_indx_[k]);
  //     mc_rtc::log::info("//");
  // }
  // mc_rtc::log::info("Vx : {}", Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Vx_.data(), Vx_.size()));
  // mc_rtc::log::info("Vy : {}", Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Vy_.data(), Vy_.size()));
  // mc_rtc::log::info("omega : {}", Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Omega_.data(), Omega_.size()));
  // for (int k = 0 ; k < 4 ; k ++){
  //   mc_rtc::log::info("Ptraj : {}" , P_traj_[k].z());
  // }

  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  double ProcessTime = time_span.count();
  // mc_rtc::log::success("StepTiming gened in: " + std::to_string(ProcessTime) + " ms");
}

Steps_timings_output FootStepGen::Get_constrained_Ts(const Eigen::VectorXd & Ts_candidate,
                                                     const std::vector<Eigen::Vector2d> & StepsTimings_Upper_Lower_cstr)
{
  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();

  int nTs = Ts_candidate.rows();
  int n_eq_cstr = t_steps_inputs_.size();
  Eigen::MatrixXd M = Eigen::MatrixXd::Identity(nTs, nTs);
  Steps_timings_output Output;
  Eigen::MatrixXd Delta = Eigen::MatrixXd::Identity(nTs, nTs);
  Eigen::VectorXd Upper_cstr = Eigen::VectorXd::Ones(nTs) * Ts_max_;
  Eigen::VectorXd Lower_cstr = Eigen::VectorXd::Ones(nTs) * Ts_min_;
  Aeq = Eigen::MatrixXd::Zero(n_eq_cstr, nTs);
  beq = Eigen::VectorXd::Zero(n_eq_cstr);

  for(int i = 1; i < nTs; i++)
  {
    if(StepsTimings_Upper_Lower_cstr[i].x() == -1)
    {
      Delta(i, i - 1) = -1;
    }
    else
    {
      Upper_cstr(i) = StepsTimings_Upper_Lower_cstr[i](0);
      Lower_cstr(i) = StepsTimings_Upper_Lower_cstr[i](1);
    }
  }
  for(int i = 0; i < t_steps_inputs_.size(); i++)
  {
    Aeq(i, i) = 1;
    beq(i) = t_steps_inputs_[i];
  }
  Aineq.resize(2 * Delta.rows(), nTs);
  bineq.resize(2 * Delta.rows());
  Aineq << Delta, -Delta;
  bineq << Upper_cstr, -Lower_cstr;
  Q_ = Eigen::MatrixXd::Identity(nTs, nTs) * 1e-12 + (M.transpose() * M);
  p_ = (-M.transpose() * Ts_candidate);
  // std::cout << Aeq << std::endl;
  // std::cout << beq << std::endl;
  // std::cout << Ts_candidate << std::endl;
  Output.Ts = solveQP();
  Output.QPsuccess = QPsuccess;
  Output.loss = (Output.Ts - Ts_candidate).norm();

  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  double ProcessTime = time_span.count();
  // mc_rtc::log::success("Time QP gened in: {}  ms", ProcessTime);

  return Output;
}

Footsteps_plan FootStepGen::compute_plan()
{
  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();

  GetStepsTimings();

  Theta_f_.resize(F_ + 1, 1);
  Theta_f_(0) = plan_.support_foot().ori();


  if(F_ <= 1)
  {
    mc_rtc::log::error_and_throw<std::runtime_error>("[Steps Generation] No step to compute");
    return plan_;
  }

  // Solving QP 1 for Orientation

  Eigen::VectorXd Dtheta_upper_lim(Eigen::VectorXd::Ones(F_) * max_theta);
  Eigen::VectorXd Dtheta_lower_lim(Eigen::VectorXd::Ones(F_) * -max_theta);
  Aeq = Eigen::MatrixXd::Zero(F_, F_);
  beq = Eigen::VectorXd::Zero(F_);
  Eigen::MatrixXd Delta(Eigen::MatrixXd::Identity(F_, F_));
  Eigen::VectorXd Dtheta_vel(Eigen::VectorXd::Zero(F_)); // Integration of omega between two steps

  for(int k = 0; k < F_; ++k)
  {
    // std::cout << "Theta k " << k << std::endl;

    if(k == 0)
    {
      Dtheta_vel(k) = P_traj_[StepsTimings_indx_[k]].ori() - Theta_f_(0) + Theta_f_(0);
      Dtheta_upper_lim(k) += Theta_f_(0);
      Dtheta_lower_lim(k) += Theta_f_(0);
    }

    else
    {

      Dtheta_vel(k) = P_traj_[StepsTimings_indx_[k]].ori() - P_traj_[StepsTimings_indx_[k - 1]].ori();

      if(std::abs(Dtheta_vel(k)) > M_PI)
      {
        Dtheta_vel(k) -= (Dtheta_vel(k) / std::abs(Dtheta_vel(k))) * 2 * M_PI;
      }

      Delta(k, k - 1) = -1;
    }

    if(FootSteps_indx_[k] != -1)
    { 
      double theta_cstr = steps_inputs_[FootSteps_indx_[k]].ori();
      if(IntegrateVelProfile(StepsTimings_indx_[k]).ori() - steps_inputs_[FootSteps_indx_[k]].ori() > M_PI)
      {
        theta_cstr -=
            steps_inputs_[FootSteps_indx_[k]].ori() / std::abs(steps_inputs_[FootSteps_indx_[k]].ori()) * 2 * M_PI;
      }
      Aeq(k, k) = 1;
      beq(k) = theta_cstr;
    }
  }

  Aineq.resize(2 * Delta.rows(), F_);
  bineq.resize(2 * Delta.rows());
  Aineq << Delta, -Delta;
  bineq << Dtheta_upper_lim, -Dtheta_lower_lim;
  // mc_rtc::log::info("theta Aineq {}",Aineq);
  // mc_rtc::log::info("theta bineq {}",bineq);

  Eigen::MatrixXd M = Eigen::MatrixXd::Identity(F_, F_);
  for(int i = 1; i < F_; i++)
  {
    M(i, i - 1) = -1;
  }
  Eigen::VectorXd b = Dtheta_vel;
  Q_ = Eigen::MatrixXd::Identity(F_, F_) * 1e-12 + (M.transpose() * M);
  p_ = (-M.transpose() * b);
  Eigen::VectorXd theta = solveQP();
  if(!QPsuccess)
  {
    mc_rtc::log::error("Step theta failed");
  }

  for(int k = 0; k < theta.size(); k++)
  {
    if(std::abs(theta(k)) < 1e-4)
    {
      theta(k) = 0;
    }
    if(std::abs(theta(k)) > M_PI)
    {
      Theta_f_(k + 1) = theta(k) - (theta(k) / std::abs(theta(k))) * 2 * M_PI;
    }
    else
    {
      Theta_f_(k + 1) = theta(k);
    }
  }

  // mc_rtc::log::info("Theta out {}",Theta_f_);

  // std::cout << "Theta" << std::endl;

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

  std::vector<Eigen::VectorXd> cstr_vec;
  std::vector<Eigen::MatrixXd> Normal_Vec;

  std::vector<Eigen::Vector2d> Dpos_vel(F_); // Integral of the speed between two steps

  Eigen::Matrix2d R_Theta_0;
  double theta_0 = P_traj_[0].ori();
  R_Theta_0 = sva::RotZ(-theta_0).block(0, 0, 2, 2);

  Dpos_vel[0] = (P_traj_[StepsTimings_indx_[0]].pose() - P_traj_[0].pose()) 
              + sgn * R_Theta_0 * Eigen::Vector2d{0, l_ / 2} 
              + plan_.support_foot().pose(); 

  Admissible_Region Kinematic_Rectangle(Eigen::Vector3d{0, 0, Theta_f_(0)}, Eigen::Vector3d{d_h_x, d_h_y, 0});

  Eigen::VectorXd bcstr = Kinematic_Rectangle.Offset
                          + Kinematic_Rectangle.Polygone_Normals.transpose()
                                * (plan_.support_foot().pose() + sgn * R_Theta_0 * Eigen::Vector2d{0, l_});

  Normal_Vec.push_back(Kinematic_Rectangle.Polygone_Normals.transpose());
  cstr_vec.push_back(bcstr);

  for(int k = 1; k < F_; k++)
  {

    double theta_k = P_traj_[StepsTimings_indx_[k - 1]].ori();

    Eigen::Vector3d Integral = P_traj_[StepsTimings_indx_[k]].vec3_pose() - P_traj_[StepsTimings_indx_[k - 1]].vec3_pose();
    Dpos_vel[k] = Eigen::Vector2d{Integral.x(), Integral.y()};

    Eigen::Matrix2d R_Theta_k;
    R_Theta_k = sva::RotZ(-theta_k).block(0, 0, 2, 2);

    Dpos_vel[k] += sgn * (1 - 2 * ((k % 2))) * R_Theta_k * Eigen::Vector2d{0, l_ / 2};

    Kinematic_Rectangle = Admissible_Region(Eigen::Vector3d{0, 0, Theta_f_(k)}, Eigen::Vector3d{d_h_x, d_h_y, 0});
    R_Theta_k = sva::RotZ(-Theta_f_(k - 1)).block(0, 0, 2, 2);

    bcstr = Kinematic_Rectangle.Offset
            + sgn * (1 - 2 * ((k % 2))) * Kinematic_Rectangle.Polygone_Normals.transpose() * R_Theta_k
                  * Eigen::Vector2d{0, l_};

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

  int N_cstr = 0;
  for(int k = 0; k < Normal_Vec.size(); k++)
  {
    N_cstr += (int)Normal_Vec[k].rows();
  }

  Aineq = Eigen::MatrixXd::Zero(N_cstr, 2 * F_);
  bineq = Eigen::VectorXd::Zero(N_cstr);

  int step_indx = 0;
  int cstr_index = 0;
  for(int i_ineq = 0; i_ineq < Normal_Vec.size(); i_ineq++)
  {

    Eigen::MatrixXd n_vec = Normal_Vec[i_ineq];
    Eigen::VectorXd ineq = cstr_vec[i_ineq];

    for(int cstr = 0; cstr < n_vec.rows(); cstr++)
    {

      Aineq(cstr_index + cstr, step_indx) = n_vec(cstr, 0);
      Aineq(cstr_index + cstr, step_indx + 1) = n_vec(cstr, 1);

      bineq(cstr_index + cstr) = ineq(cstr);
    }

    step_indx += 2;
    cstr_index += n_vec.rows();
  }

  Aineq = Aineq * Delta;
  M = Eigen::MatrixXd::Identity(2 * F_, 2 * F_);
  for(int i = 1; i < F_; i++)
  {
    M(2 * i, 2 * (i - 1)) = -1;
    M(2 * i + 1, 2 * (i - 1) + 1) = -1;
  }

  b = Eigen::VectorXd::Zero(2 * F_);
  for(int i = 0; i < F_; i++)
  {
    b(2 * i) = Dpos_vel[i].x();
    b(2 * i + 1) = Dpos_vel[i].y();
  }

  Q_ = Eigen::MatrixXd::Identity(2 * F_, 2 * F_) * 1e-12 + (M.transpose() * M);
  p_ = (-M.transpose() * b);
  // mc_rtc::log::info("Step Aineq {}",Aineq);
  // mc_rtc::log::info("Step bineq {}",bineq);
  Eigen::VectorXd XY(solveQP());
  if(!QPsuccess)
  {
    mc_rtc::log::error("Step QP failed");
  }

  // mc_rtc::log::info("Step out {}",XY);

  for(int k = 0; k < F_; k++)
  {
    double xf = XY(2 * k);
    double yf = XY(2 * k + 1);
    plan_.push_back(Footstep(Eigen::Vector2d{xf , yf} , Theta_f_(k),StepsTimings_[k],Eigen::Vector2d{0.1,0.1}));
  }
  // mc_rtc::log::success("Position OK");
  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  double ProcessTime = time_span.count();
  // mc_rtc::log::success("FootSteps computed in : {} ms",ProcessTime);
  return plan_;
}

std::vector<ref_traj_point> FootStepGen::GetRefTrajectory( ref_traj_point & P_s_0, ref_traj_point & P_s_1)
{

  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();

  Eigen::Vector2d init_pose = P_s_0.pose();
  Eigen::Vector2d init_ori = {cos(P_s_0.ori()), sin(P_s_0.ori())};
  Eigen::Vector2d target_pose = P_s_1.pose();
  Eigen::Vector2d target_ori = {cos(P_s_1.ori()), sin(P_s_1.ori())};

  std::vector<ref_traj_point> Output;

  path.reset(init_pose, init_ori, target_pose, target_ori);
  double traj_lenght = path.arcLength(0, 1);
  // mc_rtc::log::info("called p0 {} and p1 {} ; L_traj : {}",P_s_0,P_s_1,traj_lenght);
  double dl = v_ * delta_;
  int N = (int)std::round(traj_lenght / dl);

  if((init_pose - target_pose).norm() < 5e-2)
  {
    mc_rtc::log::warning("Short Trajectory");
    return std::vector<ref_traj_point>{P_s_0, P_s_1};
  }
  // N = 20;

  for(int k = 0; k < N + 1; k++)
  {
    double t = ((double)k) / ((double)N);
    Eigen::Vector2d pos_t = path.pos(t);
    // mc_rtc::log::info("t : {} ; p : {} ",t,pos_t);
    Eigen::Vector2d ori_t = path.tangent(t);
    double theta = atan2(ori_t.y(), ori_t.x());

    // Vx_.push_back( (pos_t.x() - P_traj_[P_traj_.size() - 1].x()) /delta_);
    // Vy_.push_back( (pos_t.y() - P_traj_[P_traj_.size() - 1].y()) /delta_);
    // Omega_.push_back( ( theta - P_traj_[P_traj_.size() - 1].z()) /delta_);

    Output.push_back(ref_traj_point( pos_t , theta));
  }

  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  double ProcessTime = time_span.count();
  mc_rtc::log::success("Ref trajectory computed in : {} at N : {} ", ProcessTime, N);

  return Output;
}

Eigen::VectorXd FootStepGen::solveQP()
{

  Eigen::QuadProgDense QP;

  int Nvar = Q_.rows();
  int NIneqConstr = Aineq.rows();
  int NEqConstr = Aeq.rows();
  // QP.tolerance(1e-4);
  QP.problem(Nvar, NEqConstr, NIneqConstr);
  QPsuccess = QP.solve(Q_, p_, Aeq, beq, Aineq, bineq);

  return QP.result();
}

ref_traj_point FootStepGen::IntegrateVelProfile(int k_end)
{

  int k_start = 0;
  if(k_end < P_traj_.size())
  {
    return P_traj_[k_end];
  }

  Eigen::Vector3d Output(plan_.support_foot().vec3_pose()); Output.z() = plan_.support_foot().ori();
  if(P_traj_.size() == 0)
  {

    if(supportFoot == "RightFoot")
    {
      Output += l_ / 2 * Eigen::Vector3d{-sin(Output.z()), cos(Output.z()), 0};
    }
    else
    {
      Output -= l_ / 2 * Eigen::Vector3d{-sin(Output.z()), cos(Output.z()), 0};
    }
  }
  else
  {
    Output = P_traj_.back().vec3_pose();
    k_start = P_traj_.size() - 1;
  }

  for(int i = k_start; i < std::min(k_end, (int) v_inputs_.size()); i++)
  {
    Output.z() += v_inputs_[i].angular().z() * delta_;
    Output += sva::RotZ(-Output.z()) * v_inputs_[i].linear() * delta_;
  }

  return ref_traj_point(Output.segment(0,2), Output.z());
}

int FootStepGen::Get_ki(int k, int kfoot)
{
  std::chrono::high_resolution_clock::time_point t_clock = std::chrono::high_resolution_clock::now();
  int ki = k;
  ref_traj_point P_ki(P_traj_[ki]);
  double dist_i = (P_ki.pose() - steps_inputs_[kfoot].pose()).norm();
  ref_traj_point P_kip1(P_traj_[ki + 1]);
  if(ki + 1 < P_traj_.size() - 1)
  {
    P_kip1 = P_traj_[ki + 1];
  }
  ref_traj_point P_kim1 = P_traj_[ki - 1];

  double dist_ip1 = (P_kip1.pose() - steps_inputs_[kfoot].pose()).norm(); 
  double dist_im1 = (P_kim1.pose() - steps_inputs_[kfoot].pose()).norm(); 
  if(dist_ip1 < dist_im1)
  {
    while(dist_ip1 < dist_i && (double)(ki + 1 - k) * delta_ < 0.5 * Ts_min_)
    {
      ki += 1;
      dist_i = dist_ip1;
      if(ki + 1 > (int)(Tp_ / delta_))
      {
        break;
      }
      P_kip1 = P_traj_[std::min(ki + 1, (int)P_traj_.size() - 1)];
      dist_ip1 = (P_kip1.pose() - steps_inputs_[kfoot].pose()).norm(); 
    }
  }
  else
  {
    while(dist_im1 < dist_i && (double)(k - (ki - 1)) * delta_ < Ts_min_ && (double)(ki - 1) * delta_ > Ts_min_)
    {
      ki -= 1;
      dist_i = dist_im1;
      P_kim1 = P_traj_[std::min(ki - 1, (int)P_traj_.size() - 1)];
      dist_im1 = (P_kim1.pose() - steps_inputs_[kfoot].pose()).norm(); 
    }
  }
  std::chrono::duration<double, std::milli> time_span = std::chrono::high_resolution_clock::now() - t_clock;
  double ProcessTime = time_span.count();
  // mc_rtc::log::success("Get ki time: " + std::to_string(ProcessTime) + " ms");
  return ki;
}

} //namespace footsteps_planner
} //namespace mc_plugin