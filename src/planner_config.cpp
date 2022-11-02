#include "footsteps_planner.h"
namespace mc_plugin
{
namespace footsteps_planner
{

FootStepGen::FootStepGen()
{
  Ts_ = (Ts_max_ + Ts_min_) / 2;
  P_ = (int)(Tp_ / delta_);
}
FootStepGen::FootStepGen(const mc_rtc::Configuration & config)
{
  reconfigure(config);
}

void FootStepGen::reconfigure(const mc_rtc::Configuration & config)
{
  mc_rtc::log::info("footsteps_planner::init called with configuration:\n{}", config.dump(true, true));
  if(config.has("Ts_limit"))
  {
    Eigen::Vector2d Ts_range = config("Ts_limit");
    Ts_min_ = Ts_range(0);
    Ts_max_ = Ts_range(1);
  }
  if(config.has("feet_distance"))
  {
    l_ = config("feet_distance");
  }
  if(config.has("Tp"))
  {
    Tp_ = config("Tp");
  }
  if(config.has("delta"))
  {
    delta_ = config("delta");
  }
  if(config.has("kinematics_cstr"))
  {
    Eigen::Vector2d kin_cstr = config("kinematics_cstr");
    d_h_x = kin_cstr.x();
    d_h_y = kin_cstr.y();
  }
  if(config.has("mean_speed"))
  {
    v_ = config("mean_speed");
  }
  if(config.has("robot_height"))
  {
    robot_height_ = config("robot_height");
  }
  if(config.has("max_rotation"))
  {
    max_theta = config("max_rotation");
  }
  if(config.has("offset_angle_deg"))
  {
    theta_offset_ = config("offset_angle_deg");
    theta_offset_ *= mc_rtc::constants::PI / 180.;
  }


  Ts_ = (Ts_max_ + Ts_min_) / 2;
  P_ = (int)(Tp_ / delta_);
}

} // namespace footsteps_planner
} // namespace mc_plugin
