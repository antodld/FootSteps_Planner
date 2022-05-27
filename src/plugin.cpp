#include "plugin.h"

namespace mc_plugin
{

footsteps_planner_plugin::~footsteps_planner_plugin() = default;

void footsteps_planner_plugin::init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config)
{
  controller.controller().datastore().make<mc_rtc::Configuration>("footsteps_planner::planner_config");
  controller.controller().datastore().make<std::vector<sva::MotionVecd>>("footsteps_planner::input_vel");
  controller.controller().datastore().make<std::vector<sva::PTransformd>>("footsteps_planner::input_steps");
  controller.controller().datastore().make<std::vector<double>>("footsteps_planner::input_time_steps");
  controller.controller().datastore().make<std::string>("footsteps_planner::support_foot_name");
  controller.controller().datastore().make<sva::PTransformd>("footsteps_planner::support_foot_pose");
  controller.controller().datastore().make<std::vector<sva::PTransformd>>("footsteps_planner::output_steps");
  controller.controller().datastore().make<std::vector<double>>("footsteps_planner::output_time_steps");

  controller.controller().datastore().make_call("footstep_planner::compute_plan" , [this](mc_rtc::DataStore* datastore){compute_footsteps_plan(datastore);});

  planner_ = mc_plugin::footsteps_planner::FootStepGen(config);

  mc_rtc::log::info("footsteps_planner::init called with configuration:\n{}", config.dump(true, true));
}

void footsteps_planner_plugin::reset(mc_control::MCGlobalController & controller)
{
  mc_rtc::log::info("footsteps_planner::reset called");
}

void footsteps_planner_plugin::before(mc_control::MCGlobalController &)
{
  mc_rtc::log::info("footsteps_planner::before");
}

void footsteps_planner_plugin::after(mc_control::MCGlobalController & controller)
{
  mc_rtc::log::info("footsteps_planner::after");
}

void footsteps_planner_plugin::compute_footsteps_plan(mc_rtc::DataStore* datastore)
{

  mc_plugin::footsteps_planner::Footstep support_footstep(support_foot_pose_,0,Eigen::Vector2d::Ones() * 0.1);
  std::vector<mc_plugin::footsteps_planner::Footstep> input_footsteps;
  for (int k = 0 ; k < input_steps_.size() ; k++)
  {
    input_footsteps.push_back(mc_plugin::footsteps_planner::Footstep(input_steps_[k],0,Eigen::Vector2d::Ones() * 0.1));
  }
  planner_.Init(support_foot_name_,
                support_footstep,
                input_v_,
                input_t_steps_,
                input_footsteps);
  
  planner_.compute_plan();

  datastore->assign<std::vector<sva::PTransformd>>("footsteps_planner::output_steps",planner_.footsteps_plan().steps_PTpose());
  datastore->assign<std::vector<double>>("footsteps_planner::output_time_steps",planner_.footsteps_plan().steps_timings());

}

void footsteps_planner_plugin::gui(mc_control::MCGlobalController & controller)
{
  
  controller.controller().gui()->addElement(
      {"footstep"},
      mc_rtc::gui::Trajectory("Trajectory", mc_rtc::gui::Color(1., 1., 0.),
                              [this]() -> std::vector<Eigen::Vector3d>  { return planner_.Ref_Traj(); }),

      mc_rtc::gui::Polygon(
      "Steps", mc_rtc::gui::Color(1., 0.3, 0.),
      [this]() -> const std::vector<std::vector<Eigen::Vector3d>>  { return this->planner_.footsteps_plan().get_steps_corners(); })                      
                              
      );


}

mc_control::GlobalPlugin::GlobalPluginConfiguration footsteps_planner_plugin::configuration()
{
  mc_control::GlobalPlugin::GlobalPluginConfiguration out;
  out.should_run_before = false;
  out.should_run_after = false;
  out.should_always_run = true;
  return out;
}
} // namespace mc_plugin

EXPORT_MC_RTC_PLUGIN("footsteps_planner_plugin", mc_plugin::footsteps_planner_plugin)
