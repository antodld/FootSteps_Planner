#include "plugin.h"

namespace mc_plugin
{

footsteps_planner_plugin::~footsteps_planner_plugin() = default;

void footsteps_planner_plugin::init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config)
{
  std::cout << "[footsteps_planner_plugin] initialize" << std::endl;
  controller.controller().datastore().make<mc_rtc::Configuration>("footsteps_planner::planner_config");
  controller.controller().datastore().make<std::vector<sva::MotionVecd>>("footsteps_planner::input_vel");
  controller.controller().datastore().make<std::vector<sva::PTransformd>>("footsteps_planner::input_steps");
  controller.controller().datastore().make<std::vector<double>>("footsteps_planner::input_time_steps");
  controller.controller().datastore().make<std::string>("footsteps_planner::support_foot_name");
  controller.controller().datastore().make<sva::PTransformd>("footsteps_planner::support_foot_pose");
  controller.controller().datastore().make<std::vector<sva::PTransformd>>("footsteps_planner::output_steps");
  controller.controller().datastore().make<std::vector<double>>("footsteps_planner::output_time_steps");

  auto & ctl = controller.controller();
  controller.controller().datastore().make_call("footstep_planner::compute_plan",
                                                [this, &ctl]() { compute_footsteps_plan(ctl); });
  controller.controller().datastore().make_call(
      "footstep_planner::configure",
      [this](const mc_rtc::Configuration & config) { planner_ = mc_plugin::footsteps_planner::FootStepGen(config); });

  if(controller.controller().config().has("footsteps_planner"))
  {
    planner_ = mc_plugin::footsteps_planner::FootStepGen(controller.controller().config()("footsteps_planner"));
  }
  else
  {
    planner_ = mc_plugin::footsteps_planner::FootStepGen(config);
  }

  gui(controller);
}

void footsteps_planner_plugin::reset(mc_control::MCGlobalController & controller)
{
  mc_rtc::log::info("footsteps_planner::reset called");
}

void footsteps_planner_plugin::before(mc_control::MCGlobalController & controller)
{
  mc_rtc::log::info("footsteps_planner::before");
}

void footsteps_planner_plugin::after(mc_control::MCGlobalController & controller)
{
  mc_rtc::log::info("footsteps_planner::after");
}

void footsteps_planner_plugin::compute_footsteps_plan(mc_control::MCController & controller)
{

  auto & datastore = controller.datastore();
  auto & support_foot_pose = datastore.get<sva::PTransformd>("footsteps_planner::support_foot_pose");
  auto & input_footsteps_pose = datastore.get<std::vector<sva::PTransformd>>("footsteps_planner::input_steps");
  auto & input_v = datastore.get<std::vector<sva::MotionVecd>>("footsteps_planner::input_vel");
  auto & input_t_steps = datastore.get<std::vector<double>>("footsteps_planner::input_time_steps");
  auto & support_foot_name = datastore.get<std::string>("footsteps_planner::support_foot_name");

  mc_plugin::footsteps_planner::Footstep support_footstep(support_foot_pose, 0, Eigen::Vector2d::Ones() * 0.1);
  std::vector<mc_plugin::footsteps_planner::Footstep> input_footsteps;

  for(int k = 0; k < input_footsteps_pose.size(); k++)
  {
    input_footsteps.push_back(
        mc_plugin::footsteps_planner::Footstep(input_footsteps_pose[k], 0, Eigen::Vector2d::Ones() * 0.1));
  }

  planner_.init(support_foot_name, support_footstep, input_v, input_t_steps, input_footsteps);

  planner_.compute_plan();

  datastore.assign<std::vector<sva::PTransformd>>("footsteps_planner::output_steps",
                                                  planner_.footsteps_plan().steps_PTpose());
  datastore.assign<std::vector<double>>("footsteps_planner::output_time_steps",
                                        planner_.footsteps_plan().steps_timings());
}

void footsteps_planner_plugin::gui(mc_control::MCGlobalController & controller)
{

  controller.controller().gui()->addElement(
      {"Footsteps Planner"},
      mc_rtc::gui::Trajectory("Trajectory", mc_rtc::gui::Color(1., 1., 0.),
                              [this]() -> std::vector<Eigen::Vector3d> { return planner_.Ref_Traj(); })

      // mc_rtc::gui::Polygon("Steps", mc_rtc::gui::Color(0., 1., 0.),
      //                      [this]() -> std::vector<std::vector<Eigen::Vector3d>> {
      //                        return this->planner_.footsteps_plan().get_steps_corners();
      //                      })

      // mc_rtc::gui::Point3D(
      // "Steps", mc_rtc::gui::Color(0., 1., 0.),
      // [this]() -> std::vector<Eigen::Vector3d>  { return this->planner_.footsteps_plan().steps_pose(); })

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
