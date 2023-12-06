#include "plugin.h"

namespace mc_plugin
{

footsteps_planner_plugin::~footsteps_planner_plugin() = default;

void footsteps_planner_plugin::init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config)
{
  config_.load(config);
  reset(controller);
}

void footsteps_planner_plugin::reset(mc_control::MCGlobalController & controller)
{

  controller.controller().datastore().make<std::vector<sva::MotionVecd>>("footsteps_planner::input_vel");
  controller.controller().datastore().make<std::vector<sva::PTransformd>>("footsteps_planner::input_ref_pose");
  controller.controller().datastore().make<std::vector<double>>("footsteps_planner::input_time_steps");
  controller.controller().datastore().make<std::string>("footsteps_planner::support_foot_name");
  controller.controller().datastore().make<sva::PTransformd>("footsteps_planner::support_foot_pose");
  controller.controller().datastore().make<std::vector<sva::PTransformd>>("footsteps_planner::output_steps");
  controller.controller().datastore().make<std::vector<double>>("footsteps_planner::output_time_steps");

  auto & ctl = controller.controller();
  controller.controller().datastore().make_call("footsteps_planner::compute_plan",
                                                [this, &ctl]() { compute_footsteps_plan(ctl); });
  controller.controller().datastore().make_call(
      "footsteps_planner::configure",
      [this](const mc_rtc::Configuration & config) {config_.load(config) ; planner_ = mc_plugin::footsteps_planner::FootStepGen(config_); });

  if(controller.controller().config().has("footsteps_planner"))
  {
    config_.load(controller.controller().config()("footsteps_planner"));
    if(controller.controller().config()("footsteps_planner").has("centered_trajectory"))
    {
      centered_ref_trajectory_ = controller.controller().config()("footsteps_planner")("centered_trajectory");
    }
  }

  planner_ = mc_plugin::footsteps_planner::FootStepGen(config_);
  
  gui(controller);
}

void footsteps_planner_plugin::compute_footsteps_plan(mc_control::MCController & controller)
{

  auto & datastore = controller.datastore();
  auto support_foot_pose = datastore.get<sva::PTransformd>("footsteps_planner::support_foot_pose");
  auto ref_pose = datastore.get<std::vector<sva::PTransformd>>("footsteps_planner::input_ref_pose");
  auto input_v = datastore.get<std::vector<sva::MotionVecd>>("footsteps_planner::input_vel");
  auto input_t_steps = datastore.get<std::vector<double>>("footsteps_planner::input_time_steps");
  auto support_foot_name = datastore.get<std::string>("footsteps_planner::support_foot_name");

  const Eigen::Vector2d size = Eigen::Vector2d::Ones() * 0.1;
  mc_plugin::footsteps_planner::Footstep support_footstep(support_foot_pose, 0, size);
  std::vector<mc_plugin::footsteps_planner::Footstep> input_ref_pose;

  for(size_t k = 0; k < ref_pose.size(); k++)
  {
    input_ref_pose.push_back(
        mc_plugin::footsteps_planner::Footstep(ref_pose[k], 0, size));
  }
  
  planner_.init(support_foot_name, support_footstep, input_v, input_t_steps, input_ref_pose);

  planner_.compute_plan();

  datastore.assign<std::vector<sva::PTransformd>>("footsteps_planner::output_steps",
                                                  planner_.footsteps_plan().steps_PTpose());
  datastore.assign<std::vector<double>>("footsteps_planner::output_time_steps",
                                        planner_.footsteps_plan().steps_timings());
}

void footsteps_planner_plugin::gui(mc_control::MCGlobalController & controller)
{
  controller.controller().gui()->addElement(
      {"Footsteps Planner", "Configuration"},
      mc_rtc::gui::Form(
          "Configure", [this](const mc_rtc::Configuration & conf) { planner_.reconfigure(conf); },
          mc_rtc::gui::FormArrayInput("Ts_limit", false,
                                      [this]() -> std::array<double, 2> {
                                        return {planner_.Ts_min_, planner_.Ts_max_};
                                      }),
          mc_rtc::gui::FormNumberInput("feet_distance", false, [this]() { return planner_.l_; }),
          mc_rtc::gui::FormNumberInput("Tp", false, [this]() { return planner_.Tp_; }),
          mc_rtc::gui::FormNumberInput("delta", false, [this]() { return planner_.delta_; }),
          mc_rtc::gui::FormArrayInput("kinematics_cstr", false,
                                      [this]() -> std::array<double, 2> {
                                        return {planner_.d_h_x, planner_.d_h_y};
                                      }),
          mc_rtc::gui::FormNumberInput("mean_speed", false, [this]() { return planner_.v_; }),
          mc_rtc::gui::FormNumberInput("robot_height", false, [this]() { return planner_.robot_height_; }),
          mc_rtc::gui::FormNumberInput("max_rotation", false, [this]() { return planner_.max_theta; }),
          mc_rtc::gui::FormNumberInput("offset_angle_deg", false,
                                       [this]() { return planner_.theta_offset_ * 180 / mc_rtc::constants::PI; })));


  controller.controller().gui()->addElement(
      {"Footsteps Planner","Visuals"}, mc_rtc::gui::Trajectory("Trajectory", mc_rtc::gui::Color(1., 1., 0.),
                                                    [this]() -> std::vector<Eigen::Vector3d> {
                                                      return planner_.ref_traj(centered_ref_trajectory_);
                                                    }),
      mc_rtc::gui::Point3D(
      "Intersec", mc_rtc::gui::Color(0., 1., 0.),
      [this]() -> Eigen::Vector3d { return Eigen::Vector3d{this->planner_.intersec.x(),this->planner_.intersec.y(),0}; } ,
        [this](const Eigen::Vector3d & p) {this->planner_.intersec = p.segment(0,2);}),

      mc_rtc::gui::Label("r",[this]() -> double {return this->planner_.r.norm();}),

      mc_rtc::gui::Arrow("Dist",[this]() -> Eigen::Vector3d {Eigen::Vector3d p = Eigen::Vector3d::Zero(); p.segment(0,2) = this->planner_.intersec - this->planner_.r; return p;},
                                [this]() -> Eigen::Vector3d {Eigen::Vector3d p = Eigen::Vector3d::Zero(); p.segment(0,2) = this->planner_.intersec; return p;})

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
