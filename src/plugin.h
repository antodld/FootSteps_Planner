/*
 * Copyright 2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once
#include <mc_control/GlobalPlugin.h>
#include <mc_control/GlobalPluginMacros.h>
#include <mc_control/mc_controller.h>
#include "footsteps_planner.h"
#include <vector>

namespace mc_plugin
{

struct footsteps_planner_plugin : public mc_control::GlobalPlugin
{
  void init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config) override;

  void reset(mc_control::MCGlobalController & controller) override;

  void before(mc_control::MCGlobalController &) override;

  void after(mc_control::MCGlobalController & controller) override;

  void compute_footsteps_plan(mc_control::MCController & controller);

  void gui(mc_control::MCGlobalController & controller);

  mc_control::GlobalPlugin::GlobalPluginConfiguration configuration() override;

  ~footsteps_planner_plugin() override;

private:
  // std::vector<sva::MotionVecd> input_v_;
  // std::vector<sva::PTransformd> input_steps_;
  // std::vector<double> input_t_steps_;
  // std::string support_foot_name_;
  // sva::PTransformd support_foot_pose_;

  bool centered_ref_trajectory_ = false;

  footsteps_planner::FootStepGen planner_;

  std::vector<sva::PTransformd> output_steps_;
  std::vector<double> output_t_steps_;
};

} // namespace mc_plugin
