cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME footsteps_planner_plugin)
set(PROJECT_DESCRIPTION "")
set(PROJECT_URL "")
set(PROJECT_VERSION 0.0.1)
set(PROJECT_USE_CMAKE_EXPORT TRUE)
set(CMAKE_CXX_STANDARD 11)

option(AUTOLOAD_${PROJECT_NAME}_PLUGIN "Automatically load ${PROJECT_NAME} plugin" OFF)

project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})

# Find mc_rtc if we are not building inside mc_rtc itself
if(NOT TARGET mc_rtc::mc_control)
  find_package(mc_rtc REQUIRED)
endif()

find_package(lipm_walking_controller REQUIRED)
include_directories(${lipm_walking_controller_DIR}/../../../include)

# Configure and install the plugin configuration
configure_file(etc/footsteps_planner_plugin.in.yaml "${CMAKE_CURRENT_BINARY_DIR}/etc/footsteps_planner_plugin.yaml")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/etc/footsteps_planner_plugin.yaml" DESTINATION "${MC_PLUGINS_RUNTIME_INSTALL_PREFIX}/etc")



set(plugin_SRC
  src/plugin.cpp
  src/footsteps_planner.cpp
  src/planner_config.cpp
  
  )

set(plugin_HDR
  src/plugin.h
  src/footsteps_planner.h
)

# The add_plugin macro is provided by mc_rtc
add_plugin(${PROJECT_NAME} "${plugin_SRC}" "${plugin_HDR}")
target_link_libraries(${PROJECT_NAME} PUBLIC mc_rtc::mc_control)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME footsteps_planner_plugin)







