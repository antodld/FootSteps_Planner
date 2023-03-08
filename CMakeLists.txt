cmake_minimum_required(VERSION 3.1)

set(PROJECT_NAME footsteps_planner_plugin)
set(PROJECT_DESCRIPTION "")
set(PROJECT_URL "")
set(PROJECT_VERSION 0.0.1)
set(PROJECT_USE_CMAKE_EXPORT FALSE)
set(CMAKE_CXX_STANDARD 11)

include(cmake/base.cmake)

project(${PROJECT_NAME} LANGUAGES CXX VERSION ${PROJECT_VERSION})
option(AUTOLOAD_${PROJECT_NAME}_PLUGIN "Automatically load ${PROJECT_NAME} plugin" OFF)

# Find mc_rtc if we are not building inside mc_rtc itself
if(NOT TARGET mc_rtc::mc_control)
  find_package(mc_rtc REQUIRED)
endif()

# Configure and install the plugin configuration
configure_file(etc/${PROJECT_NAME}.in.yaml "${CMAKE_CURRENT_BINARY_DIR}/etc/${PROJECT_NAME}.yaml")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/etc/${PROJECT_NAME}.yaml" DESTINATION "${MC_PLUGINS_RUNTIME_INSTALL_PREFIX}/etc")

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






