mc_rtc Foot steps planner plugin
==

This plugin generate velocity based footsteps for walking pattern generation 

## Dependencie

[mc_rtc](https://github.com/jrl-umi3218/mc_rtc)

## Installation

```bash
cd Footsteps_Planner
mkdir build && cd build
cmake ..
make 
sudo make install
```
## Using the Foot step Planner

Add to your controller configuration file

```
Plugins: [footsteps_planner_plugin]
```
