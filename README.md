# Overview
Particle Filter project for the Localization course in Udacity's Self-Driving Car Nanodegree.

## Project Introduction
Your robot has been kidnapped and transported to a new location! 
Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, 
and lots of (noisy) sensor and control data.

In this project you will implement a 2 dimensional particle filter in C++. 
Your particle filter will be given a map and some initial localization information 
(analogous to what a GPS would provide). At each time step your filter will also get 
observation and control data. 

## Running the Code
Once you have this repository on your machine, `cd` into the repository's root directory and run the following commands from the command line:

```
> ./clean.sh
> ./build.sh
> ./run.sh
```

If everything worked you should see something like the following output:

```
Time step: 2444
Cumulative mean weighted error: x .1 y .1 yaw .02
Runtime (sec): 38.187226
Success! Your particle filter passed!

```

## Visualization
Provided notebook [data-exploration.ipynb](./data-exploration.ipynb) visualises the data
and creates the following animation of how filter works (click to see animation):

[![Particle Filter Visualisation](./particle_filter.png)](./particle_animation.mp4)


# Implementing the Particle Filter
The directory structure of this repository is as follows:

```
root
|   build.sh
|   clean.sh
|   CMakeLists.txt
|   README.md
|   run.sh
|
|___data
|   |   control_data.txt
|   |   gt_data.txt
|   |   map_data.txt
|   |
|   |___observation
|       |   observations_000001.txt
|       |   ... 
|       |   observations_002444.txt
|   
|___src
    |   helper_functions.h
    |   main.cpp
    |   map.h
    |   particle_filter.cpp
    |   particle_filter.h
```

## Inputs to the Particle Filter
You can find the inputs to the particle filter in the `data` directory. 

#### The Map*
`map_data.txt` includes the position of landmarks (in meters) on an arbitrary Cartesian coordinate system. Each row has three columns
1. x position
2. y position
3. landmark id

> * Map data provided by 3D Mapping Solutions GmbH.


#### Control Data
`control_data.txt` contains rows of control data. Each row corresponds to the control data for the corresponding time step. The two columns represent
1. vehicle speed (in meters per second)
2. vehicle yaw rate (in radians per second)

#### Observation Data
The `observation` directory includes around 2000 files. Each file is numbered according to the timestep in which that observation takes place. 

These files contain observation data for all "observable" landmarks. Here observable means the landmark is sufficiently close to the vehicle. Each row in these files corresponds to a single landmark. The two columns represent:
1. x distance to the landmark in meters (right is positive) RELATIVE TO THE VEHICLE. 
2. y distance to the landmark in meters (forward is positive) RELATIVE TO THE VEHICLE.

