/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>


#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set the number of particles.
  num_particles = 100;
  weights.resize(num_particles, 0.);
  particles.resize(num_particles);

  // Initialize all particles to first position (based on estimates of
  //   x, y, theta and their uncertainties from GPS) and all weights to 1.
  // Add random Gaussian noise to each particle.
  default_random_engine gen;
  normal_distribution<double> N_x(0, std[0]);
  normal_distribution<double> N_y(0, std[1]);
  normal_distribution<double> N_theta(0, std[2]);
  for (int i=0; i<num_particles; i++) {
    weights[i] = 1.;

    Particle p;
    p.id = i;
    p.x = x + N_x(gen);
    p.y = y + N_y(gen);
    p.theta = theta + N_theta(gen);
    p.weight = weights[i];

    particles[i] = p;
  }
  is_initialized = true;
}

// Add measurements to each particle and add random Gaussian noise.
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  random_device rd;
  default_random_engine gen(rd());
  normal_distribution<float> N_x(0, std_pos[0]);
  normal_distribution<float> N_y(0, std_pos[1]);
  normal_distribution<float> N_theta(0, std_pos[2]);

  Particle p;
  // optimize calculation outside the loop
  float vy = velocity / yaw_rate;
  float tdt = yaw_rate * delta_t;
  // loop over all particles
  for (int i=0; i<num_particles; i++) {
    p = particles[i];

    float thetaold = p.theta;
    // do theta prediction first to add noise. then that noise affects x and y
    // assuming noise is not proportional to delta_t. not clear from instructions if it should be
    p.theta += tdt                             + N_theta(gen);
    p.x += vy * (sin(p.theta) - sin(thetaold)) + N_x(gen);
    p.y += vy * (cos(thetaold) - cos(p.theta)) + N_y(gen);

    particles[i] = p; // overwrite particle with new state
  }

}



// x_map and y_map will have x and y translated to map coordinates
void CarToMap(double x, double y, double x_car_map, double y_car_map, double phi_car_map, double& x_map, double & y_map)
{
  // using http://planning.cs.uiuc.edu/node99.html
  x_map = x * cos(phi_car_map) + y * sin(phi_car_map) + x_car_map;
  y_map = -x * sin(phi_car_map) + y * cos(phi_car_map) + y_car_map;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
  //   implement this method and use it as a helper during the updateWeights phase.

}

// Update the weights of each particle using a mult-variate Gaussian distribution
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations, Map map_landmarks) {

  // loop over all particles
  Particle p;
  for (int i=0; i<num_particles; i++) {
    p = particles[i];

    // translate observations from particle coordinates to map coordinates
    std::vector<LandmarkObs> observations_map_coord; // observations in map coordinates
    for (int j=0; j<observations.size(); j++)
    {
      LandmarkObs obs = observations[j]; // copy observation as new object
      CarToMap(obs.x, obs.y, p.x, p.y, p.theta, obs.x, obs.y); // we are overwriting obs.x and obs.y inside function
      observations_map_coord.push_back(obs);
    }

    // associate observations with landmarks
    // overview of the algorithm:
    // 1. find all map landmarks withing sensor_range around particle, save in close_landmarks
    // 2. build distance matrix (obs -> close_landmark)
    // 3. use Hungarian algorithm to associate observation to landmarks
    std::vector<LandmarkObs> close_landmarks;
    for (int i=0; i<map_landmarks.landmark_list.size(); i++)
    {
      Map::single_landmark_s l = map_landmarks.landmark_list[i];
      double x_diff = l.x_f - p.x;
      double y_diff = l.y_f - p.y;
      double d = sqrt(x_diff*x_diff + y_diff*y_diff);
      if (d<=sensor_range) {
        LandmarkObs obs; // new object
        obs.x = l.x_f;
        obs.y = l.y_f;
        close_landmarks.push_back(obs);
      }
    }

    //dataAssociation()

    // calculate probabilities
    float prob = 1.0; // probability

    particles[i] = p;
  }

}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::write(std::string filename) {
  // You don't need to modify this file.
  std::ofstream dataFile;
  dataFile.open(filename, std::ios::app);
  for (int i = 0; i < num_particles; ++i) {
    dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
  }
  dataFile.close();
}
