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

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
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
    // do it first to add noise. then that noise affects x and y
    p.theta += tdt                                  + N_theta(gen);
    p.x += vy * (sin(p.theta) - sin(thetaold)) + N_x(gen);
    p.y += vy * (cos(thetaold) - cos(p.theta)) + N_y(gen);

    particles[i] = p;
  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.

}

// x_map and y_map will have x and y translated to map coordinates
void CarToMap(double x, double y, double x_car_map, double y_car_map, double phi_car_map, double& x_map, double & y_map)
{
  // using http://planning.cs.uiuc.edu/node99.html
  x_map = x * cos(phi_car_map) - y * sin(phi_car_map) + x_car_map;
  y_map = x * sin(phi_car_map) + y * cos(phi_car_map) + y_car_map;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
    std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution

  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation 
  //   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
  //   for the fact that the map's y-axis actually points downwards.)
  //   http://planning.cs.uiuc.edu/node99.html

  // loop over all particles
  Particle p;
  for (int i=0; i<num_particles; i++) {
    p = particles[i];
    // translate observations from particle coordinates to map coordinates
    std::vector<LandmarkObs> observations_map;
    float prob = 1.0; // probability
    for (int j=0; j<observations.size(); j++)
    {
      LandmarkObs obs = observations[j];
      CarToMap(obs.x, obs.y, p.x, p.y, p.theta, obs.x, obs.y); // we are overwriting obs.x and obs.y inside function
      observations_map.push_back(obs);
    }

    float thetaold = p.theta;
    // do it first to add noise. then that noise affects x and y
//    p.theta += tdt                                  + N_theta(gen);
//    p.x += vy * (sin(p.theta) - sin(thetaold)) + N_x(gen);
//    p.y += vy * (cos(thetaold) - cos(p.theta)) + N_y(gen);

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
