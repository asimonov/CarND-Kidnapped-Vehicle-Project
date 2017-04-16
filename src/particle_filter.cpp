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
#include <map>


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

// Apply movement equation to each particle and add random Gaussian noise.
void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  random_device rd;
  default_random_engine gen(rd());
  normal_distribution<double> N_x(0, std_pos[0]);
  normal_distribution<double> N_y(0, std_pos[1]);
  normal_distribution<double> N_theta(0, std_pos[2]);

  Particle p;
  // optimize calculation outside the loop
  double vy = velocity / yaw_rate;
  double tdt = yaw_rate * delta_t;
  // loop over all particles
  for (int i=0; i<num_particles; i++) {
    p = particles[i];

    double thetanew = p.theta + tdt;
    if (fabs(yaw_rate)>1e-6) {
      p.x += vy * (sin(thetanew) - sin(p.theta)) + N_x(gen);
      p.y += vy * (cos(p.theta) - cos(thetanew)) + N_y(gen);
    } else {
      p.x += velocity * delta_t * cos(p.theta) + N_x(gen);
      p.y += velocity * delta_t * sin(p.theta) + N_y(gen);
    }
    p.theta = thetanew + N_theta(gen);

    particles[i] = p; // overwrite particle with new state
  }

}




void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
// unused. all calculations are done in updateWeights
}


// x_map and y_map will have x and y translated to map coordinates
void CarToMap(double x, double y, double x_car_map, double y_car_map, double phi_car_map, double& x_map, double & y_map)
{
  // using http://planning.cs.uiuc.edu/node99.html
  // Assuming map x axis points right, y axis points up.
  // Assuming car x axis points forward, y axis points left.
  x_map = x * cos(phi_car_map) - y * sin(phi_car_map) + x_car_map;
  y_map = x * sin(phi_car_map) + y * cos(phi_car_map) + y_car_map;
}


// calculate 2D-gaussian probability density with 0 correlation
double Gaussian2DNoCorrelation(double x, double mu_x, double y, double mu_y, double sigma_x, double sigma_y)
{
  return (1./(2.*M_PI*sigma_x*sigma_y)) * exp( - ( pow(x-mu_x, 2)/(2*sigma_x*sigma_x) + pow(y-mu_y, 2)/(2*sigma_y*sigma_y) ) );
}


// calculate distance between 2 LandmarkObs objects
double ObservationsDistance(LandmarkObs& obs1, LandmarkObs& obs2)
{
  return sqrt( pow(obs1.x - obs2.x, 2) + pow(obs1.y - obs2.y, 2) );
}


// Update the weights of each particle using a multi-variate Gaussian distribution
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations, Map map_landmarks) {

  // if there are no observations do nothing as we have no new information to incorporate into weights
  if (observations.size()==0)
    return;

  // loop over all particles
  Particle p;
  double sum_weights = 0.0; // sum of calculated weights for normalisation
  for (int i=0; i<num_particles; i++) {
    p = particles[i];

    // translate observations from particle coordinates to map coordinates
    std::vector<LandmarkObs> observations_map_coord; // observations in map coordinates
    LandmarkObs obs; // temp object
    for (int o=0; o<observations.size(); o++)
    {
      obs = observations[o]; // copy observation as new object
      CarToMap(obs.x, obs.y, p.x, p.y, p.theta, obs.x, obs.y); // we are overwriting obs.x and obs.y inside function
      observations_map_coord.push_back(obs); // copies obs as new object
    }

    // Associate observations with landmarks using 'nearest neighbour'
    // The alternative would be to use Hungarian assignment algorithm (http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html)
    // TODO: I may implement Hungarian algorithm later if I have time.
    // The 'nearest neighbour' overview is:
    // 1. find all map landmarks withing sensor_range around particle, save in close_landmarks
    // 2. loop over landmarks and find closest observations to each.
    std::vector<LandmarkObs> close_landmarks;
    std::vector<double> close_landmark_distances;
    LandmarkObs particle_position; // new object to convert particle position format to observation format
    particle_position.x = p.x;
    particle_position.y = p.y;
    for (int l=0; l<map_landmarks.landmark_list.size(); l++)
    {
      obs.x = map_landmarks.landmark_list[l].x_f; // reusing temp object obs
      obs.y = map_landmarks.landmark_list[l].y_f;
      double d = ObservationsDistance(obs, particle_position);
      if (d <= sensor_range) {
        close_landmarks.push_back(obs); // push_back copies obs as new object
        close_landmark_distances.push_back(d);
      }
    }
    // now sort the landmarks by distance from the particle
    // using approach: http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
    vector<int> sort_idx(close_landmark_distances.size());
    std::iota(sort_idx.begin(), sort_idx.end(), 0); // original indices 0, 1, 2...
    // sort using lambda function
    std::sort(sort_idx.begin(), sort_idx.end(),
              [&close_landmark_distances](int i1, int i2){return close_landmark_distances[i1]<close_landmark_distances[i2];});
    // do the nearest neighbour and calculate probabilities at the same time
    // cummulative probability for this particle for all observations
    double prob = 1.0;
    for (int l=0; l<close_landmarks.size(); l++)
    {
      LandmarkObs lm = close_landmarks[sort_idx[l]/*using sort by distance*/];
      // note: observations are removed from observations list (own copy in this function)
      double dist, min_dist = numeric_limits<double>::max();
      int min_o = -1;
      for (int o=0; o<observations_map_coord.size(); o++)
      {
        dist = ObservationsDistance(lm, observations_map_coord[o]);
        if (dist < min_dist)
        {
          min_dist = dist;
          min_o = o;
        }
      }
      if (min_o>-1)
      {
        // we matched landmark to observation. calculate probability and remove that observation from further consideration
        prob *= Gaussian2DNoCorrelation(observations_map_coord[min_o].x, lm.x,
                                        observations_map_coord[min_o].y, lm.y,
                                        std_landmark[0], std_landmark[1]);
        observations_map_coord.erase(observations_map_coord.begin()+min_o); // this is horrible in terms of performance as it has to shift remaining elements
      } else
      {
        // this is the case we have run out of observations to match.
        // one way here is to assign minimum probability to this landmarks. but that's not good. it will produce zero total probability
        // instead we just ignore this landmark. may be it was occluded or there is fog and sensor cannot see it or whatever.
        // let's just focus on landmarks we can match
      }
    } // loop over close landmarks

    // Use cummulative probability as the particle weight.
    // But if there are no landmarks close to this particle (within range) then weight should be zero.
    // So we assume the landmark density is such that there must be a landmark within our sensor range.
    // If it is not the case due to poor map the car should rely on other ways to localize. visual odometry, gps etc.
    // If map is detailed then our particle is really off the map and should be eliminated, hence zero weight.
    if (close_landmarks.size()) {
      p.weight = prob;
      sum_weights += p.weight;
    } else
      p.weight = 0.;

    particles[i] = p;
  } // loop over particles

  if (sum_weights < 1e-9)
    cout << "VERY SMALL WEIGHTS" << endl;

  if (sum_weights == 0)
    cout << "ZERO WEIGHTS SUM. PROBABLY KAPUT!" << endl;

  //normalize weights
  double max_weight = 0.;
  for (int i=0; i<num_particles; i++) {
    particles[i].weight /= (sum_weights > 1e-14 ? sum_weights : 1e-14);
    if (particles[i].weight > max_weight)
      max_weight = particles[i].weight;
    weights[i] = particles[i].weight;
  }

  cout << "max weight: " << max_weight << endl;
}

void ParticleFilter::resample() {
  // Resample particles with replacement with probability proportional to their weight.
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::discrete_distribution<> d(weights.begin(), weights.end()); // initialize discrete distribution over particles with weights that we calculated
  // sample from the distribution
  std::vector<Particle> particles_new(num_particles);
  std::map<int, int> m; // this is to count how many particles survived
  for (int i=0; i<num_particles; i++)
  {
    int sample = d(gen);
    ++m[sample];
    particles_new[i] = particles[sample];
  }
  cout << "number of particles surviving: " << m.size() << endl;

  particles = particles_new;
}

// write particles to the specified file
void ParticleFilter::write(std::string filename) {
  std::ofstream dataFile;
  dataFile.open(filename, std::ios::app);
  for (int i = 0; i < num_particles; ++i) {
    dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << " " << particles[i].weight <<"\n";
  }
  dataFile.close();
}
