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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 100;
	is_initialized = true;
	weights = vector<double>(num_particles);
	fill(weights.begin(), weights.end(), 0);
	particles = vector<Particle>(num_particles);
	Particle part;
	part.x = x;
	part.y = y;
	part.theta = theta;
	fill(particles.begin(),particles.end(),part);
	// Create normal distribution of noise parameters
	default_random_engine gen;
	normal_distribution<double> dist_x(0, std[0]);
	normal_distribution<double> dist_y(0, std[1]);
	normal_distribution<double> dist_theta(0, std[2]);
	for(auto &particle : particles)
	{
		particle.x +=  dist_x(gen);
		particle.y +=  dist_y(gen);
		particle.theta +=  dist_theta(gen);
	}

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	for(auto &particle : particles){
		particle.x += (velocity/yaw_rate) * (sin(particle.theta + yaw_rate*delta_t) - sin(particle.theta));
		particle.y += (velocity/yaw_rate) * (-cos(particle.theta + yaw_rate*delta_t) + cos(particle.theta));
		particle.theta += yaw_rate*delta_t;
	}

	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	for(auto &particle : particles)
	{
		particle.x +=  dist_x(gen);
		particle.y +=  dist_y(gen);
		particle.theta +=  dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	
	for(auto &particle : particles)
	{
		//Bringing map land marks  into Landmarkobs coordinate format
		std::vector<LandmarkObs> land_mark_obs;
		for(auto lm : map_landmarks.landmark_list){
			//checking if the values are in range
			if((fabs(lm.x_f - particle.x) <= sensor_range) && (fabs(lm.y_f - particle.y) <= sensor_range)){
				land_mark_obs.push_back(LandmarkObs{lm.id_i,lm.x_f,lm.y_f});
			}
		}
		//
		std::vector<LandmarkObs> transformed_obs;
		for(auto obs : observations){
			double xmap = particle.x + cos(particle.theta)*obs.x - sin(particle.theta)*obs.y;
			double ymap = particle.y + sin(particle.theta)*obs.x + cos(particle.theta)*obs.y;
			transformed_obs.push_back(LandmarkObs{obs.id,xmap,ymap});
		}

		//Set associations between sensor observations and map lanmarks
		dataAssociation(land_mark_obs,transformed_obs);

		particle.weight = 1.0;
		double as_x,as_y;
		for(auto obs : transformed_obs){
			for(auto lm : land_mark_obs){
				if(obs.id == lm.id){
					as_x = lm.x;
					as_y = lm.y;
					break;
				}
			}
			double std_x = std_landmark[0];
			double std_y = std_landmark[1];
			//observations weight
			double obs_wt = (1/(2*M_PI*std_x*std_y))* exp(-(pow(as_x - obs.x, 2) / (2 * pow(std_x, 2)) + (pow(as_y - obs.y, 2) / (2 * pow(std_y, 2)))));
			particle.weight *= obs_wt;
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
