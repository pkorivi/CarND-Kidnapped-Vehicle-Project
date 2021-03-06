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
	num_particles = 20;
	default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	for(int i =0;i<num_particles;i++){
		Particle part;
		part.id = i;
		part.x = dist_x(gen);
		part.y = dist_y(gen);
		part.theta = dist_theta(gen);
		part.weight = 1.0;
		particles.push_back(part);
		weights.push_back(1.0);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	for(auto &particle : particles){
		//No chnage in yaw, thus dont consider yaw rate while calculating deltax
		if(fabs(yaw_rate) < 0.0001){
			particle.x += velocity*delta_t*cos(particle.theta);
			particle.y += velocity*delta_t*sin(particle.theta); 	
		}
		else{
			particle.x += (velocity/yaw_rate) * (sin(particle.theta + yaw_rate*delta_t) - sin(particle.theta));
			particle.y += (velocity/yaw_rate) * (-cos(particle.theta + yaw_rate*delta_t) + cos(particle.theta));
			particle.theta += yaw_rate*delta_t;
		}
		//Random Noise addition
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
	//Simple nearest neighbour implementation
	for(auto &obs : observations){
		double min_dist = numeric_limits<double>::max();
		int map_id = -1;
		for(auto pred : predicted){
			double dis = dist(obs.x, obs.y, pred.x, pred.y);
			if( dis< min_dist ){
				min_dist = dis;
                map_id = pred.id;
			}
		}
		obs.id = map_id;
	}
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
		//Transform observations to map coordinates 
		std::vector<LandmarkObs> transformed_obs;
		for(auto obs : observations){
			double xmap = particle.x + cos(particle.theta)*obs.x - sin(particle.theta)*obs.y;
			double ymap = particle.y + sin(particle.theta)*obs.x + cos(particle.theta)*obs.y;
			transformed_obs.push_back(LandmarkObs{obs.id,xmap,ymap});
		}

		//Set associations between sensor observations and map lanmarks
		dataAssociation(land_mark_obs,transformed_obs);

		//reset weight
		particle.weight = 1.0;
		for(auto obs : transformed_obs){
			double as_x=0,as_y=0;
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
			double norm_fact = 1/(2*M_PI*std_x*std_y);
			double exponent = (pow(as_x - obs.x, 2) / (2 * pow(std_x, 2)) + (pow(as_y - obs.y, 2) / (2 * pow(std_y, 2))));
			double obs_wt = norm_fact*exp(-exponent);
			//check to see that the values dont go too low
			if(obs_wt < 0.000001){
				particle.weight *=0.000001;	
			}
			else{
				particle.weight *= obs_wt;
			}
		}
	}
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//Gather weights for samling wheel
	vector<double> weights;
	for(auto particle : particles){
		weights.push_back(particle.weight);
	}
	//maxweight of weights
	double max_wt = *max_element(weights.begin(),weights.end());

	default_random_engine gen;
	//random distribution between 0-max_Wt
	uniform_real_distribution<double> uni_real_dist(0.0, max_wt);
	//random start point for resampling wheel
	uniform_int_distribution<int>  uni_int_dist(0, num_particles - 1);
	auto rand_index = uni_int_dist(gen);
	

	double beta = 0.0;
	vector<Particle> new_particles;
	//spin---
	for(int i=0;i<num_particles;i++){
		beta += uni_real_dist(gen)*2;
		while(beta > weights[rand_index]){
			beta -= weights[rand_index];
			rand_index = (rand_index+1)%num_particles;
		}
		new_particles.push_back(particles[rand_index]);
	}
	particles = new_particles;
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
    return particle;
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
