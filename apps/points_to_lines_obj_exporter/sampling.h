#ifndef SAMPLING_H
#define SAMPLING_H

#include "utils.h"

#include <vector>
//#include <random>

		namespace sampling {
			std::vector<point> apply_distance_optimization_sampling(std::vector<point> & input_cluster, uint num_remaining_points);
			std::vector<point> apply_random_gridbased_sampling(std::vector<point> & input_cluster, std::mt19937 g);
		}	

#endif //SAMPLING_H