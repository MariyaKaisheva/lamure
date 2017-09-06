#include <lamure/npr/binning.h>
#include <lamure/npr/utils.h>

#include <queue>
#include <list>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


namespace npr {
namespace binning {
bool 
  evaluate_similarity(bin const& bin_A, 
                      bin const& bin_B){

 	//coefficents used by computation of Jaccard index
 	uint m_11 = 0; //both binary values = 1
 	uint m_01 = 0; //first binary value = 0, second binary value = 1
 	uint m_10 = 0; //first binary value = 1, second binary value = 0
 	uint m_00 = 0; //both binary values = 0
 	uint num_cells = bin_A.binary_image_.size(); //bin_A and bin_B have binary_image vectors of equal size

 	for(uint cell_iterator = 0; cell_iterator < num_cells; ++cell_iterator){

 		if( (bin_A.binary_image_[cell_iterator] > 0) && (bin_B.binary_image_[cell_iterator] > 0) ){
 			++m_11;
 		}else if( (bin_A.binary_image_[cell_iterator] == 0) && (bin_B.binary_image_[cell_iterator] != 0) ){
 			++m_01;
 		}else if( (bin_A.binary_image_[cell_iterator] != 0) && (bin_B.binary_image_[cell_iterator] == 0) ){
 			++m_10;
 		}else if( (bin_A.binary_image_[cell_iterator] == 0) && ( bin_B.binary_image_[cell_iterator] == 0) ){
 			++m_00;
 		}
 	}
 
 	float jaccard_enumerator = m_11 ;
 	float jaccard_denominator = (m_01 + m_10 + m_11);

 	float jaccard_index = 0.0;

 	if( 0.0 != jaccard_denominator ) {
 		jaccard_index = jaccard_enumerator / jaccard_denominator;
 	}


    float const splitting_depth_dependent_sensitivity = 0.35;

 	//std::cout << "Jaccard Index: " << jaccard_index << std::endl;
  	if(jaccard_index > splitting_depth_dependent_sensitivity){
 		return true; //the 2 binary images are similar => second bin can be removed
 	}else{
 		return false;
 	}
}

//function returns TRUE if distance between 2 bins is beyond desired threshold
bool evaluate_if_distance_is_too_large(bin const& bin_A, 
                                       bin const& bin_B,
                                       float max_distance){

    auto pos_A = bin_A.pos_along_slicing_axis_;
    auto pos_B = bin_B.pos_along_slicing_axis_;

    float epsilon_factor = 1.001; // to mitigate for roundoff errors if max_distance = n * min_distance
    //std::cout << "!!! Distance between bins: " << std::fabs(pos_A - pos_B) << "; MAX: " << max_distance << "\n";  
    return std::fabs(pos_A - pos_B) * epsilon_factor >= max_distance; 
}

std::vector<bin> 
 generate_all_bins(std::vector<xyzall_surfel_t> const& all_surfels, 
                   float initial_bound_value, uint& max_num_layers,
                   float max_distance_between_two_neighbouring_bins,
                   bool verbose){
 	
 	std::vector<bin> bins;
 

        bool static_binnig = false; 

        if(static_binnig){
            auto num_elements = all_surfels.size();
            auto height = all_surfels.at(num_elements-1).pos_coordinates[1] - all_surfels.at(0).pos_coordinates[1]; //all_surfels vector is already sorted at this point 
            float const offset = height / max_num_layers;
            float current_y_min = all_surfels.at(0).pos_coordinates[1];
            float current_y_max = current_y_min + offset;

            for(uint i = 0; i < max_num_layers; ++i) {
                float current_y_mean = (current_y_min + current_y_max) / 2.0;
                /*if(initial_bound_value >= current_y_max - current_y_mean) {
                  throw  std::runtime_error("density thershold might be too low");
                } */

                bins.emplace_back(all_surfels, initial_bound_value, initial_bound_value, current_y_mean);
                current_y_min = current_y_max;
                current_y_max += offset;
            }
        }else{ //dynamic binnig 
                
            auto bounding_corners = utils::compute_bounding_corners(all_surfels);
                
            uint grid_resolution = 80; //num cells pro dim for generation of binary_image for bin comparison

                std::list<bin> working_list_of_bins;
                auto bound_value = (fabs(bounding_corners.max_y - bounding_corners.min_y) / max_num_layers) /2.0;
                float new_bin_location = bounding_corners.min_y + bound_value;
                for(uint bin_counter = 0; bin_counter < max_num_layers; ++bin_counter){                    
                    working_list_of_bins.emplace_back(all_surfels, bound_value, bound_value, new_bin_location);
                    working_list_of_bins.back().evaluate_content_to_binary(bounding_corners, grid_resolution);
                    new_bin_location += 2.0 * bound_value;
                }

                std::list<bin>::iterator it1,it2;
                it1 = it2 = working_list_of_bins.begin();

                ++it2;//std::advance(it2, 1);
                bins.push_back(*it1);
                while(it2 != --working_list_of_bins.end()){
                    bool bins_are_similar = evaluate_similarity(*it1, *it2/*, true*/);
                    bool distance_exceeds_max_distance_threshold = evaluate_if_distance_is_too_large(*it1, *it2, max_distance_between_two_neighbouring_bins); 
                    
                    //keep data if bins are too dissimilar, or if distance between then is NOT too large
                    if(bins_are_similar && (!distance_exceeds_max_distance_threshold)){

                        auto pos_A = it1->pos_along_slicing_axis_;
                        auto pos_B = it2->pos_along_slicing_axis_;

                        it2 = working_list_of_bins.erase(it2);

                    }
                    else{
                        ++it1;
                        ++it2;
                        bins.push_back(*it1);
                        
                    }
                }

                bins.push_back(*it2);
        }

    
    if(verbose){
        std::cout << "Starting with " << max_num_layers << " slicing layers, ";
        std::cout << "and finally returning " << bins.size() << " of them after binning\n";
    }
    
 	return bins;
}


} //namespace binning
} //namespace npr