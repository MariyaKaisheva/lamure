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
  evaluate_similarity(bin const& bin_A, bin const& bin_B, bool merge_binning){

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


    float splitting_depth_dependent_sensitivity = 0.3;

    if(!merge_binning){
        float base_sensitivity = 0.38;
        splitting_depth_dependent_sensitivity = base_sensitivity;
        uint32_t max_split_depth = std::max(bin_A.bin_depth, bin_B.bin_depth);

        for(uint split_depth_idx = 0; split_depth_idx < max_split_depth; ++split_depth_idx) {
            splitting_depth_dependent_sensitivity *= 0.58;
        } 
    }    
    
 	//std::cout << "Jaccard Index: " << jaccard_index << std::endl;
  	if(jaccard_index > splitting_depth_dependent_sensitivity){
 		return true; //the 2 binary images are similar => no additional bin should be created between them
 	}else{
 		return false;
 	}
}

std::vector<bin> 
 generate_all_bins(std::vector<xyzall_surfel_t> const& all_surfels, float const initial_bound_value, uint& max_num_layers){
 	
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
                if(initial_bound_value >= current_y_max - current_y_mean) {
                  throw  std::runtime_error("density thershold might be too low");
                } 

                bins.emplace_back(all_surfels, initial_bound_value, initial_bound_value, current_y_mean);
                current_y_min = current_y_max;
                current_y_max += offset;
            }
        }else{ //dynamic binnig 
                
            auto bounding_corners = utils::compute_bounding_corners(all_surfels);
                
            uint grid_resolution = 80; //num cells pro dim for generation of binary_image for bin comparison

            #if 0 //split-based adaptive binning
            
                float current_pos_along_slicing_axis = bounding_corners.min_y + initial_bound_value;

                size_t current_slice_count = 6;

                float initial_offset = fabs(bounding_corners.max_y - bounding_corners.min_y - 2 * initial_bound_value) / (current_slice_count-1);


                std::queue<evaluation_job> working_queue;

                //clreate 3 starting layers and corresponding image planes
                for (uint bin_idx = 0; bin_idx < current_slice_count; ++bin_idx){
                    bins.emplace_back(all_surfels, initial_bound_value, initial_bound_value, current_pos_along_slicing_axis);
                    current_pos_along_slicing_axis += initial_offset;
                    bins[bin_idx].evaluate_content_to_binary(bounding_corners, grid_resolution); //creates binary image representation of given bin
                }


                for(uint32_t bin_idx = 0; bin_idx < current_slice_count - 1; ++bin_idx) {
                    working_queue.push(evaluation_job(bin_idx, bin_idx+1));
                }

                if(max_num_layers < current_slice_count){
                    max_num_layers = 250;
                }

                while( (!working_queue.empty()) && (current_slice_count < max_num_layers)  ){
                    auto current_job = working_queue.front();
                    working_queue.pop();

                    auto top_id = current_job.top_bin_id_;
                    auto bottom_id = current_job.bottom_bin_id_;

                    auto & top_bin = bins[top_id];
                    auto & bottom_bin = bins[bottom_id];

                    //std::cout << "Evaluating top bin: " << top_id << " against bottom bin: " << bottom_id << "\n";

                    bool are_bins_similar = evaluate_similarity(top_bin, bottom_bin, false);
                    if (!are_bins_similar){


                        //create && insert new bin
                        float new_bin_location = bottom_bin.pos_along_slicing_axis_ + (top_bin.pos_along_slicing_axis_ - bottom_bin.pos_along_slicing_axis_) / 2.0;
                        auto new_upper_bound = top_bin.lower_bound_size_ / 2.0;
                        auto new_lower_bound = bottom_bin.upper_bound_size_ / 2.0;
                        bins.emplace_back(all_surfels, new_upper_bound, new_lower_bound, new_bin_location);

                        size_t new_slice_idx = bins.size() - 1;
                        if(bins[new_slice_idx].content_.empty()) {
                            bins.pop_back();
                            continue;
                        }

                        //remove redundant surfel && update bound values
                        top_bin.shrink_to_half_lower_bound(); 
                        bottom_bin.shrink_to_half_upper_bound();


                        bins[new_slice_idx].evaluate_content_to_binary(bounding_corners, grid_resolution); //creates binary image representation of given bin
                        bins[new_slice_idx].bin_depth = std::max(top_bin.bin_depth, bottom_bin.bin_depth) + 1;
                        //create respective evaluation jobs for the new bin
                        auto new_job_1 = evaluation_job(top_id, new_slice_idx);
                        auto new_job_2 = evaluation_job(new_slice_idx, bottom_id);
                        working_queue.push(new_job_1);
                        working_queue.push(new_job_2);

                        ++current_slice_count;

                        //std::cout << current_slice_count << "\n";
                    }
                }

                //remove unneeded data
                for(auto & current_bin : bins) {
                    current_bin.clear_binary_image();

                    //project all surfels to central XZ plane of the bin
                    for (auto & current_surfel : current_bin.content_) {
                        current_surfel.pos_coordinates[1] = current_bin.pos_along_slicing_axis_;
                    }
                }
                
            #else//merge-based dynamic binnig 

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

                std::advance(it2, 1);
                bins.push_back(*it1);
                while(it2 != working_list_of_bins.end()){
                    bool are_bins_similar = evaluate_similarity(*it1, *it2, true);
                    if(are_bins_similar){
                        it2 = working_list_of_bins.erase(it2);

                    }
                    else{
                        ++it1;
                        ++it2;
                        bins.push_back(*it1);
                    }
                }


            #endif
        }

    max_num_layers =  bins.size();
    //std::cout << "Finally returning " << bins.size() << " bins\n";
 	return bins;
}


} //namespace binning
} //namespace npr