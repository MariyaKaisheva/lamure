#include "binning.h"

#include "utils.h"

#include <queue>

bool binning::
  evaluate_similarity(bin const& bin_A, bin const& bin_B){

 	std::cout << "Starting to eval similarity\n";

 	//coefficents used by computation of Jaccard index
 	uint m_11 = 0; //both binary values = 1
 	uint m_01 = 0; //first binary value = 0, second binary value = 1
 	uint m_10 = 0; //first binary value = 1, second binary value = 0
 	uint m_00 = 0; //both binary values = 0
 	uint num_cells = bin_A.binary_image_.size(); //bin_A and bin_B have binary_image vectors of equal size

 	std::cout << "Size of binary image a: "  << num_cells << "\n";
 	std::cout << "Size of binary image b: "  << bin_B.binary_image_.size() << "\n";

 	for(uint cell_iterator = 0; cell_iterator < num_cells; ++cell_iterator){

 		uint32_t total_num_points_in_cells = bin_A.binary_image_[cell_iterator] + bin_B.binary_image_[cell_iterator];

 		if( (bin_A.binary_image_[cell_iterator] > 0) && (bin_B.binary_image_[cell_iterator] > 0) ){
 			++m_11;
 			//m_11 += total_num_points_in_cells;
 		}else if( (bin_A.binary_image_[cell_iterator] == 0) && (bin_B.binary_image_[cell_iterator] != 0) ){
 			++m_01;
 			//m_01 += total_num_points_in_cells;
 		}else if( (bin_A.binary_image_[cell_iterator] != 0) && (bin_B.binary_image_[cell_iterator] == 0) ){
 			++m_10;
 			//m_10 += total_num_points_in_cells;
 		}else if( (bin_A.binary_image_[cell_iterator] == 0) && ( bin_B.binary_image_[cell_iterator] == 0) ){
 			++m_00;
 			//m_00 += total_num_points_in_cells;
 		}
 	}
 	std::cout << "DENOMINATOR: " << m_01 + m_10 + m_11 << "\n";

 	float jaccard_enumerator = m_11 ;
 	float jaccard_denominator = (m_01 + m_10 + m_11);

 	float jaccard_index = 0.0;

 	if( 0.0 != jaccard_denominator ) {
 		jaccard_index = jaccard_enumerator / jaccard_denominator;
 	}

 	float base_sensitivity = 0.38;

 	float splitting_depth_dependent_sensitivity = base_sensitivity;

 	uint32_t max_split_depth = std::max(bin_A.bin_depth, bin_B.bin_depth);
 	
 	for(uint split_depth_idx = 0; split_depth_idx < max_split_depth; ++split_depth_idx) {
 		splitting_depth_dependent_sensitivity *= 0.58;
 	}

 	std::cout << "Jaccard Index: " << jaccard_index << "\n";
  	if(jaccard_index > splitting_depth_dependent_sensitivity){
 		return true; //the 2 binary images are similar => no additional bin should be created between them
 	}else{
 		return false;
 	}
}

std::vector<bin> binning::
 generate_all_bins(std::vector<xyzall_surfel_t> const& all_surfels, float const initial_bound_value){
 	std::vector<bin> bins;
 	auto bounding_corners = utils::compute_bounding_corners(all_surfels);
 	float current_pos_along_slicing_axis = bounding_corners.min_y + initial_bound_value;


 	size_t current_slice_count = 6;

 	float initial_offset = fabs(bounding_corners.max_y - bounding_corners.min_y - 2 * initial_bound_value) / (current_slice_count-1);
 	uint grid_resolution = 80; //num cells pro dim for generation of binary_image for bin comparison

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



 	while( (!working_queue.empty()) && (current_slice_count < 500)  ){
 		auto current_job = working_queue.front();
 		working_queue.pop();

 		auto top_id = current_job.top_bin_id_;
 		auto bottom_id = current_job.bottom_bin_id_;

 		auto & top_bin = bins[top_id];
 		auto & bottom_bin = bins[bottom_id];

 		std::cout << "Evaluating top bin: " << top_id << " against bottom bin: " << bottom_id << "\n";

 		bool are_bins_similar = evaluate_similarity(top_bin, bottom_bin);
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
 			std::cout << current_slice_count << "\n";
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

 	std::cout << "Finally returning " << bins.size() << " bins\n";

 	return bins;
}