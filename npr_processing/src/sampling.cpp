#include <lamure/npr/sampling.h>

#include <lamure/types.h>

#include <algorithm>
#include <limits>
  
namespace npr {
namespace sampling {
  std::vector<point> 
  apply_distance_optimization_sampling(std::vector<point> & input_cluster, uint num_remaining_points){

    uint input_size = input_cluster.size();

    if( num_remaining_points >= input_size ) {
      return input_cluster;
    }

    std::vector<point> sampled_cluster;
    //sort points to  make sure we beggin with an outer most location
    std::sort(input_cluster.begin(), input_cluster.end(), [](point const& A, point const& B) {return A.pos_coordinates_[0] < B.pos_coordinates_[0];});


    std::vector<bool> is_used_vec (input_size, false);

    //sample the starting point
    sampled_cluster.push_back(input_cluster[0]);
    is_used_vec[0] = true;


    for(uint counter = 1; counter < num_remaining_points; ++counter) {
      float max_distance_to_candidate_point = -1.0;
      int candidate_point_index = -1;

      for(uint cluster_index = 0; cluster_index < input_size; ++cluster_index) {
        float min_distance_to_sampled_point = std::numeric_limits<float>::max();
        auto const& current_point = input_cluster[cluster_index];

        if(!is_used_vec[cluster_index]){
            auto start_point_coord = lamure::vec3f(current_point.pos_coordinates_[0], current_point.pos_coordinates_[1], current_point.pos_coordinates_[2]);

            // inner loop: get minimum distance of current point to "border point" of already selected set
            for(auto const& sampled_point : sampled_cluster){
              auto end_point_coord = lamure::vec3f(sampled_point.pos_coordinates_[0], sampled_point.pos_coordinates_[1], sampled_point.pos_coordinates_[2]);
              float current_distance = utils::compute_distance(start_point_coord, end_point_coord);

              min_distance_to_sampled_point = std::min(min_distance_to_sampled_point, current_distance);

            }

            // check if the minimum distance was increased. in case it was: the candidate is better.x_index
            if(max_distance_to_candidate_point < min_distance_to_sampled_point){
              max_distance_to_candidate_point = min_distance_to_sampled_point;
              candidate_point_index = cluster_index;
            }

        }
      }

      //candidate_point_index =  std::max( int64_t(candidate_point_index), int64_t(0) );
      auto const& selected_point = input_cluster[candidate_point_index];
      is_used_vec[candidate_point_index] = true;
      sampled_cluster.push_back(selected_point);

    }

    return sampled_cluster;
  }

  
  std::vector<point> 
  apply_random_gridbased_sampling(std::vector<point> & input_cluster, std::mt19937 g){

      //sampling step (reduction of points to chosen % of original amout of points per cluster)
      float const sampling_rate = 1.0; //persentage points to remain after sampling
      unsigned const min_number_points_in_cell = 2;

      float min_x = std::numeric_limits<float>::max();
      float max_x = std::numeric_limits<float>::lowest();
      float min_z = std::numeric_limits<float>::max();
      float max_z = std::numeric_limits<float>::lowest();
      for (auto & current_point : input_cluster){
        min_x = std::min(min_x, current_point.pos_coordinates_[0]);
        max_x = std::max(max_x, current_point.pos_coordinates_[0]);
        min_z = std::min(min_z, current_point.pos_coordinates_[2]);
        max_z = std::max(max_z, current_point.pos_coordinates_[2]);
      }             

      auto num_cells_x_direction = 10;
      auto num_cells_z_direction = 10;
      auto cell_length = (max_z - min_z) / num_cells_z_direction;
      auto cell_width = (max_x - min_x) / num_cells_x_direction;

      std::vector<std::vector<point const*>> vector_of_cells (num_cells_x_direction * num_cells_z_direction);

      //split all cluter points into respective grid cells
      for (auto const& current_point : input_cluster){ //TODO check constness
        auto x_index = std::min(num_cells_x_direction - 1, std::max(0, int( (current_point.pos_coordinates_[0] - min_x) / cell_width)));
        auto z_index = std::min(num_cells_z_direction - 1, std::max(0, int( (current_point.pos_coordinates_[2]- min_z) / cell_length)));
        int64_t cell_index = z_index * num_cells_x_direction + x_index;
        auto& current_cell = vector_of_cells[cell_index];
        current_cell.push_back(&current_point);
      }
     // auto num_points_in_cluster = current_cluster.size();
      //int total_num_remaining_points = num_points_in_cluster * sampling_rate;
      std::vector<point> sampled_cluster;
      for (auto& current_cell : vector_of_cells) {

        std::shuffle(current_cell.begin(), current_cell.end(), g);
        uint64_t num_point_to_remain_in_cell = current_cell.size()*sampling_rate;

        if( min_number_points_in_cell < current_cell.size() ) {
          current_cell.resize(num_point_to_remain_in_cell);
          current_cell.shrink_to_fit(); 
        } else {
          //std::cout << "Did not reduce num points of this cell. " << "\n";
        }
        for(auto& current_point : current_cell){
          sampled_cluster.push_back(*current_point);
        }
      }
      return sampled_cluster;
    }
} //namespace sampling
} //namespace npr