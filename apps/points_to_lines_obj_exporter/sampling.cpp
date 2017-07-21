#include "sampling.h"

  

    void sampling::
    apply_distance_optimization_sampling(){

    }

    std::vector<point> sampling::
    apply_random_gridbased_sampling(std::vector<point> & input_cluster, std::mt19937 g){

        //sampling step (reduction of points to chosen % of original amout of points per cluster)
        float const sampling_rate = 0.7; //persentage points to remain after sampling
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

        std::vector<std::vector<point*>> vector_of_cells (num_cells_x_direction * num_cells_z_direction);

        //split all cluter points into respective grid cells
        for (auto & current_point : input_cluster){ //TODO check constness
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