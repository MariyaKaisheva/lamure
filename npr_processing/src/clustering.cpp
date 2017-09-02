#include <lamure/npr/clustering.h>
#include <lamure/npr/utils.h>

namespace npr {
namespace clustering {

uint32_t const degree = 3;

std::pair<float, point>
find_nearest_neighbour (point const& start_point, std::vector<point> const& all_points) {
  float current_min = std::numeric_limits<float>::max();
  point current_nearest_neighbour;

  bool is_valid = false;

  for (auto& p : all_points){
    if((!p.is_used_)){
      float distance = utils::compute_distance(lamure::vec3f(p.pos_coordinates_[0], p.pos_coordinates_[1], p.pos_coordinates_[2]), 
                                        lamure::vec3f(start_point.pos_coordinates_[0], start_point.pos_coordinates_[1], start_point.pos_coordinates_[2]));
      if( distance > 0.0 && distance <= current_min ){
        current_min = distance;
        current_nearest_neighbour = p;

        is_valid = true;
      }
    }
  }

  if( !is_valid ) {
    current_min = std::numeric_limits<float>::max();
  }

  std::pair<float, point> result(current_min, current_nearest_neighbour);
  return result;
}


std::vector<point*>
find_near_neighbours(float max_distance, point const& start_point, std::vector<point>& all_points){
  std::vector<point*> neighbourhood;
  for(auto& p : all_points){

    float distance = utils::compute_distance(lamure::vec3f(p.pos_coordinates_[0], p.pos_coordinates_[1], p.pos_coordinates_[2]), 
                                      lamure::vec3f(start_point.pos_coordinates_[0], start_point.pos_coordinates_[1], start_point.pos_coordinates_[2]));
    if(distance <= max_distance ){
      neighbourhood.push_back(&p);
    }
  }
  
  return neighbourhood; 
}

void
expand_cluster(point const& current_point, std::vector<point*>& neighbourhood, std::vector<point>& all_points_per_layer, clusters_t& current_cluster, float eps, uint8_t min_points) {
  current_cluster.push_back(current_point);

  auto size = neighbourhood.size();

  for(uint32_t point_itr = 0; point_itr < size; ++point_itr) {
    auto& neighbour_point = neighbourhood.at(point_itr);  
    if(!neighbour_point->is_used_){ //if P' is not visited
      neighbour_point->is_used_ = true; //mark P' as visited
      auto next_neighbourhood = find_near_neighbours(eps, *neighbour_point, all_points_per_layer); //region Query (P', eps)
      if(next_neighbourhood.size() >= min_points){ 
        for(auto& next_neighbour_point : next_neighbourhood){
          if(!next_neighbour_point->is_used_) {
            if( std::find(neighbourhood.begin(), neighbourhood.end(), next_neighbour_point) == neighbourhood.end() ){
              neighbourhood.push_back(next_neighbour_point);
            }
          }
        }

      }
      size = neighbourhood.size();

      if( !neighbour_point->is_member_of_cluster_ ) {
        neighbour_point->is_member_of_cluster_ = true;
        current_cluster.push_back(*neighbour_point);
      }
    }
  }
}

std::vector<clusters_t>
create_DBSCAN_clusters (bins_t const& all_surfels_per_layer, float eps, uint8_t min_points){
  //create vector of points, containg all surfels per layer in point format
  std::vector<point> all_points_per_layer;
  uint32_t id_counter = 0; 
  for (auto const& s: all_surfels_per_layer){
    point current_point( const_cast<float*>(s.pos_coordinates), id_counter, s.r_, s.g_, s.b_);
    all_points_per_layer.push_back(current_point);
    ++id_counter;
  }

  std::vector<clusters_t> clusters_vec;


  //uint32_t current_cluster_id = 0;

  //iterate over all points 
  for (auto& current_point : all_points_per_layer) {
    if(current_point.is_used_){
      continue; //point has been already evaluated
    } 
    else{
      current_point.is_used_ = true;
      auto neighbourhood = find_near_neighbours(eps, current_point, all_points_per_layer);

      neighbourhood.reserve(all_surfels_per_layer.size());

      if(neighbourhood.size() < min_points){
        continue; //point is considered noise and is not added to any cluster
      }
      else{
        clusters_t current_cluster;
        expand_cluster(current_point, neighbourhood, all_points_per_layer, current_cluster, eps, min_points);
        if(current_cluster.size() > degree + 1){
          clusters_vec.push_back(current_cluster);
        }/*else{
          std::cout << "!!!Skipping cluster with too few points \n";
        }*/
      }
    }
  }

  return clusters_vec;
}


std::vector<clusters_t>
create_clusters (bins_t const& all_surfels_per_layer){
   
  //create vector of points, containg all surfels per layer in point format
  std::vector<point> all_points_per_layer;
  int32_t id_counter = 0; 
  for (auto const& s: all_surfels_per_layer){
    point current_point( const_cast<float*>(s.pos_coordinates), id_counter, s.r_, s.g_, s.b_, false);
    all_points_per_layer.push_back(current_point);
    ++id_counter;
  }

  //auto centroid = compute_average_position_per_layer(all_surfels_per_layer);
  std::vector<point> point_cluster;
  std::vector<clusters_t> clusters_vector;
  float distance_threshold = 60.0; 
  
  float last_measured_distance = 0;  //TODO think about this value

  //int cluster_counter  = 0; 

  //distribute points to clusters based on distance to closest available neigbour point
  for(uint point_iterator = 0; point_iterator < all_points_per_layer.size(); ++point_iterator){
      int point_index = point_iterator;
     // bool next_point_is_used = false;
      do {
        point* current_point = &all_points_per_layer[point_index];
        if(!current_point->is_used_){
            point_cluster.push_back(*current_point);
            current_point->is_used_ = true;

            //reset distance record for each new cluster; no distance has been stored so far
            if(point_cluster.size() < 2){
              last_measured_distance = std::numeric_limits<float>::max() / (2.0 * distance_threshold); 
            }
            
            auto search_result =  find_nearest_neighbour(*current_point, all_points_per_layer);
            auto distance_to_nearest_unused_point = search_result.first; 
            auto next_point = search_result.second; //function find_nearest_neighbour() makes sure next point is not used so far
            
            
            if(distance_to_nearest_unused_point <= (last_measured_distance * distance_threshold)){  
              point_cluster.push_back(next_point);
             // next_point.is_used_ = true;
              last_measured_distance  = distance_to_nearest_unused_point; 
              point_index = next_point.id_;
            }
            else{
             // next_point.is_used_ = false;
              //push current cluster to common container and start a new
              clusters_vector.push_back(point_cluster);
             // ++cluster_counter; 
              point_cluster.clear();
            }
            //next_point_is_used = next_point.is_used_;
        } else {
          break;
        }
      } while(true); //!next_point_is_used is always true 
  }  

  //std::cout << cluster_counter << std::endl;
  return clusters_vector; 
}

} //namespace clustering
} //namespace npr