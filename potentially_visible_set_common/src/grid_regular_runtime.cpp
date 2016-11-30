#include "lamure/pvs/grid_regular_runtime.h"

#include <fstream>
#include <stdexcept>
#include <string>
#include <climits>

namespace lamure
{
namespace pvs
{

grid_regular_runtime::
grid_regular_runtime() : grid_regular_runtime(1, 1.0, scm::math::vec3d(0.0, 0.0, 0.0))
{
}

grid_regular_runtime::
grid_regular_runtime(const size_t& number_cells, const double& cell_size, const scm::math::vec3d& position_center)
{
	file_path_pvs_ = "";
	create_grid(number_cells, cell_size, position_center);
}

grid_regular_runtime::
~grid_regular_runtime()
{
	cells_.clear();
}

size_t grid_regular_runtime::
get_cell_count() const
{
	return cells_.size();
}

view_cell* grid_regular_runtime::
get_cell_at_index(const size_t& index)
{
	view_cell* current_cell = &cells_.at(index);
	if(current_cell->contains_visibility_data() || file_path_pvs_ == "")
	{
		return current_cell;
	}

	// One line per model.get_cell_at_index
	node_t single_model_bytes = 0;

	for(model_t model_index = 0; model_index < ids_.size(); ++model_index)
	{
		node_t num_nodes = ids_.at(model_index);

		// If the number of node IDs is not dividable by 8 there is one additional character.
		node_t addition = 0;
		if(num_nodes % CHAR_BIT != 0)
		{
			addition = 1;
		}

		node_t line_size = (num_nodes / CHAR_BIT) + addition;
		single_model_bytes += line_size;
	}

	std::fstream file_in;
	file_in.open(file_path_pvs_, std::ios::in | std::ios::binary);

	if(!file_in.is_open())
	{
		throw std::invalid_argument("invalid file path: " + file_path_pvs_);
	}

	file_in.seekg(index * single_model_bytes);

	for(model_t model_index = 0; model_index < ids_.size(); ++model_index)
	{
		node_t num_nodes = ids_.at(model_index);

		// If the number of node IDs is not dividable by 8 there is one additional character.
		node_t addition = 0;
		if(num_nodes % CHAR_BIT != 0)
		{
			addition = 1;
		}

		node_t line_size = (num_nodes / CHAR_BIT) + addition;

		for(node_t character_index = 0; character_index < line_size; ++character_index)
		{
			char current_byte = 0x00;
			file_in.read(&current_byte, 1);
			
			for(unsigned short bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
			{
				bool visible = ((current_byte >> bit_index) & 1) == 0x01;
				current_cell->set_visibility(model_index, (character_index * CHAR_BIT) + bit_index, visible);
			}
		}
	}

	file_in.close();

	return current_cell;
}

view_cell* grid_regular_runtime::
get_cell_at_position(const scm::math::vec3d& position)
{
	size_t num_cells = std::pow(cells_.size(), 1.0f/3.0f);
	double half_size = cell_size_ * (double)num_cells * 0.5f;
	scm::math::vec3d distance = position - (position_center_ - half_size);

	size_t index_x = (size_t)(distance.x / cell_size_);
	size_t index_y = (size_t)(distance.y / cell_size_);
	size_t index_z = (size_t)(distance.z / cell_size_);

	// Check calculated index so we know if the position is inside the grid at all.
	if(index_x < 0 || index_x >= num_cells ||
		index_y < 0 || index_y >= num_cells ||
		index_z < 0 || index_z >= num_cells)
	{
		return nullptr;
	}

	size_t general_index = (num_cells * num_cells * index_z) + (num_cells * index_y) + index_x;
	return get_cell_at_index(general_index);
}

void grid_regular_runtime::
save_grid_to_file(const std::string& file_path) const
{
	std::fstream file_out;
	file_out.open(file_path, std::ios::out);

	if(!file_out.is_open())
	{
		throw std::invalid_argument("invalid file path: " + file_path);
	}

	// Grid file type.
	file_out << "regular" << std::endl;

	// Number of grid cells per dimension.
	file_out << std::pow(cells_.size(), 1.0f/3.0f) << std::endl;

	// Grid size and position
	file_out << cell_size_ << std::endl;
	file_out << position_center_.x << " " << position_center_.y << " " << position_center_.z << std::endl;

	file_out.close();
}

void grid_regular_runtime::
save_visibility_to_file(const std::string& file_path, const std::vector<node_t>& ids) const
{
	std::fstream file_out;
	file_out.open(file_path, std::ios::out | std::ios::binary);

	if(!file_out.is_open())
	{
		throw std::invalid_argument("invalid file path: " + file_path);
	}

	// Iterate over view cells.
	for(size_t cell_index = 0; cell_index < cells_.size(); ++cell_index)
	{
		// Iterate over models in the scene.
		for(lamure::model_t model_id = 0; model_id < ids.size(); ++model_id)
		{
			node_t num_nodes = ids.at(model_id);
			char current_byte = 0x00;

			// Iterate over nodes in the model.
			for(lamure::node_t node_id = 0; node_id < num_nodes; ++node_id)
			{
				if(cells_.at(cell_index).get_visibility(model_id, node_id))
				{
					current_byte |= 1 << (node_id % CHAR_BIT);
				}

				// Flush character if either 8 bits are written or if the node id is the last one.
				if((node_id + 1) % CHAR_BIT == 0 || node_id == (num_nodes - 1))
				{
					file_out.write(&current_byte, 1);
					current_byte = 0x00;
				}
			}
		}
	}

	file_out.close();
}

bool grid_regular_runtime::
load_grid_from_file(const std::string& file_path)
{
	std::fstream file_in;
	file_in.open(file_path, std::ios::in);

	if(!file_in.is_open())
	{
		return false;
	}

	// Start reading the header info which is used to recreate the grid.
	std::string grid_type;
	file_in >> grid_type;
	if(grid_type != "regular")
	{
		return false;
	}

	size_t num_cells;
	file_in >> num_cells;

	double cell_size;
	file_in >> cell_size;

	double pos_x, pos_y, pos_z;
	file_in >> pos_x >> pos_y >> pos_z;

	create_grid(num_cells, cell_size, scm::math::vec3d(pos_x, pos_y, pos_z));

	file_in.close();
	return true;
}


bool grid_regular_runtime::
load_visibility_from_file(const std::string& file_path, const std::vector<node_t>& ids)
{
	std::fstream file_in;
	file_in.open(file_path, std::ios::in | std::ios::binary);

	if(!file_in.is_open())
	{
		return false;
	}

	file_in.close();

	file_path_pvs_ = file_path;
	ids_ = ids;

	return true;
}

void grid_regular_runtime::
create_grid(const size_t& num_cells, const double& cell_size, const scm::math::vec3d& position_center)
{
	cells_.clear();

	double half_size = (cell_size * (double)num_cells) * 0.5;		// position of grid is at grid center, so cells have a offset
	double cell_offset = cell_size * 0.5f;							// position of cell is at cell center

	for(size_t index_z = 0; index_z < num_cells; ++index_z)
	{
		for(size_t index_y = 0; index_y < num_cells; ++index_y)
		{
			for(size_t index_x = 0; index_x < num_cells; ++index_x)
			{
				scm::math::vec3d pos = position_center + (scm::math::vec3d(index_x , index_y, index_z) * cell_size) - half_size + cell_offset;
				cells_.push_back(view_cell_regular(cell_size, pos));
			}
		}
	}

	cell_size_ = cell_size;
	position_center_ = scm::math::vec3d(position_center);
}

}
}
