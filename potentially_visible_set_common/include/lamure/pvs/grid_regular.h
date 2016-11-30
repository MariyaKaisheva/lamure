#ifndef LAMURE_PVS_REGULAR_GRID_H
#define LAMURE_PVS_REGULAR_GRID_H

#include <vector>
#include <string>

#include "lamure/pvs/grid.h"
#include "lamure/pvs/view_cell.h"
#include "lamure/pvs/view_cell_regular.h"

#include <scm/core/math.h>

namespace lamure
{
namespace pvs
{

class grid_regular : public grid
{
public:
	grid_regular();
	grid_regular(const size_t& number_cells, const double& cell_size, const scm::math::vec3d& position_center);
	~grid_regular();

	virtual size_t get_cell_count() const;
	virtual view_cell* get_cell_at_index(const size_t& index);
	virtual view_cell* get_cell_at_position(const scm::math::vec3d& position);

	virtual void save_grid_to_file(const std::string& file_path) const;
	virtual void save_visibility_to_file(const std::string& file_path, const std::vector<node_t>& ids) const;

	virtual bool load_grid_from_file(const std::string& file_path);
	virtual bool load_visibility_from_file(const std::string& file_path, const std::vector<node_t>& ids);

protected:
	void create_grid(const size_t& num_cells, const double& cell_size, const scm::math::vec3d& position_center);

private:
	double cell_size_;
	scm::math::vec3d position_center_;

	std::vector<view_cell_regular> cells_;
};

}
}

#endif
