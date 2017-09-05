#ifndef COLOR_HASH_MAP_H
#define COLOR_HASH_MAP_H

#include <lamure/types.h>

namespace npr {

extern lamure::vec3b color_array[1022];

uint32_t id_to_color_hash( uint32_t x );

} //namespace npr

#endif //COLOR_HASH_MAP_H