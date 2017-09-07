// Copyright (c) 2014 Bauhaus-Universitaet Weimar
// This Software is distributed under the Modified BSD License, see license.txt.
//
// Virtual Reality and Visualization Research Group 
// Faculty of Media, Bauhaus-Universitaet Weimar
// http://www.uni-weimar.de/medien/vr

#version 420 core

INCLUDE common/attribute_dequantization_header.glsl

uniform mat4 inv_mv_matrix;
uniform float model_radius_scale;
uniform float point_size_factor;

out VertexData {
  vec3 pass_ms_u;
  vec3 pass_ms_v;
  vec3 pass_normal;
} VertexOut;

INCLUDE common/attribute_dequantization_functions.glsl
INCLUDE common/compute_tangent_vectors.glsl

void main() {
  // the "in" prefix is kept to be consistent with the interface of the uncompressed shaders.
  // conceptually, the decompressed attributes are the input for the next stages.
  vec3 in_position = vec3(0.0);
  vec3 in_normal   = vec3(0.0);
  float in_radius  = 0.0;

  dequantize_surfel_attributes_without_color(in_qz_pos_xy_16_16, in_qz_pos_z_normal_enum_16_16, in_rgb_777_and_radius_11, //compressed v-attributes
                                             in_position, in_normal, in_radius); //decompressed v-attributes



  // precalculate tangent vectors to establish the surfel shape
  vec3 tangent = vec3(0.0);
  vec3 bitangent = vec3(0.0);

  compute_tangent_vectors(in_normal, in_radius, tangent, bitangent);

  vec3 normal = normalize((inv_mv_matrix * vec4(in_normal, 0.0)).xyz);  

  // passed attributes: vertex shader -> geometry shader
  VertexOut.pass_ms_u = tangent;
  VertexOut.pass_ms_v = bitangent;
  VertexOut.pass_normal = normal;
  gl_Position = vec4(in_position, 1.0);
}