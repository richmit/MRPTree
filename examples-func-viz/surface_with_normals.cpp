// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      surface_with_normals.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Simple surface function graph.@EOL
 @std       C++23
 @copyright 
  @parblock
  Copyright (c) 2024, Mitchell Jay Richling <http://www.mitchr.me/> All rights reserved.

  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this list of conditions, and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions, and the following disclaimer in the documentation
     and/or other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software
     without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
  OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGE.
  @endparblock
 @filedetails   

  Surface normals may be used by many visualization tools to render smoother results.  In this example we demonstrate:

   - How to compute a surface gradient for a function plot
   - How to unitize the gradient into a surface normal
   - How to add the normal to the sample data stored by a `MRPTree`
   - How to include normals in the cell complex
   - How to increase sampling with a SDF function
   - How to increase sampling near humps by testing derivatives
   - How to balance a tree
   - How to dump a cell complex into various file types
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d5rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t damp_cos_wave(tt_t::drpt_t xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double d = x*x+y*y;
  double m = std::exp(-d/4);
  double s = std::sqrt(d);
  double z = m*cos(4*s);
  double dx = -(cos((4 * s)) * s + 4 * sin( (4 * s))) * x * exp(-x * x / 2 - y * y / 2);
  double dy = -(cos((4 * s)) * s + 4 * sin( (4 * s))) * y * exp(-x * x / 2 - y * y / 2);
  double dd =   -m*(cos(4*s)*s+8*sin(4*s));
  if (s>1.0e-5) {
    dx = dx / s;
    dy = dy / s;
    dd = dd / (4 * s);
  } else {
    dx = 1;
    dy = 1;
    dd = 1;
  }
  double nm = std::sqrt(1+dx*dx+dy*dy);
  return {z, -dx/nm, -dy/nm, 1/nm, dd};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double circle_sdf(double r, tt_t::drpt_t xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double m = x*x+y*y;
  return (r*r-m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-2.1, -2.1}, 
            { 2.1,  2.1});
  cc_t ccplx;
  
  // Make a few samples on a uniform grid
  tree.refine_grid(2, damp_cos_wave);

  // The humps need extra samples.  We know where they are, and we could sample on them with an SDF like this:
  // for(double i: {0, 1, 2, 3}) {
  //   double r = i*std::numbers::pi/4;
  //   tree.refine_leaves_recursive_cell_pred(6, damp_cos_wave, [&tree, r](int i) { return (tree.cell_cross_sdf(i, std::bind_front(circle_sdf, r))); });
  // }

  // Alternately, we can test the derivative values to identify the humps
  // tree.refine_leaves_recursive_cell_pred(6, damp_cos_wave, [&tree](tt_t::diti_t i) { return tree.cell_cross_range_level(i, 1, 0.0); });
  // tree.refine_leaves_recursive_cell_pred(6, damp_cos_wave, [&tree](tt_t::diti_t i) { return tree.cell_cross_range_level(i, 2, 0.0); });

  // Lastly we can use the directional derivative radiating from the origin
  tree.refine_leaves_recursive_cell_pred(6, damp_cos_wave, [&tree](tt_t::diti_t i) { return tree.cell_cross_range_level(i, 4, 0.0); });

  // Balance the three to the traditional level of 1 (no  cell borders a cell more than half it's size)
  tree.balance_tree(1, damp_cos_wave);

  tree.dump_tree(5);

  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FDOMAIN, 0}, 
                                 {tc_t::val_src_spc_t::FDOMAIN, 1},
                                 {tc_t::val_src_spc_t::FRANGE,  0}});

  // Note we use the single argument version of create_named_datasets() because we don't want to name elements 3, 4, & 5 (the components of the normal
  // Note if we had placed the ddiv component right after z, then we could have used the two argument version...
  ccplx.set_data_name_to_data_idx_lst({{"x",        {0}},
                                       {"y",        {1}},
                                       {"z=f(x,y)", {2}},
                                       {"ddiv",     {6}},
                                       {"NORMALS",  {3,4,5}}});

  ccplx.dump_cplx(5);

  ccplx.write_legacy_vtk("surface_with_normals.vtk", "surface_with_normals");
  ccplx.write_xml_vtk(   "surface_with_normals.vtu", "surface_with_normals");
  ccplx.write_ply(       "surface_with_normals.ply", "surface_with_normals");
}
/** @endcond */
