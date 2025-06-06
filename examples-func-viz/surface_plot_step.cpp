// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      surface_plot_step.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-16
 @brief     Surface with a step discontinuity along a curve.@EOL
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

  The function illustrated here is defined on the entire plan, but has a step discontinuity on the unit circle.  Except on the unit circle, the function's
  derivative is zero.  If we sample on a uniform grid, some of the resulting polygons will have vertexes both inside and outside the unit circle -- they cross
  over the discontinuity!  When drawn we get a continuous surface with a circular bump!  It should look like the x-y plane has a circle cut out that is
  hovering one unit above the plane.

    - How to drive up the sample rate near a particular SDF -- so that we get higher resolution where the surface meets the plane.
    - How to delete triangles that cross over a discontinuity -- via an SDF.
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d1rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

tt_t::rrpt_t hover_circle(tt_t::drpt_t xvec) {
  return (xvec[0] * xvec[0] + xvec[1] * xvec[1] < 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::src_t unit_circle_sdf(tt_t::drpt_t xvec) {
  double m = xvec[0] * xvec[0] + xvec[1] * xvec[1];
  return (1-m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-1.5, -1.5}, 
            { 1.5,  1.5});
  cc_t ccplx;

  /* Here is another way to get fine samples on the circle, but with a SDF this time. */
  tree.refine_grid(5, hover_circle);

  /* Increase sample resolution on the unit circle.  Here we do that with an SDF. */
  tree.refine_leaves_recursive_cell_pred(7, hover_circle, [&tree](int i) { return (tree.cell_cross_sdf(i, unit_circle_sdf)); });

  /* Balance the three to the traditional level of 1 (no cell borders a cell more than half it's size) */
  tree.balance_tree(1, hover_circle);

  /* Take a peek at the raw tree data */
  tree.dump_tree(10);

  /* Generate a cell complex from the tree samples */
  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FDOMAIN, 0}, 
                                 {tc_t::val_src_spc_t::FDOMAIN, 1},
                                 {tc_t::val_src_spc_t::FRANGE,  0}});

  /* The single argument form of create_named_datasets() allows us to easily name data points. */
  ccplx.create_named_datasets({"x", "y", "f(x,y)"});

  /* Take a look at the generated cell complex */
  ccplx.dump_cplx(10);

  /* Cut out the tiny triangles on the unit circle! */
  tc_t::cull_cc_cells_on_domain_sdf_boundry(ccplx, unit_circle_sdf);

  /* We can do the above directly with MR_cell_cplx::cull_cells().  If you rewrote unit_circle_sdf(), you don't need to nest lambdas... */
  // ccplx.cull_cells([&ccplx](cc_t::cell_verts_t c){ return ccplx.cell_near_sdf_boundry(c, 
  //                                                                                [](cc_t::node_data_t pd) { return (tc_t::tsdf_to_csdf(unit_circle_sdf, 
  //                                                                                                                                     pd)); 
  //                                                                                                        }); 
  //                                          });

  /* Notice how it changed after the fold */
  ccplx.dump_cplx(10);

  ccplx.write_xml_vtk("surface_plot_step.vtu", "surface_plot_step");
}
/** @endcond */
