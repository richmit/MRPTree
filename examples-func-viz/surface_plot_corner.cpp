// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      surface_plot_corner.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @brief     Surface with a sharp edge.@EOL
 @date      2024-07-16
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

  The function illustrated here is continuous on the entire plane, but has no derivative on the unit circle.  While no derivative exists on the unit circle,
  directional derivatives pointing from the origin approach infinity as we get close to the unit circle.  The derivative at the origin is zero.  Thus the
  surface is not only zero on the unit circle, but it drops to zero very quickly from it's local extrema at the origin.

  If we sample on a uniform grid, some of the resulting polygons will have vertexes both inside and outside the unit circle.  These polygons will never touch
  the x-y plane, and thus the surface will not appear to have a uniform zero set on the unit circle.  At low resolution the results are so bad they are
  difficult to interpret.  At higher resolutions we see what appears to be a jagged edge over the unit circle.  Meaning the results are visually quite wrong,
  but an astute viewer might well guess the true behavior of the function from the resulting image.  In order to correct this graph we need sample points in
  the triangulation that are on, or very near, the unit circle.  We can do that by folding and resampling the cell complex on the unit circle.

    - How to drive up the sample rate near a particular SDF -- so that we get higher resolution where the surface meets the plane.
    - How to "fold" the resulting triangles to achieve higher accuracy on the non-differentiable edge.
    - How to use `tsampf_to_cdatf()` & `tsdf_to_csdf()` to adapt functions designed for a `MR_rect_tree` for use with a `MR_cell_cplx`.
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

tt_t::rrpt_t half_sphere_hat(tt_t::drpt_t xvec) {
  double m = xvec[0] * xvec[0] + xvec[1] * xvec[1];
  return (std::sqrt(std::abs(1-m)));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::src_t unit_circle_sdf(tt_t::drpt_t xvec) {
  double m = xvec[0] * xvec[0] + xvec[1] * xvec[1];
  return (1-m);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-1.1, -1.1}, 
            { 1.1,  1.1});
  cc_t ccplx;

  /* Here is another way to get fine samples on the circle, but with a SDF this time. */
  tree.refine_grid(5, half_sphere_hat);

  /* Increase sample resolution on the unit circle.  Here we do that with an SDF. */
  tree.refine_leaves_recursive_cell_pred(7, half_sphere_hat, [&tree](int i) { return (tree.cell_cross_sdf(i, unit_circle_sdf)); });

  /* Balance the three to the traditional level of 1 (no cell borders a cell more than half it's size) */
  tree.balance_tree(1, half_sphere_hat);

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

  /* Fold all triangles that cross the unit circle! */
  ccplx.triangle_folder([](cc_t::node_data_t x){return tc_t::tsampf_to_cdatf(half_sphere_hat, x); }, 
                        [](cc_t::node_data_t x){return tc_t::tsdf_to_csdf(unit_circle_sdf,    x); });

  /* Notice how it changed after the fold */
  ccplx.dump_cplx(10);

  ccplx.write_xml_vtk("surface_plot_corner.vtu", "surface_plot_corner");
}
/** @endcond */
