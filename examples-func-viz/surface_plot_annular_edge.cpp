// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      surface_plot_annular_edge.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-16
 @brief     Surface with an undefined annular region.@EOL
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

  In example surface_plot_edge.cpp we had a surface plot with a large, undefined region easily discovered by Maple & Nathematica.  In this example we have a
  surface with a smaller, annular undefined region that is harder for automatic software to get right.  The function in question is as follows:

    @f[ \sqrt{\sqrt{\vert 1 - x^2 - y^2\vert} - \frac{3}{20}} @f]

  This example differs from surface_plot_edge.cpp in that we direct additional sampling near the unit circle.  Other than that, it's just about identical.
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
tt_t::rrpt_t annular_hat(tt_t::drpt_t xvec) {
  double v = std::sqrt(std::abs(1 - xvec[0] * xvec[0] - xvec[1] * xvec[1])) - 0.15;
  if (v < 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  } else {
    return std::sqrt(v);
  }
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

  // Sample a uniform grid across the domain
  tree.refine_grid(5, annular_hat);

  /* annular_hat produces NaNs in a thin annulus near the unit circle.*/
  tree.refine_leaves_recursive_cell_pred(7, annular_hat, [&tree](int i) { return (tree.cell_cross_sdf(i, unit_circle_sdf)); });

  /* Balance the three to the traditional level of 1 (no cell borders a cell more than half it's size) */
  tree.balance_tree(1, annular_hat);

  tree.dump_tree(10);

  /* By passing annular_hat() to the construct_geometry_fans() we enable broken edges (an edge with one good point and one NaN) to be repaired. */
  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FDOMAIN, 0}, 
                                 {tc_t::val_src_spc_t::FDOMAIN, 1},
                                 {tc_t::val_src_spc_t::FRANGE,  0}},
                                annular_hat
                               );

  ccplx.create_named_datasets({"x", "y", "f(x,y)"},
                              {{"NORMALS", {0, 1, 2}}});

  ccplx.dump_cplx(10);

  ccplx.write_xml_vtk("surface_plot_annular_edge.vtu", "surface_plot_annular_edge");
}
/** @endcond */
