// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      recipe-surf-plot-norm.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-10-02
 @brief     Cookbook recipie: recipe-surf-plot-norm.@EOL
 @keywords  surface plot 2d 3d regular sample
 @std       C++23
 @see
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
  double dd = -m*(cos(4*s)*s+8*sin(4*s));
  double nm = std::sqrt(1+dx*dx+dy*dy);
  if (s>1.0e-5) {
    dx = dx / s;
    dy = dy / s;
    dd = dd / (4 * s);
  } else {
    dx = 1;
    dy = 1;
    dd = 1;
    nm = 1;
  }
  return {z, -dx/nm, -dy/nm, 1/nm, dd};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-2.1, -2.1},
            { 2.1,  2.1});
  cc_t ccplx;

  tree.refine_grid(4, damp_cos_wave);
  tree.dump_tree(5);

  tree.refine_leaves_recursive_cell_pred(6, damp_cos_wave, [&tree](tt_t::diti_t i) { return tree.cell_cross_range_level(i, 4, 0.0); });
  tree.balance_tree(1, damp_cos_wave);
  tree.dump_tree(5);

  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FDOMAIN, 0},
                                 {tc_t::val_src_spc_t::FDOMAIN, 1},
                                 {tc_t::val_src_spc_t::FRANGE,  0}});
  ccplx.set_data_name_to_data_idx_lst({{"x",        {0}},
                                       {"y",        {1}},
                                       {"z=f(x,y)", {2}},
                                       {"ddiv",     {6}},
                                       {"NORMALS",  {3,4,5}}});

  ccplx.write_xml_vtk("recipe-surf-plot-norm.vtu", "recipe-surf-plot-norm");
}
/** @endcond */
