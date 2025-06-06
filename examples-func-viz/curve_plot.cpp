// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      curve_plot.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-18
 @brief     A simple curve plot.@EOL
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

  Univariate function plots are the bread-and-butter of the plotting world.  Normally a simple, uniformly spaced, sequence is enough to get the job 
  done quite nicely.  Still, a few things can come up:

   - Jump discontinuities & Vertical asymptotes: Resolved with higher sampling near the discontinuities and a cutting edge (TBD)
   - Isolated, non-differentiable points:  Resolved with higher sampling near the points and a folding edge (TBD)
   - Undefined intervals:  Resolved with higher sampling near the edges and NaN edge repair
   - Regions of high oscillation: Resolved with higher sampling on the regions
   - Extrema: Resolved with higher sampling near the extrema

  Note that most of the items above are listed TBD.  A few features need to be added to `MR_rt_to_cc`. ;)  Note the TODO comments in the body of `main()`.
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b1d1rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t f(tt_t::drpt_t x) { 
  double ret = (x<0?-1:1)*std::pow(std::abs(x), 1/3.0) * std::sqrt((x+1.5)*(x+1.5)-1) * (x-2);
  if (x>2)
    ret = 2+std::sin(20*x);
  if (ret < -3)
    ret = -3;
  if (ret > 3.2)
    ret = 3.2;
  return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree(-3, 3);
  cc_t ccplx;

  // Sample a uniform grid across the domain
  tree.refine_grid(5, f);

  // Refine near NaN
  tree.refine_recursive_if_cell_vertex_is_nan(10, f);
  // TODO: Add NaN edge repair when implemented in MR_rt_to_cc

  // Refine near vertical tangent line
  tree.refine_leaves_recursive_cell_pred(10, f, [&tree](tt_t::diti_t i) { return (tree.cell_near_domain_point(0.0, 1.0e-2, i)); });
  // TODO: Use derivative test for this

  // Step discontinuities at 2.
  tree.refine_leaves_recursive_cell_pred(10, f, [&tree](tt_t::diti_t i) { return (tree.cell_near_domain_point(2.0, 1.0e-2, i)); });
  // TODO: Add cell cut when implemented in MR_rt_to_cc

  // Non differentiable point near x=-2.619185320
  tree.refine_leaves_recursive_cell_pred(11, f, [&tree](tt_t::diti_t i) { return (tree.cell_near_domain_point(-2.619185320, 1.0e-2, i)); });
  // TODO: Add folding edge when implemented in MR_rt_to_cc

  // High oscillation from [2,3]
  tree.refine_leaves_recursive_cell_pred(10, f, [&tree](tt_t::diti_t i) { return (tree.diti_to_drpt(i) >= 2.0); });

  // Extrema near -0.2171001290
  tree.refine_leaves_recursive_cell_pred(10, f, [&tree](tt_t::diti_t i) { return (tree.cell_near_domain_point(-0.2171001290, 1.0e-2, i)); });
  // TODO: Use derivative test for this

  // Extrema near 0.8775087009
  tree.refine_leaves_recursive_cell_pred(8, f, [&tree](tt_t::diti_t i) { return (tree.cell_near_domain_point(0.8775087009, 1.0e-2, i)); });
  // TODO: Use derivative test for this

  tree.dump_tree(10);

  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                1,
                                {{tc_t::val_src_spc_t::FDOMAIN,   0  }, 
                                 {tc_t::val_src_spc_t::FRANGE,    0  },
                                 {tc_t::val_src_spc_t::CONSTANT, 0.0}},
                                f
                               );

  // Note the first argument need not name *every* data element, just the first ones.
  ccplx.create_named_datasets({"x", "f(x)"});

  ccplx.dump_cplx(10);

  ccplx.write_xml_vtk("curve_plot.vtu", "curve_plot");
}
/** @endcond */
