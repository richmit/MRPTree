// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      surface_branch_glue.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-16
 @brief     Mirroring & gluing surfaces together.@EOL
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

  In surface_plot_edge.cpp we encountered the unit sphere defined by the zeros of @f$1^2=x^2+y^2+z^2@f$, and the related function 
  @f$f(x,y)=\sqrt{1-x^2-y^2}@f$ obtained by "solving" for @f$z@f$.  Note that if @f$z=f(x, y)@f$, then both @f$z@f$ and @f$-z@f$ satisfy the original equation.
  While the square root function is positive by definition, we might wish to think of  @f$f(x, y)@f$ as a multi-valued function with two branches -- a
  positive one and a negative one.

  In simple cases like this, where the two branches are reflections across an axis plane, we can use `MR_cell_cplx::mirror()` to mirror the geometry and seal up
  any holes.  This is really the only change from surface_plot_edge.cpp.

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
tt_t::rrpt_t half_sphere(tt_t::drpt_t xvec) {
  double m = xvec[0] * xvec[0] + xvec[1] * xvec[1];
  if (m > 1) {
    return std::numeric_limits<double>::quiet_NaN();
  } else {
    return std::sqrt(1-m);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-1.2, -1.2}, 
            { 1.2,  1.2});
  cc_t ccplx;

  // Sample a uniform grid across the domain
  tree.refine_grid(5, half_sphere);

  /* Refine near the edge */
  tree.refine_recursive_if_cell_vertex_is_nan(6, half_sphere);

  tree.dump_tree(10);

  /* By passing half_sphere() to the construct_geometry_fans() we enable broken edges (an edge with one good point and one NaN) to be repaired. */
  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FDOMAIN, 0}, 
                                 {tc_t::val_src_spc_t::FDOMAIN, 1},
                                 {tc_t::val_src_spc_t::FRANGE,  0}},
                                half_sphere
                               );

  ccplx.create_named_datasets({"x", "y", "f(x,y)"},
                              {{"NORMALS", {0, 1, 2}}});

  /* This is the magic.  We add new cells with the third element of each point data vector negated. */
  ccplx.mirror({0, 0, 1});

  ccplx.dump_cplx(10);

  ccplx.write_xml_vtk("surface_branch_glue.vtu", "surface_branch_glue");
  ccplx.write_legacy_vtk("surface_branch_glue.vtk", "surface_branch_glue");
  ccplx.write_ply("surface_branch_glue.ply", "surface_branch_glue");
}
/** @endcond */
