// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      ear_surface_glue.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-08-18
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

  This example illustrates the same object as ear_surface.cpp, but uses a different strategy to generate the triangulation.  For reference, the surface we
  are interested in is defined by the zeros of the following polynomial:

  @f[ x^2-y^2 z^2+z^3 @f]

  In ear_surface.cpp we used a mesh of cubes covering a 3D region, and extracted the triangulation as a contour surface.  In this example, we solve this equation
  for @f$x@f$ giving us two surfaces:

  @f[ +z\cdot\sqrt{y^2-z} @f]
  @f[ -z\cdot\sqrt{y^2-z} @f]

  As in the example surface_plot_edge.cpp, we use the NaN solver to flesh the surface out so that it hits the X plane.  And as in surface_branch_glue.cpp, we
  then glue the two surfaces together.
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
tt_t::rrpt_t ear_yz(tt_t::drpt_t xvec) {
  tt_t::src_t y = xvec[0];
  tt_t::src_t z = xvec[1];
  tt_t::src_t m = y*y-z;
  if (m < 0) {
    return std::numeric_limits<double>::quiet_NaN();
  } else {
    return std::sqrt(m)*z;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree;
  cc_t ccplx;

  // Sample a uniform grid across the domain
  tree.refine_grid(7, ear_yz);

  /* Refine near the edge */
  tree.refine_recursive_if_cell_vertex_is_nan(8, ear_yz);

  tree.dump_tree(10);

  /* By passing half_sphere() to the construct_geometry_fans() we enable broken edges (an edge with one good point and one NaN) to be repaired. */
  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FRANGE,  0},
                                 {tc_t::val_src_spc_t::FDOMAIN, 0}, 
                                 {tc_t::val_src_spc_t::FDOMAIN, 1}},
                                ear_yz
                               );

  ccplx.create_named_datasets({"y", "z", "x=f(x,z)"});

  /* This is the magic.  We add new cells with the third element of each point data vector negated. */
  ccplx.mirror({0, 0, 1}, 1.0e-5);

  ccplx.dump_cplx(10);

  ccplx.write_xml_vtk("ear_surface_glue.vtu", "ear_surface_glue");
}
/** @endcond */
