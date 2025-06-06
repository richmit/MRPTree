// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      parametric_surface_with_defects.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Parametric surface with defects.@EOL
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

  This example illustrates some of the things that can go wrong when generating parametric surfaces.  We dump two version of the tessellation -- one with
  quads and one with triangles.  This allows us to better illustrate how some defects show up.
   - Quads that are not plainer.
     Look closely at the rectangular tessellation, and note the "rectangles" appear to be broken in across the diagonal -- at least that's how they appear in
     most tools including Paraview & meshlab.
   - At the poles, the rectangular cells of the tree map to three distinct points instead of four.
     This means for the rectangular tessellation, the rectangles at the poles are degenerate!  Then converting from tree to cell complex, these quads are
     removed because they are degenerate. This is not an issue for the triangular tessellation (FANS).
   - The v=0 edge meets up with the v=1 edge.  
     Because we have `chk_point_unique` set to true for the cell complex object, the duplicate points are "welded" together when the points are added.  This
     results in a seamless edge.
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <numbers>                                                       /* C++ math constants      C++20    */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d3rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t par_sphere(tt_t::drpt_t xvec) {
  double u = std::numbers::pi/4 * xvec[0] + std::numbers::pi/4;
  double v = std::numbers::pi   * xvec[1] + std::numbers::pi;
  return { std::sin(u)*std::cos(v),
           std::sin(u)*std::sin(v),
           std::cos(u)
         };
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree;
  cc_t ccplx;

  /* Uniform sampling */
  tree.refine_grid(6, par_sphere);


  /* First we dump a tessellation made of triangles */
  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FRANGE, 0},
                                 {tc_t::val_src_spc_t::FRANGE, 1},
                                 {tc_t::val_src_spc_t::FRANGE, 2}});
  ccplx.create_named_datasets({"u", "v", "x(u,v)", "y(u,v)", "z(u,v)"});
  ccplx.dump_cplx(5);
  ccplx.write_xml_vtk("parametric_surface_with_defects-tri.vtu", "parametric_surface_with_defects-tri");

  /* Next we dump a tessellation made of rectangles */
  ccplx.clear(); // We need to clear out the old contents first!
  tc_t::construct_geometry_rects(ccplx,
                                 tree,
                                 2,
                                 {{tc_t::val_src_spc_t::FRANGE, 0},
                                  {tc_t::val_src_spc_t::FRANGE, 1},
                                  {tc_t::val_src_spc_t::FRANGE, 2}});
  ccplx.create_named_datasets({"u", "v", "x(u,v)", "y(u,v)", "z(u,v)"});
  ccplx.dump_cplx(5);
  ccplx.write_xml_vtk("parametric_surface_with_defects-rect.vtu", "parametric_surface_with_defects-rect");
}
/** @endcond */
