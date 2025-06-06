// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      ear_surface.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Sampling for an implicit surface.@EOL
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

  This example is very similar to implicit_surface.cpp; however, instead of extracting a surface from a quad tessellation of a hexahedrona we extract the
  surface from a tessellation of a pyramids.  Why use pyramids instead of hexahedrona?  In the example implicit_surface.cpp all of the underlying cells are
  the same size, and thus no gaps occur in the extracted level set (the surface).  In this example we have cells that vary a great deal in size.

  This example demonstrates scratch made cell predicates, and uses them to increse sample resolution on parts of the surface that are particularly difficlut
  to extract.

  This surface is defined by the zeros of the following polynomial

  @f[ x^2-y^2 z^2+z^3 @f]
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b3d1rT            tt_t;
typedef mjr::MRccT9                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t isf(tt_t::drpt_t xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double z = xvec[2];
  return x*x-y*y*z*z+z*z*z;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t besdf(tt_t::drpt_t xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double z = xvec[2];
  return x*(z-y*y); // Inner ear edge
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree;
  cc_t ccplx;

  /* Initial uniform sample */
  tree.refine_grid(3, isf);

  /* Refine near surface */
  tree.refine_leaves_recursive_cell_pred(6, isf, [&tree](tt_t::diti_t i) { return (tree.cell_cross_sdf(i, isf)); });

  /* Refine between the ears.  Could have also just refined on x=0 plane.  */
  tree.refine_leaves_recursive_cell_pred(8, isf, [&tree](tt_t::diti_t i) { auto x = tree.diti_to_drpt(i); return (std::abs(x[1])<0.5) && (tree.cell_cross_sdf(i, besdf)); });

  /* Refine on the x-y plane */
  tree.refine_leaves_recursive_cell_pred(7, isf, [&tree](tt_t::diti_t i) { return (tree.cell_near_domain_level(i, 2, 0.0, 1.0e-6)); });

  /* Balance the tree */
  tree.balance_tree(1, isf);

  tree.dump_tree(5);

  /* Convert our tree to a cell complex.  Note that we use an SDF to export only cells that contain our surface */
  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                tree.get_leaf_cells_pred(tree.ccc_get_top_cell(), [&tree](tt_t::diti_t i) { return (tree.cell_cross_sdf(i, isf)); }),
                                3,
                                {{tc_t::val_src_spc_t::FDOMAIN, 0}, 
                                 {tc_t::val_src_spc_t::FDOMAIN, 1},
                                 {tc_t::val_src_spc_t::FDOMAIN, 2}});

  /* Name the data points */
  ccplx.create_named_datasets({"x", "y", "z", "f(x,y,z)"});
  
  /* Display some data about the cell complex */
  ccplx.dump_cplx(5);

  /* Write out our cell complex */
  ccplx.write_xml_vtk("ear_surface.vtu", "ear_surface");
}
/** @endcond */
