// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      parametric_curve_3d.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Parametric curve as the intersection of two parametric surfaces.@EOL
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

  This program produces an interesting visualization of an object known as the twisted cubic.  In parametric form, the curve may be expressed as

  @f[ f(t)=[t, t^2, t^3] @f]

  Alternately the curve is also the intersection of two surfaces in @f$\mathbb{R}^3@f$:

  @f[ y=f_2(x, z)=y^2 @f]
  @f[ z=f_3(x, y)=x^3 @f]

  The "typical" way to graph a surface like @f$f_2@f$ is to transform it into pseudo-parametric form.  In Maple that might look like this

  \verbatim                                                                                                              
  plot3d([u, u^2, v], u=-1..1, v=-1..1):
  \endverbatim

  We could do that with `MRPTree`, but it is easier to simply map the variables when we use `construct_geometry_fans()`.

  Another interesting use of `MRPTree` in this example is the way we have transformed each surface function into an SDF to drive up sample resolution near the
  surface intersection.  This would allow us to use a tool like Paraview to compute an approximation to the the intersection.  Just in case the reader is
  not using a tool that can extract a nice surface intersection, I have also dumped the curve out in a 3rd `.VTU` file.

*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b1d3rT              tt1_t;
typedef mjr::MRccT5                    cc1_t;
typedef mjr::MR_rt_to_cc<tt1_t, cc1_t> tc1_t;

typedef mjr::tree15b2d1rT              tt2_t;
typedef mjr::MRccT5                    cc2_t;
typedef mjr::MR_rt_to_cc<tt2_t, cc2_t> tc2_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt1_t::rrpt_t twisted_cubic_crv(tt1_t::drpt_t t) {
  return { t, t*t, t*t*t };
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt2_t::rrpt_t twisted_cubic_srf1(tt2_t::drpt_t xzvec) {
  tt2_t::src_t x = xzvec[0];
  return x*x;
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt2_t::src_t twisted_cubic_srf1_sdf(tt2_t::drpt_t xzvec) {
  tt2_t::src_t z = xzvec[1];
  return (twisted_cubic_srf1(xzvec)-z);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt2_t::rrpt_t twisted_cubic_srf2(tt2_t::drpt_t xyvec) {
  tt2_t::src_t x = xyvec[0];
  return x*x*x;
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt2_t::src_t twisted_cubic_srf2_sdf(tt2_t::drpt_t xyvec) {
  tt2_t::src_t y = xyvec[1];
  return (twisted_cubic_srf2(xyvec)-y);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt1_t crv_tree;
  cc1_t crv_ccplx;
  crv_tree.refine_grid(8, twisted_cubic_crv);
  tc1_t::construct_geometry_fans(crv_ccplx,
                                 crv_tree,
                                 1,
                                 {{tc1_t::val_src_spc_t::FRANGE, 0},
                                  {tc1_t::val_src_spc_t::FRANGE, 1},
                                  {tc1_t::val_src_spc_t::FRANGE, 2}});
  crv_ccplx.create_named_datasets({"t", "x(t)", "y(t)", "z(t)"});
  crv_ccplx.dump_cplx(5);
  crv_ccplx.write_xml_vtk("parametric_curve_3d-crv.vtu", "parametric_curve_3d-crv");

  tt2_t srf1_tree;
  cc2_t srf1_ccplx;
  srf1_tree.refine_grid(5, twisted_cubic_srf1);
  srf1_tree.refine_leaves_recursive_cell_pred(6, twisted_cubic_srf1, [&srf1_tree](tt2_t::diti_t i) { return srf1_tree.cell_cross_sdf(i, twisted_cubic_srf2_sdf); });
  srf1_tree.balance_tree(1, twisted_cubic_srf1);
  tc2_t::construct_geometry_fans(srf1_ccplx,
                                 srf1_tree,
                                 2,
                                 {{tc2_t::val_src_spc_t::FDOMAIN, 0},
                                  {tc2_t::val_src_spc_t::FRANGE,  0},
                                  {tc2_t::val_src_spc_t::FDOMAIN, 1}});
  srf1_ccplx.create_named_datasets({"u", "v", "x(u,v)", "y(u,v)", "z(u,v)"});
  srf1_ccplx.dump_cplx(5);
  srf1_ccplx.write_xml_vtk("parametric_curve_3d-srf1.vtu", "parametric_curve_3d-srf1");

  tt2_t srf2_tree;
  cc2_t srf2_ccplx;
  srf2_tree.refine_grid(5, twisted_cubic_srf2);
  srf2_tree.refine_leaves_recursive_cell_pred(6, twisted_cubic_srf2, [&srf2_tree](tt2_t::diti_t i) { return srf2_tree.cell_cross_sdf(i, twisted_cubic_srf1_sdf); });
  srf2_tree.balance_tree(1, twisted_cubic_srf2);
  tc2_t::construct_geometry_fans(srf2_ccplx,
                                 srf2_tree,
                                 2,
                                 {{tc2_t::val_src_spc_t::FRANGE, 0},
                                  {tc2_t::val_src_spc_t::FRANGE, 1},
                                  {tc2_t::val_src_spc_t::FRANGE, 2}});
  srf2_ccplx.create_named_datasets({"u", "v", "x(u,v)", "y(u,v)", "z(u,v)"});
  srf2_ccplx.dump_cplx(5);
  srf2_ccplx.write_xml_vtk("parametric_curve_3d-srf2.vtu", "parametric_curve_3d-srf2");
}
/** @endcond */
