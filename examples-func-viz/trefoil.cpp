// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      trefoil.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Trefoil parametric surface.@EOL
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

  This example doesn't really demonstrate anything not found in the other examples.  It's just a neat surface. ;)
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
typedef mjr::tree15b2d6rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t trefoil(tt_t::drpt_t xvec) {
  double u = xvec[0] * std::numbers::pi;
  double v = xvec[1] * std::numbers::pi;
  double r = 5;
  double x = r * std::sin(3 * u) / (2 + std::cos(v));
  double y = r * (std::sin(u) + 2 * std::sin(2 * u)) / (2 + std::cos(v + std::numbers::pi * 2 / 3));
  double z = r / 2 * (std::cos(u) - 2 * std::cos(2 * u)) * (2 + std::cos(v)) * (2 + std::cos(v + std::numbers::pi * 2 / 3)) / 4;
  double dxdu = (3*r*std::cos(3*u))/(std::cos(v)+2);
  double dxdv = (r*std::sin(3*u)*std::sin(v))/(std::cos(v)+2)/(std::cos(v)+2);
  double dydu = (r*(4*std::cos(2*u)+std::cos(u)))/(std::cos(v+(2*std::numbers::pi)/3)+2);
  double dydv = (r*(2*std::sin(2*u)+std::sin(u))*std::sin(v+(2*std::numbers::pi)/3))/
    ((std::cos(v+(2*std::numbers::pi)/3)+2)*(std::cos(v+(2*std::numbers::pi)/3)+2));
  double dzdu = (r*(4*std::sin(2*u)-std::sin(u))*(std::cos(v)+2)*(std::cos(v+(2*std::numbers::pi)/3)+2))/8;
  double dzdv = (-(r*(std::cos(u)-2*std::cos(2*u))*(std::cos(v)+2)*std::sin(v+(2*std::numbers::pi)/3))/8) -
    (r*(std::cos(u)-2*std::cos(2*u))*std::sin(v)*(std::cos(v+(2*std::numbers::pi)/3)+2))/8;
  double nx   = dydu*dzdv-dydv*dzdu;
  double ny   = dxdv*dzdu-dxdu*dzdv;
  double nz   = dxdu*dydv-dxdv*dydu;
  double nm   = std::sqrt(nx*nx+ny*ny+nz*nz);
  nm = (nm > 0 ? nm : 1);
  nx = nx / nm;
  ny = ny / nm;
  nz = nz / nm;
  return {x, y, z, nx, ny, nz};
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree;

  cc_t ccplx;

  tree.refine_grid(7, trefoil);

  tree.dump_tree(20);

  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FRANGE,  0},
                                 {tc_t::val_src_spc_t::FRANGE,  1},
                                 {tc_t::val_src_spc_t::FRANGE,  2}});

  ccplx.create_named_datasets({"u", "v", "x(u,v)", "y(u,v)", "z(u,v)", "nx", "ny", "nz"},
                              {{"NORMALS", {5, 6, 7}}});

  ccplx.write_xml_vtk("trefoil.vtu", "trefoil");
}
/** @endcond */
