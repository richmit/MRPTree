// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      hello_world_mraster.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Minimal example for MR_rect_tree/MR_cell_cplx/MR_rt_to_cc and MRaster.@EOL
 @keywords  surface plot 2d 3d color
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
  This example is almost identical to hello_world.cpp except:
   - We use MRaster to generate a color pallet (MRaster has a great many pallets: https://richmit.github.io/mraster/ColorSchemes.html)
   - We directly embed the colors into the resulting MR_cell_cplx 
   - We identify the color vector in the .VTU file.
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"
#include "MRcolor.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d4rT            tt_t;
typedef mjr::MRccT5                  cc_t;
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;
typedef mjr::color3c64F              ct_t; // This is a 3 channel color with 64-bit floating point channels (i.e. a 192-bit color type)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t damp_cos_wave(tt_t::drpt_t xvec) {
  double x = xvec[0];
  double y = xvec[1];
  double d = x*x+y*y;
  double z = std::exp(-d/4)*std::cos(4*std::sqrt(d));
  ct_t   c = ct_t::csPLYviridis::c((z+0.87)/1.87);
  return {z, c.getRed(), c.getGreen(), c.getBlue()};  // We can grab the raw channel values with functions like getRed() because the channels are doubles.  
}                                                     // If we had used an integer channel color (like color3c8b), we would have needed to use getRed_dbl().


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-2.1, -2.1}, 
            { 2.1,  2.1});
  cc_t ccplx;

  tree.refine_grid(7, damp_cos_wave);

  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FDOMAIN, 0}, 
                                 {tc_t::val_src_spc_t::FDOMAIN, 1},
                                 {tc_t::val_src_spc_t::FRANGE,  0}});
  ccplx.create_named_datasets({"x", "y", "f(x,y)", "c_r(x,y)", "c_g(x,y)", "c_b(x,y)"},
                             {{"COLORS", {3, 4, 5}}});

  ccplx.write_xml_vtk("hello_MRaster.vtu", "hello_MRaster");
}
/** @endcond */
