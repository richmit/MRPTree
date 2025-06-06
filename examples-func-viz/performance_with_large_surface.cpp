// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      performance_with_large_surface.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-23
 @brief     Stress test with a large surface object.@EOL
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

  Just a nice parametric surface without any weirdness.  Some things demonstrated:
   - How to time various operations. 
     - Try with a large mesh (use a 9 in `refine_grid`).
     - Try reducing the number of data variables stored in the cell complex
     - Try removing the normal vector from the output
     - Try both `MRccT5` & `MRccF5` for `cc_t`
   - How to include a synthetic value that can be used for color mapping --  @f$ c(u,v) @f$ can be used to render stripes on the surface.
   - How to compute a normal to a parametric surface.  If the surface is defined by 
     @f[ \vec{f}(u,v)=(x(u,v), y(u,v), z(u,v)) @f]
     Then the normal vector is given by:
     @f[ \vec{n}=\frac{\partial \vec{f}}{\partial u}\times\frac{\partial \vec{f}}{\partial v} @f]
     Note that @f$ \frac{\partial \vec{f}}{\partial u} = \left[ \frac{\partial x}{\partial u}, \frac{\partial y}{\partial u}, \frac{\partial z}{\partial u} \right] @f$ &
      @f$ \frac{\partial \vec{f}}{\partial v} = \left[ \frac{\partial x}{\partial v}, \frac{\partial y}{\partial v}, \frac{\partial z}{\partial v} \right] @f$
   - How to include a normal in the cell complex
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <chrono>                                                        /* time                    C++11    */
#include <numbers>                                                       /* C++ math constants      C++20    */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"
#include "MR_cell_cplx.hpp"
#include "MR_rt_to_cc.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d15rT           tt_t;
typedef mjr::MRccT5                  cc_t;   // Replace with mjr::MRccF5, and compare bridge performance.
typedef mjr::MR_rt_to_cc<tt_t, cc_t> tc_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t stripy_shell(tt_t::drpt_t xvec) {
  double u    = std::numbers::pi   * xvec[0] + std::numbers::pi + 0.1; // U transformed from unit interval
  double v    = std::numbers::pi/2 * xvec[1] + std::numbers::pi/2;     // V transformed from unit interval
  double x    = u*std::sin(u)*std::cos(v);                             // X
  double y    = u*std::cos(u)*std::cos(v);                             // Y
  double z    = u*std::sin(v);                                         // Z
  double c    = std::fmod(u*sin(v), 2);                                // Stripes
  double dxdu = std::sin(u)*std::cos(v)+u*std::cos(u)*std::cos(v);     // dX/du
  double dxdv = -u*std::sin(u)*std::sin(v);                            // dX/dv
  double dydu = std::cos(u)*std::cos(v)-u*std::sin(u)*std::cos(v);     // dY/du
  double dydv = -u*std::cos(u)*std::sin(v);                            // dY/dv
  double dzdu = std::sin(v);                                           // dZ/du
  double dzdv = u*std::cos(v);                                         // dZ/dv
  double nx   = dydu*dzdv-dydv*dzdu;                                   // normal_X     This noraml      
  double ny   = dxdv*dzdu-dxdu*dzdv;                                   // normal_Y     will not be of 
  double nz   = dxdu*dydv-dxdv*dydu;                                   // normal_Z     unit length
  return {x, y, z, c, dxdu, dxdv, dydu, dydv, dzdu, dzdv, nx, ny, nz};
}                          

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  std::chrono::time_point<std::chrono::system_clock> start_time = std::chrono::system_clock::now();
  tt_t tree;
  cc_t ccplx;
  std::chrono::time_point<std::chrono::system_clock> construct_time = std::chrono::system_clock::now();

  tree.refine_grid(7, stripy_shell);
  std::chrono::time_point<std::chrono::system_clock> sample_time = std::chrono::system_clock::now();

  tree.dump_tree(10);
  std::chrono::time_point<std::chrono::system_clock> tdump_time = std::chrono::system_clock::now();

  tc_t::construct_geometry_fans(ccplx,
                                tree,
                                2,
                                {{tc_t::val_src_spc_t::FRANGE,  0},
                                 {tc_t::val_src_spc_t::FRANGE,  1},
                                 {tc_t::val_src_spc_t::FRANGE,  2}});
  std::chrono::time_point<std::chrono::system_clock> fan_time = std::chrono::system_clock::now();

  ccplx.create_named_datasets({"u", "v", 
                               "x(u,v)", "y(u,v)", "z(u,v)",
                               "c(u,v)", 
                               "dx(u,v)/du", "dx(u,v)/dv", "dy(u,v)/du", "dy(u,v)/dv", "dz(u,v)/du", "dz(u,v)/dv",
                               "nx", "ny", "nz"}, 
                              {{"NORMALS", {12, 13, 14}}});
  std::chrono::time_point<std::chrono::system_clock> dat_anno_time = std::chrono::system_clock::now();

  ccplx.dump_cplx(10);
  std::chrono::time_point<std::chrono::system_clock> cdump_time = std::chrono::system_clock::now();

  ccplx.write_xml_vtk("performance_with_large_surface.vtu", "performance_with_large_surface");
  std::chrono::time_point<std::chrono::system_clock> write_time = std::chrono::system_clock::now();

  std::cout << "construct_time time .. " << static_cast<std::chrono::duration<double>>(construct_time-start_time)   << " sec" << std::endl;
  std::cout << "sample_time time ..... " << static_cast<std::chrono::duration<double>>(sample_time-construct_time)  << " sec" << std::endl;
  std::cout << "tree dump time ....... " << static_cast<std::chrono::duration<double>>(tdump_time-sample_time)      << " sec" << std::endl;
  std::cout << "bridge time .......... " << static_cast<std::chrono::duration<double>>(fan_time-tdump_time)         << " sec" << std::endl;
  std::cout << "dataset anno time .... " << static_cast<std::chrono::duration<double>>(dat_anno_time-fan_time)      << " sec" << std::endl;
  std::cout << "complex dump time .... " << static_cast<std::chrono::duration<double>>(cdump_time-dat_anno_time)    << " sec" << std::endl;
  std::cout << "write_vtk time ....... " << static_cast<std::chrono::duration<double>>(write_time-cdump_time)       << " sec" << std::endl;
  std::cout << "Total Run _time ...... " << static_cast<std::chrono::duration<double>>(write_time-start_time)       << " sec" << std::endl;
}
/** @endcond */
