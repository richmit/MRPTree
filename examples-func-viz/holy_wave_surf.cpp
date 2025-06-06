// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      holy_wave_surf.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-08-15
 @brief     Directly use MR_cell_cplx to randomly place triangles on a surface.@EOL
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
  Just a fun example directly using `MR_cell_cplx` to drop randomly placed triangles on a surface.
*/
/*******************************************************************************************************************************************************.H.E.**/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <random>                                                        /* C++ random numbers      C++11    */
#include <chrono>                                                        /* time                    C++11    */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_cell_cplx.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double f(double x, double y) {
  double d = x*x+y*y;
  double z = std::exp(-d/4)*std::cos(4*std::sqrt(d));
  return z;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  std::cout << "PROGRAM: START" << std::endl;
  std::chrono::time_point<std::chrono::system_clock> program_start_time = std::chrono::system_clock::now();

  mjr::MRccT5 aPoly;

  aPoly.create_dataset_to_point_mapping({0, 1, 2});
  aPoly.create_named_datasets({"x", "y", "z", "za", "xyMag", "xyzMag", "dir"});

  const int    num_tri =  20000;
  const double tri_siz =  0.07;
  const double x_min   = -2.1;
  const double x_max   =  2.1;
  const double y_min   = -2.1;
  const double y_max   =  2.1;

  std::random_device rd;
  std::mt19937 rEng(rd());
  std::uniform_real_distribution<double> x_uniform_dist_float(x_min, x_max);
  std::uniform_real_distribution<double> y_uniform_dist_float(y_min, y_max);
  std::uniform_real_distribution<double> siz_uniform_dist_float(-tri_siz, tri_siz);

  std::cout << "SAMPLE: START" << std::endl;  
  std::chrono::time_point<std::chrono::system_clock> sample_start_time = std::chrono::system_clock::now();
  for(int i=0; i<num_tri; i++) {
    double xc = x_uniform_dist_float(rEng);
    double yc = y_uniform_dist_float(rEng);

    double x1 = xc + siz_uniform_dist_float(rEng);
    double y1 = yc + siz_uniform_dist_float(rEng);
    double z1 = f(x1, y1);
    double x2 = xc + siz_uniform_dist_float(rEng);
    double y2 = yc + siz_uniform_dist_float(rEng);
    double z2 = f(x2, y2);
    double x3 = xc + siz_uniform_dist_float(rEng);
    double y3 = yc + siz_uniform_dist_float(rEng);
    double z3 = f(x3, y3);

    double xa = (x1+x2+x3)/3;
    double ya = (y1+y2+y3)/3;
    double za = (z1+z2+z3)/3;

    double xy_mag  = std::sqrt(xa*xa+ya*ya);
    double xyz_mag = std::sqrt(xa*xa+ya*ya+za*za);
    
    mjr::MRccT5::node_idx_t p1i = aPoly.add_node({x1, y1, z1, za, xy_mag, xyz_mag});
    mjr::MRccT5::node_idx_t p2i = aPoly.add_node({x2, y2, z2, za, xy_mag, xyz_mag});
    mjr::MRccT5::node_idx_t p3i = aPoly.add_node({x3, y3, z3, za, xy_mag, xyz_mag});

    aPoly.add_cell(mjr::MRccT5::cell_kind_t::TRIANGLE, {p1i, p2i, p3i});
  }

  std::cout << "SAMPLE: Total Points: " << aPoly.node_count() << std::endl;
  std::cout << "SAMPLE: Total Cells: " << aPoly.num_cells() << std::endl;

  std::chrono::duration<double> sample_run_time = std::chrono::system_clock::now() - sample_start_time;
  std::cout << "SAMPLE: Total Runtime " << sample_run_time.count() << " sec" << std::endl;
  std::cout << "SAMPLE: END" << std::endl;  

  std::cout << "XML WRITE: START" << std::endl;
  std::chrono::time_point<std::chrono::system_clock> xwrite_start_time = std::chrono::system_clock::now();
  int xwrite_result = aPoly.write_xml_vtk("holy_wave_surf.vtu", "holy_wave_surf");
  std::chrono::duration<double> xwrite_run_time = std::chrono::system_clock::now() - xwrite_start_time;
  std::cout << "XML WRITE: Total Runtime " << xwrite_run_time.count() << " sec" << std::endl;
  std::cout << "XML WRITE: END -- " << (xwrite_result != 0 ? "BAD" : "GOOD") << std::endl;

  std::chrono::duration<double> program_run_time = std::chrono::system_clock::now() - program_start_time;
  std::cout << "PROGRAM: Total Runtime " << program_run_time.count() << " sec" << std::endl;
  std::cout << "PROGRAM: END" << std::endl;

  return 0;
}

