// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      hello_world_tree_adaptive.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Minimal example for MR_rect_tree with adaptive sampling.@EOL
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

  This program writes a simple tabular file containing the tree's data -- one sample per line with space separated coordinates and values.  These
  may be interrupted as points in 3D, and graphed by tools like GNU Plot.  On Windows running MSYS2, the following command sequence will produce a
  graph on the screen:

     make hello_world_tree_adaptive
     ./hello_world_tree_adaptive.exe
     gnuplot.exe ../examples/hello_world_tree_adaptive.gp

*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MR_rect_tree.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef mjr::tree15b2d5rT tt_t;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
tt_t::rrpt_t damp_cos_wave(tt_t::drpt_t xvec) {
  double x  = xvec[0];
  double y  = xvec[1];
  double d  = x*x+y*y;
  double m  = std::exp(-d/4);
  double s  = std::sqrt(d);
  double z  = m*cos(4*s);
  double dd = -m*(cos(4*s)*s+8*sin(4*s));
  if (s>1.0e-5) 
    dd = dd / (4 * s);
  else 
    dd = 1;
  return {z, dd};
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main() {
  tt_t tree({-2.1, -2.1}, 
            { 2.1,  2.1});
  
  // Make a few samples on a uniform grid
  tree.refine_grid(2, damp_cos_wave);

  // Use the directional derivative radiating from the origin to deep sample on the humps
  tree.refine_leaves_recursive_cell_pred(6, damp_cos_wave, [&tree](tt_t::diti_t i) { return tree.cell_cross_range_level(i, 1, 0.0); });

  // Balance the three to the traditional level of 1 (no  cell borders a cell more than half it's size)
  tree.balance_tree(1, damp_cos_wave);

  tree.dump_tree_datafile("hello_world_tree_adaptive.tab");
}
/** @endcond */
