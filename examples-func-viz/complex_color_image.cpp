// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      complex_color_image.cpp
 @author    Mitch Richling <https://www.mitchr.me>
 @brief     2D complex function plot.@EOL
 @std       C++20
 @copyright
  @parblock
  Copyright (c) 1988-2015, Mitchell Jay Richling <https://www.mitchr.me> All rights reserved.

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
*/
/*******************************************************************************************************************************************************.H.E.**/
/** @cond exj */

//--------------------------------------------------------------------------------------------------------------------------------------------------------------
#include "ramCanvas.hpp"

//--------------------------------------------------------------------------------------------------------------------------------------------------------------
using cplx = std::complex<double>;
using ct_t = mjr::ramCanvas3c8b::colorType;

//--------------------------------------------------------------------------------------------------------------------------------------------------------------
cplx f(cplx z) {
  if ( (std::abs(z-1.0) > 1.0e-5) && (std::abs(z+1.0) > 1.0e-5) ) 
    return 1.0/(z+1.0) + 1.0/(z-1.0);
  else
    return 0;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------
int main(void) {
  std::chrono::time_point<std::chrono::system_clock> startTime = std::chrono::system_clock::now();
  const double ar       = 16/9.0; // Aspect ratio
  const int    hdLevel  = 4;      // 1=FHD/2 2=FHD, 4=4k, 8=8k
  mjr::ramCanvas3c8b theRamCanvas(960*hdLevel, 540*hdLevel, -2.2*ar, 2.2*ar, -2.2, 2.2);
  ct_t aColor;

  for(int y=0;y<theRamCanvas.getNumPixY();y++)  {
    if (0==(y % (hdLevel*20)))
        std::cout << "LINE: " << y << " of " << (540*hdLevel) << std::endl;
    for(int x=0;x<theRamCanvas.getNumPixX();x++) {
      cplx fz = f(cplx(theRamCanvas.int2realX(x), theRamCanvas.int2realY(y)));

      // // Schemes designed for coloring complex functions
      // aColor.csSet<ct_t::cs2dRichardson< 10.0, 10.0, 10.0, 1>>(fz);
      // aColor.csSet<ct_t::cs2dThallerHSL< 10.0, 10.0, 10.0, 1>>(fz);
      // aColor.csSet<ct_t::cs2dThallerHSV< 10.0, 10.0, 10.0, 1>>(fz);
      // aColor.csSet<ct_t::cs2dThallerHSVm<10.0, 10.0, 10.0, 1>>(fz);

      // A common choice for color scheme is ct_t::csCColdeRainbow:
      //aColor.csSet<ct_t::cs2dIdxPalArg<ct_t::csCColdeRainbow, 3, 10.0, 10.0, 2.0, 1>>(fz); 

      aColor.csSet<ct_t::cs2dIdxPalArg<ct_t::csCColdeRainbow, 3, 5.0, 20.0, 2.0, 1>>(fz); 

      theRamCanvas.drawPoint(x, y, aColor);
    }
  }
  theRamCanvas.writeTIFFfile("complex_color_image.tiff");
  std::chrono::duration<double> runTime = std::chrono::system_clock::now() - startTime;
  std::cout << "Total Runtime " << runTime.count() << " sec" << std::endl;
}
/** @endcond */
