// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      tree_corners.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-14
 @brief     Unit tests for MR_rect_tree.@EOL
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
*/
/*******************************************************************************************************************************************************.H.E.**/

#include <gtest/gtest.h>
#include "MR_rect_tree.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(tree_corners, NA) {
// What we are testing:
//   - Number of corners returned
//   - Correct corners returned
//   - Order of corners on the list

  mjr::tree7b1d1rT tree1;

  EXPECT_EQ(tree1.ccc_get_corners(0x40).size(),         2 );
  EXPECT_EQ(tree1.ccc_get_corners(0x40)[0],             0x00 );
  EXPECT_EQ(tree1.ccc_get_corners(0x40)[1],             0x80 );

  EXPECT_LT(tree1.ccc_get_corners(0x40)[0], tree1.ccc_get_corners(0x40)[1]);

  EXPECT_EQ(tree1.ccc_get_corners(0x40, 0, -1).size(),  1 );
  EXPECT_EQ(tree1.ccc_get_corners(0x40, 0, -1)[0],      0x00 );

  EXPECT_EQ(tree1.ccc_get_corners(0x40, 0,  1).size(),  1 );
  EXPECT_EQ(tree1.ccc_get_corners(0x40, 0,  1)[0],      0x80 );

  mjr::tree7b2d1rT tree2;

  EXPECT_EQ(tree2.ccc_get_corners(0x4040).size(),         4 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040)[0],             0x0000 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040)[1],             0x0080 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040)[2],             0x8000 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040)[3],             0x8080 );

  for(int i=0; i<3; ++i) 
    EXPECT_LT(tree2.ccc_get_corners(0x4040)[i], tree2.ccc_get_corners(0x4040)[i+1]);

  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 0, -1).size(),  2 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 0, -1)[0],      0x0000 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 0, -1)[1],      0x8000 );

  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 0,  1).size(),  2 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 0,  1)[0],      0x0080 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 0,  1)[1],      0x8080 );

  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 1, -1).size(),  2 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 1, -1)[0],      0x0000 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 1, -1)[1],      0x0080 );

  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 1,  1).size(),  2 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 1,  1)[0],      0x8000 );
  EXPECT_EQ(tree2.ccc_get_corners(0x4040, 1,  1)[1],      0x8080 );

  mjr::tree7b3d1rT tree3;

  EXPECT_EQ(tree3.ccc_get_corners(0x404040).size(),         8 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040)[0],             0x000000 ); // 0-1 1-1 2-1  AXIS+DIR
  EXPECT_EQ(tree3.ccc_get_corners(0x404040)[1],             0x000080 ); // 0+1 1-1 2-1
  EXPECT_EQ(tree3.ccc_get_corners(0x404040)[2],             0x008000 ); // 0-1 1+1 2-1
  EXPECT_EQ(tree3.ccc_get_corners(0x404040)[3],             0x008080 ); // 0+1 1+1 2-1
  EXPECT_EQ(tree3.ccc_get_corners(0x404040)[4],             0x800000 ); // 0-1 1-1 2+1
  EXPECT_EQ(tree3.ccc_get_corners(0x404040)[5],             0x800080 ); // 0+1 1-1 2+1
  EXPECT_EQ(tree3.ccc_get_corners(0x404040)[6],             0x808000 ); // 0-1 1+1 2+1
  EXPECT_EQ(tree3.ccc_get_corners(0x404040)[7],             0x808080 ); // 0+1 1+1 2+1

  for(int i=0; i<7; ++i) 
    EXPECT_LT(tree3.ccc_get_corners(0x404040)[i], tree3.ccc_get_corners(0x404040)[i+1]);

  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0, -1).size(),  4 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0, -1)[0],      0x000000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0, -1)[1],      0x008000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0, -1)[2],      0x800000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0, -1)[3],      0x808000 );

  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0,  1).size(),  4 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0,  1)[0],      0x000080 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0,  1)[1],      0x008080 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0,  1)[2],      0x800080 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 0,  1)[3],      0x808080 );

  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1, -1).size(),  4 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1, -1)[0],      0x000000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1, -1)[1],      0x000080 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1, -1)[2],      0x800000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1, -1)[3],      0x800080 );
                                                                          
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1,  1).size(),  4 );  
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1,  1)[0],      0x008000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1,  1)[1],      0x008080 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1,  1)[2],      0x808000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 1,  1)[3],      0x808080 );

  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2, -1).size(),  4 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2, -1)[0],      0x000000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2, -1)[1],      0x000080 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2, -1)[2],      0x008000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2, -1)[3],      0x008080 );
                                                                          
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2,  1).size(),  4 );  
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2,  1)[0],      0x800000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2,  1)[1],      0x800080 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2,  1)[2],      0x808000 );
  EXPECT_EQ(tree3.ccc_get_corners(0x404040, 2,  1)[3],      0x808080 );

  mjr::tree7b4d1rT tree4;

  EXPECT_EQ(tree4.ccc_get_corners(0x40404040).size(),  16 );

  for(int i=0; i<15; ++i) 
    EXPECT_LT(tree4.ccc_get_corners(0x40404040)[i], tree4.ccc_get_corners(0x40404040)[i+1]);

  for(int dir=-1; dir<2; dir+=2) 
    for(int dim=0; dim<4; ++dim) 
      EXPECT_EQ(tree4.ccc_get_corners(0x40404040, dim, dir).size(),  8 );

  for(int dir=-1; dir<2; dir+=2) 
    for(int dim=0; dim<4; ++dim) 
      for(int i=0; i<7; ++i) 
        EXPECT_LT(tree4.ccc_get_corners(0x40404040, dim, dir)[i], tree4.ccc_get_corners(0x40404040, dim, dir)[i+1]);

}
