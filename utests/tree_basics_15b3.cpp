// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      tree_basics_15b3.cpp
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
TEST(tree_basics_15b3, NA) {

  mjr::tree15b3d1rT tree({-1.0, -1.0, -1.0},
                         { 1.0,  1.0,  1.0});

  EXPECT_EQ(tree.ccc_get_top_cell(),                              0x400040004000 );
  EXPECT_EQ(tree.ccc_cell_level(0x400040004000),                  0              );

  EXPECT_EQ(tree.ccc_cell_quarter_width(0x400040004000),          0x2000         ); 
  EXPECT_EQ(tree.ccc_cell_half_width(0x400040004000),             0x4000         ); 
  EXPECT_EQ(tree.ccc_cell_full_width(0x400040004000),             0x8000         ); 

  EXPECT_EQ(tree.ccc_cell_get_corner_min(0x400040004000),         0x000000000000 );
  EXPECT_EQ(tree.ccc_cell_get_corner_max(0x400040004000),         0x800080008000 );

  for(int i=0; i<2; i++) {
    EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x000000000000), i), -1.0,       0.00001);  
    EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x200020002000), i), -0.5,       0.00001);
    EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x400040004000), i),  0.0,       0.00001);
    EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x800080008000), i),  1.0,       0.00001);  

    EXPECT_NEAR(tree.dom_at(tree.get_bbox_min(), i), -1.0,       0.00001);
    EXPECT_NEAR(tree.dom_at(tree.get_bbox_max(), i),  1.0,       0.00001);
    EXPECT_NEAR(tree.dom_at(tree.get_bbox_delta(), i),  1.0/16384, 0.00001);
  }

  EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x800040002000), 0), -0.5,       0.00001);
  EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x800040002000), 1),  0.0,       0.00001);
  EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x800040002000), 2),  1.0,       0.00001);

  EXPECT_EQ(tree.cuc_get_crd(0xCCC1BBB1AAA1, 0),                  0xAAA1         );
  EXPECT_EQ(tree.cuc_get_crd(0xCCC1BBB1AAA1, 1),                  0xBBB1         );
  EXPECT_EQ(tree.cuc_get_crd(0xCCC1BBB1AAA1, 2),                  0xCCC1         );

  EXPECT_EQ(tree.cuc_inc_crd(0xCCC1BBB1AAA1, 0, 0x1),             0xCCC1BBB1AAA2 );
  EXPECT_EQ(tree.cuc_dec_crd(0xCCC1BBB1AAA1, 0, 0x1),             0xCCC1BBB1AAA0 );

  EXPECT_EQ(tree.cuc_inc_crd(0xCCC1BBB1AAA1, 1, 0x1),             0xCCC1BBB2AAA1 );
  EXPECT_EQ(tree.cuc_dec_crd(0xCCC1BBB1AAA1, 1, 0x1),             0xCCC1BBB0AAA1 );

  EXPECT_EQ(tree.cuc_inc_crd(0xCCC1BBB1AAA1, 2, 0x1),             0xCCC2BBB1AAA1 );
  EXPECT_EQ(tree.cuc_dec_crd(0xCCC1BBB1AAA1, 2, 0x1),             0xCCC0BBB1AAA1 );

  EXPECT_EQ(tree.cuc_dec_all_crd(0xCCC1BBB1AAA1, 0x1),            0xCCC0BBB0AAA0 );
  EXPECT_EQ(tree.cuc_inc_all_crd(0xCCC1BBB1AAA1, 0x1),            0xCCC2BBB2AAA2 );

  EXPECT_EQ(tree.cuc_set_all_crd(0xAAA1),                         0xAAA1AAA1AAA1 );

}
