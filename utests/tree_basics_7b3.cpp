// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      tree_basics_7b3.cpp
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
TEST(tree_basics_7b3, NA) {

  mjr::tree7b3d1rT tree({-1.0, -1.0, -1.0},
                        { 1.0,  1.0,  1.0});

  EXPECT_EQ(tree.ccc_get_top_cell(),                        0x404040 );
  EXPECT_EQ(tree.ccc_cell_level(0x404040),                  0        );

  EXPECT_EQ(tree.ccc_cell_quarter_width(0x404040),          0x20     ); 
  EXPECT_EQ(tree.ccc_cell_half_width(0x404040),             0x40     ); 
  EXPECT_EQ(tree.ccc_cell_full_width(0x404040),             0x80     ); 

  EXPECT_EQ(tree.ccc_cell_get_corner_min(0x404040),         0x000000 );
  EXPECT_EQ(tree.ccc_cell_get_corner_max(0x404040),         0x808080 );

  for(int i=0; i<2; i++) {
    EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x000000), i), -1.0,    0.00001);  
    EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x202020), i), -0.5,    0.00001);
    EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x404040), i),  0.0,    0.00001);
    EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x808080), i),  1.0,    0.00001);  

    EXPECT_NEAR(tree.dom_at(tree.get_bbox_min(), i), -1.0,    0.00001);
    EXPECT_NEAR(tree.dom_at(tree.get_bbox_max(), i),  1.0,    0.00001);
    EXPECT_NEAR(tree.dom_at(tree.get_bbox_delta(), i),  1.0/64, 0.00001);
  }

  EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x804020), 0), -0.5,    0.00001);
  EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x804020), 1),  0.0,    0.00001);
  EXPECT_NEAR(tree.dom_at(tree.diti_to_drpt(0x804020), 2),  1.0,    0.00001);

  EXPECT_EQ(tree.cuc_get_crd(0xC1B1A1, 0),                  0xA1     );
  EXPECT_EQ(tree.cuc_get_crd(0xC1B1A1, 1),                  0xB1     );
  EXPECT_EQ(tree.cuc_get_crd(0xC1B1A1, 2),                  0xC1     );

  EXPECT_EQ(tree.cuc_inc_crd(0xC1B1A1, 0, 0x1),             0xC1B1A2 );
  EXPECT_EQ(tree.cuc_dec_crd(0xC1B1A1, 0, 0x1),             0xC1B1A0 );

  EXPECT_EQ(tree.cuc_inc_crd(0xC1B1A1, 1, 0x1),             0xC1B2A1 );
  EXPECT_EQ(tree.cuc_dec_crd(0xC1B1A1, 1, 0x1),             0xC1B0A1 );

  EXPECT_EQ(tree.cuc_inc_crd(0xC1B1A1, 2, 0x1),             0xC2B1A1 );
  EXPECT_EQ(tree.cuc_dec_crd(0xC1B1A1, 2, 0x1),             0xC0B1A1 );

  EXPECT_EQ(tree.cuc_dec_all_crd(0xC1B1A1, 0x1),            0xC0B0A0 );
  EXPECT_EQ(tree.cuc_inc_all_crd(0xC1B1A1, 0x1),            0xC2B2A2 );

  EXPECT_EQ(tree.cuc_set_all_crd(0xA1),                     0xA1A1A1 );

}
