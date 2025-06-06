// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      check_cell_segment.cpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Unit tests for MR_cell_cplx.@EOL
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
#include "MR_cell_cplx.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(FUN_check_cell, segment) {

  mjr::MRccT5 aPoly;

  aPoly.add_node({0.0, 0.0, 0.0});
  aPoly.add_node({1.0, 1.0, 1.0});
  aPoly.add_node({1.0, 2.0, 3.0});
  aPoly.add_node({2.0, 2.0, 2.0});

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------

  EXPECT_EQ(aPoly.check_cell_vertexes(mjr::MRccT5::cell_kind_t::SEGMENT, mjr::MRccT5::cell_verts_t({0,-1})),       mjr::MRccT5::cell_stat_t::NEG_PNT_IDX);

  EXPECT_EQ(aPoly.check_cell_vertexes(mjr::MRccT5::cell_kind_t::SEGMENT, mjr::MRccT5::cell_verts_t({0,10})),       mjr::MRccT5::cell_stat_t::BIG_PNT_IDX);

  EXPECT_EQ(aPoly.check_cell_vertexes(mjr::MRccT5::cell_kind_t::SEGMENT, mjr::MRccT5::cell_verts_t({0,1,2})),      mjr::MRccT5::cell_stat_t::TOO_MANY_PNT);

  EXPECT_EQ(aPoly.check_cell_vertexes(mjr::MRccT5::cell_kind_t::SEGMENT, mjr::MRccT5::cell_verts_t({0})),          mjr::MRccT5::cell_stat_t::TOO_FEW_PNT);

  EXPECT_EQ(aPoly.check_cell_vertexes(mjr::MRccT5::cell_kind_t::SEGMENT, mjr::MRccT5::cell_verts_t({0, 0})),       mjr::MRccT5::cell_stat_t::DUP_PNT);

  EXPECT_EQ(aPoly.check_cell_vertexes(mjr::MRccT5::cell_kind_t::SEGMENT, mjr::MRccT5::cell_verts_t({0,1})),        mjr::MRccT5::cell_stat_t::GOOD);

  //--------------------------------------------------------------------------------------------------------------------------------------------------------------

  EXPECT_EQ(aPoly.check_cell_dimension(mjr::MRccT5::cell_kind_t::SEGMENT, mjr::MRccT5::cell_verts_t({0, 0})),      mjr::MRccT5::cell_stat_t::GOOD); // Always returns GOOD for SEGMENTs.

  EXPECT_EQ(aPoly.check_cell_dimension(mjr::MRccT5::cell_kind_t::SEGMENT, mjr::MRccT5::cell_verts_t({0,1})),       mjr::MRccT5::cell_stat_t::GOOD);

}
