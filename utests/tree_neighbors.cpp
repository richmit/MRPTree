// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      tree_neighbors.cpp
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main
#include <boost/test/unit_test.hpp>

#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE Main
#endif
#include <boost/test/unit_test.hpp>

#include "MR_rect_tree.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BOOST_AUTO_TEST_CASE(tree_neighbors) {
// What we are testing:
//   - Number of neighbors returned
//   - Correct neighbors returned

  mjr::tree7b1d1rT tree1;

  BOOST_TEST_CHECK(tree1.ccc_get_neighbors(0x40).size() == 0    );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbor( 0x40, 0, -1) == 0    );  // 0 means not a neighbor
  BOOST_TEST_CHECK(tree1.ccc_get_neighbor( 0x40, 0,  1) == 0    );

  BOOST_TEST_CHECK(tree1.ccc_get_neighbors(0x20).size() == 1    );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbors(0x20)[0]     == 0x60 );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbor( 0x20, 0, -1) == 0    );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbor( 0x20, 0,  1) == 0x60 );

  BOOST_TEST_CHECK(tree1.ccc_get_neighbors(0x60).size() == 1    );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbors(0x60)[0]     == 0x20 );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbor( 0x60, 0, -1) == 0x20 );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbor( 0x60, 0,  1) == 0    );

  BOOST_TEST_CHECK(tree1.ccc_get_neighbors(0x30).size() == 2    );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbors(0x30)[0]     == 0x10 );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbors(0x30)[1]     == 0x50 );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbor( 0x30, 0, -1) == 0x10 );
  BOOST_TEST_CHECK(tree1.ccc_get_neighbor( 0x30, 0,  1) == 0x50 );

  mjr::tree7b2d1rT tree2;

  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x4040).size() == 0      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x4040, 0, -1) == 0      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x4040, 0,  1) == 0      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x4040, 1, -1) == 0      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x4040, 1,  1) == 0      );

  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x2020).size() == 2      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x2020)[0]     == 0x2060 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x2020)[1]     == 0x6020 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x2020, 0, -1) == 0      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x2020, 0,  1) == 0x2060 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x2020, 1, -1) == 0      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x2020, 1,  1) == 0x6020 );

  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x6060).size() == 2      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x6060)[0]     == 0x6020 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x6060)[1]     == 0x2060 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x6060, 0, -1) == 0x6020 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x6060, 0,  1) == 0      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x6060, 1, -1) == 0x2060 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x6060, 1,  1) == 0      );

  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x1030).size() == 3      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x1030)[0]     == 0x1010 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x1030)[1]     == 0x1050 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x1030)[2]     == 0x3030 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x1030, 0, -1) == 0x1010 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x1030, 0,  1) == 0x1050 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x1030, 1, -1) == 0      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x1030, 1,  1) == 0x3030 );

  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x3030).size() == 4      );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x3030)[0]     == 0x3010 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x3030)[1]     == 0x3050 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x3030)[2]     == 0x1030 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbors(0x3030)[3]     == 0x5030 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x3030, 0, -1) == 0x3010 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x3030, 0,  1) == 0x3050 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x3030, 1, -1) == 0x1030 );
  BOOST_TEST_CHECK(tree2.ccc_get_neighbor( 0x3030, 1,  1) == 0x5030 );

  mjr::tree7b3d1rT tree3;

  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x404040).size() == 0        );

  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x202020).size() == 3        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x202020)[0]     == 0x202060 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x202020)[1]     == 0x206020 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x202020)[2]     == 0x602020 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x202020, 0, -1) == 0        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x202020, 0,  1) == 0x202060 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x202020, 1, -1) == 0        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x202020, 1,  1) == 0x206020 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x202020, 2, -1) == 0        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x202020, 2,  1) == 0x602020 );

  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x606060).size() == 3        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x606060)[0]     == 0x606020 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x606060)[1]     == 0x602060 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x606060)[2]     == 0x206060 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x606060, 0, -1) == 0x606020 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x606060, 0,  1) == 0        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x606060, 1, -1) == 0x602060 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x606060, 1,  1) == 0        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x606060, 2, -1) == 0x206060 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x606060, 2,  1) == 0        );

  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x101030).size() == 4        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x101030)[0]     == 0x101010 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x101030)[1]     == 0x101050 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x101030)[2]     == 0x103030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x101030)[3]     == 0x301030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x101030, 0, -1) == 0x101010 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x101030, 0,  1) == 0x101050 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x101030, 1, -1) == 0        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x101030, 1,  1) == 0x103030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x101030, 2, -1) == 0        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x101030, 2,  1) == 0x301030 );

  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x103030).size() == 5        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x103030)[0]     == 0x103010 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x103030)[1]     == 0x103050 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x103030)[2]     == 0x101030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x103030)[3]     == 0x105030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x103030)[4]     == 0x303030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x103030, 0, -1) == 0x103010 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x103030, 0,  1) == 0x103050 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x103030, 1, -1) == 0x101030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x103030, 1,  1) == 0x105030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x103030, 2, -1) == 0        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x103030, 2,  1) == 0x303030 );

  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x303030).size() == 6        );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x303030)[0]     == 0x303010 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x303030)[1]     == 0x303050 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x303030)[2]     == 0x301030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x303030)[3]     == 0x305030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x303030)[4]     == 0x103030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbors(0x303030)[5]     == 0x503030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x303030, 0, -1) == 0x303010 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x303030, 0,  1) == 0x303050 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x303030, 1, -1) == 0x301030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x303030, 1,  1) == 0x305030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x303030, 2, -1) == 0x103030 );
  BOOST_TEST_CHECK(tree3.ccc_get_neighbor( 0x303030, 2,  1) == 0x503030 );

  mjr::tree7b4d1rT tree4;

  BOOST_TEST_CHECK(tree4.ccc_get_neighbors(0x40404040).size() == 0 );
}