// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      MR_math.hpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-09-12
 @brief     Very simple math stuff very carefully implimented.@EOL
 @keywords
 @std       C++20
 @see
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

/*############################################################################################################################################################*/
#ifndef MJR_INCLUDE_mrmath

namespace mjr {
  namespace math {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Check that x is between r1 and r2, but not equal to them.
    @param x       Value to test.  Must be a floating point type.
    @param r1      Lower interval limit.  Must be the same type as x.
    @param r2      Upper interval limit.  Must be the same type as x.
    @param epsilon Epsilion.  Must be the saem type as x.
    @return Return true if @f$x\in (r_1+\epsilon, r_2-\epsilon)@f$, and false otherwise. */
    template <typename numType>
    requires (std::floating_point<numType>)
    inline int
    fbetween(numType x, numType r1, numType r2, numType epsilon) {
      return ( (x > (r1 + epsilon)) && (x < (r2 - epsilon)) );
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Detect if a floating point value is near zero.
    @param x            Value to test.  Must be a floating point type.
    @param zero_epsilon Epsilion to detect zero sign.  Must be the saem type as x.
    @return Returns true if x in the closed ball of radius zero_epsilon around zero, otherwise returns false. */
    template <typename numType>
    requires (std::floating_point<numType>)
    inline bool fnear_zero(numType x, numType zero_epsilon) {
      if ((x <= zero_epsilon) && (x >= -zero_epsilon)) { // logically the same as: (std::abs(x) <= zero_epsilon)
        return true;
      } else {
        return false;
      }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Detect if a floating point value is not near zero.
    @param x            Value to test.  Must be a floating point type.
    @param zero_epsilon Epsilion to detect zero sign.  Must be the saem type as x.
    @return Returns true if x is outside the closed ball of radius zero_epsilon around zero, otherwise returns false. */
    template <typename numType>
    requires (std::floating_point<numType>)
    inline bool fnot_near_zero(numType x, numType zero_epsilon) {
      if ((x > zero_epsilon) || (x < -zero_epsilon)) { // logically the same as: (std::abs(x) > zero_epsilon)
        return true;
      } else {
        return false;
      }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Detect if two floating point values are near each other.
    @param x1           Value to test.  Must be a floating point type.
    @param x2           Value to test.  Must be an integer or floating point type.
    @param zero_epsilon Epsilion to detect zero sign.  Must be the saem type as x.
    @return Returns true if @f$\vert x_1-x_2\vert\le\epsilon@f$, otherwise returns false. */
    template <typename numType>
    requires (std::floating_point<numType>)
    inline bool fnear(numType x1, numType x2, numType zero_epsilon) {
      if (std::abs(x1-x2) <= zero_epsilon) {
        return true;
      } else {
        return false;
      }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Detect if two floating point values are not near each other.
    @param x1           Value to test.  Must be a floating point type.
    @param x2           Value to test.  Must be an integer or floating point type.
    @param zero_epsilon Epsilion to detect zero sign.  Must be the saem type as x.
    @return Returns true if @f$\vert x_1-x_2\vert\gt\epsilon@f$, otherwise returns false. */
    template <typename numType>
    requires (std::floating_point<numType>)
    inline bool fnot_near(numType x1, numType x2, numType zero_epsilon) {
      if (std::abs(x1-x2) > zero_epsilon) {
        return true;
      } else {
        return false;
      }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Signum (sgn, sign) function -- -1 if input is negative, +1 if it is positive, and 0 if it is zero.
    @param x        Value to test.  Must be an integer or floating point type.
    @return The sign */
    template <typename numType>
    requires (std::integral<numType> || std::floating_point<numType>)
    inline int
    sgn(numType x) {
      if (x > numType(0)) return  1;
      if (x < numType(0)) return -1;
      return 0;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Signum with zero epsilon check.
    This function is used if we want to make sure a non-zero sign is only returned if the input value is larger than some epsilon.  i.e. we don't want tiny
    values near zero to show up as positive or negative.
    @param x            Value to test.  Must be an integer or floating point type.
    @param zero_epsilon Epsilion to detect zero sign.  Must be the saem type as x.
    @return 0 if x in the closed ball of radius zero_epsilon around zero, otherwise like sgn(x). */
    template <typename numType>
    requires (std::floating_point<numType>)
    inline int
    sgne(numType x, numType zero_epsilon) {
      if (fnear_zero(x, zero_epsilon)) {
        return 0;
      } else {
        if (x > numType(0))
          return 1;
        else
          return -1;
      }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Minimum of three numbers.
    @param x1 First number.Must be an integer or floating point type.
    @param x2 Second number.  Must be the saem type as x1.
    @param x3 Third number.  Must be the saem type as x1.
    @return Minimum of input values. */
    template <typename numType>
    inline numType min3(numType x1, numType x2, numType x3) {
      if (x1 < x2) {   // opts: x1 or x3  (x1<x2)
        if (x1 < x3) {    // opts: X1        (x3>x1<x2)
          return x1;
        } else {          // opts: X3        (x3<=x1<x2)
          return x3;
        }
      } else {         // opts: X2 or X3  (x1>=x2)
        if (x2 < x3) {    // opts: X2        (x1>=x2<x3)
          return x2;
        } else {          // opts: X3        (x1>=x2>=x3)
          return x3;
        }
      }
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** Minimum of three numbers.
    @param x1 First number.Must be an integer or floating point type.
    @param x2 Second number.  Must be the saem type as x1.
    @param x3 Third number.  Must be the saem type as x1.
    @return Minimum of input values. */
    template <typename numType>
    inline numType max3(numType x1, numType x2, numType x3) {
      if (x2 < x1) {   // opts: x1 or x3
        if (x3 < x1) {    // opts: X1
          return x1;
        } else {          // opts: X3
          return x3;
        }
      } else {         // opts: X2 or X3
        if (x3 < x2) {    // opts: X2
          return x2;
        } else {          // opts: X3
          return x3;
        }
      }
    }

  }
}

#define MJR_INCLUDE_mrmath
#endif
