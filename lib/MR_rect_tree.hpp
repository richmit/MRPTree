// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      MR_rect_tree.hpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-06-13
 @brief     Implimentation of the MR_rect_tree class.@EOL
 @keywords  quadtree octree lod
 @std       C++23
 @see       MR_rt_to_cc.hpp, MR_cell_cplx.hpp
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef MJR_INCLUDE_MR_rect_tree

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <algorithm>                                                     /* STL algorithm           C++11    */
#include <array>                                                         /* array template          C++11    */
#include <bit>                                                           /* STL bit manipulation    C++20    */
#include <climits>                                                       /* std:: C limits.h        C++11    */
#include <cmath>                                                         /* std:: C math.h          C++11    */
#include <complex>                                                       /* Complex Numbers         C++11    */
#include <concepts>                                                      /* Concepts library        C++20    */
#include <cstdint>                                                       /* std:: C stdint.h        C++11    */
#include <cstring>                                                       /* std:: C string.h        C++11    */
#include <fstream>                                                       /* C++ fstream             C++98    */
#include <functional>                                                    /* STL funcs               C++98    */
#include <iomanip>                                                       /* C++ stream formatting   C++11    */
#include <iostream>                                                      /* C++ iostream            C++11    */
#include <limits>                                                        /* C++ Numeric limits      C++11    */
#include <set>                                                           /* STL set                 C++98    */
#include <span>                                                          /* STL spans               C++20    */
#include <sstream>                                                       /* C++ string stream       C++      */
#include <string>                                                        /* C++ strings             C++11    */
#include <tuple>                                                         /* STL tuples              C++11    */
#include <type_traits>                                                   /* C++ metaprogramming     C++11    */
#include <unordered_map>                                                 /* STL hash map            C++11    */
#include <map>                                                           /* STL map                 C++11    */
#include <utility>                                                       /* STL Misc Utilities      C++11    */
#include <vector>                                                        /* STL vector              C++11    */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MRMathCPP.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put everything in the mjr namespace
namespace mjr {
/** @brief Template Class used to house an MR_rect_tree.

    Overview
    ========

    This class forms the primary component of the MR @f$2^P\mathrm{-Tree}@f$ library (also known as "MRPTree").  In terms of this library, a quadtree is a
    @f$2^2\mathrm{-Tree}@f$ and an octree is a @f$2^3\mathrm{-Tree}@f$ -- in essence this class implements a universal data structure capable of holding such
    trees of any dimension.  Contrary to the name, this is all achieved without using any tree-like data structures at all!

    The `dom_dim` template parameter corresponds to the "@f$P@f$" in the exponent.  If one wishes to implement a quadtree, then dom_dim should be set to 2.
    For octrees the value of dom_dim should be set to 3.  We are not limited to dimensions 2 & 3 -- for example, we could make 1D "bitrees" or a 4D "space
    time trees".

    This class makes the unusual assumption that the sampled data stored in the tree is a finite, real vector space.  So instead of using a template parameter
    to specify the type of the sampled data, we provide 1) a type for the vector space scalar (spc_real_t) and 2) a dimension (rng_dim).  In other words these
    trees are optimized to model model functions @f$f:\mathbb{R}^d\rightarrow\mathbb{R}^r@f$ where @f$d@f$ and @f$r@f$ *are* template parameters(`dom_dim` &
    `rng_dim` respectively).

    The `max_level` parameter will seem odd for anyone who has used a traditional quadtree/octree library.  As a practical matter trees in traditional
    libraries are limited to a shallow depths (say 10 or 12 levels); however, the code bases themselves don't usually place an explicit limit on depth.  In
    this library the maximum depth must be specified from the start.

    As an example, suppose we have a function @f$f:\mathbb{R}^2\rightarrow\mathbb{R}^3@f$ we wish to sample.  For this application we might use the following values:
      - `spc_real_t = double`
      - `dom_dim = 2`
      - `rng_dim = 3`
      - `max_level = 15`

    We could use an alternate value for `max_level`, but I think 15 is a good balance.  We can do quadtrees at this depth using 32-bit integer keys, and we
    can do octrees with 64-bit integer keys.

    Naming Conventions
    ==================

    - Methods for numeric computation on coordinate tuples (::diti_t types).  These methods have nothing to do with cell sample status -- just coordinates!
      - ccc_ -- Cell Coordinate Computation:  Coordinates are interpreted as cell coordinate centers
      - cuc_ -- Coordinate Utility Computation: Coordinate are not necessarily cell coordinate centers
    - Type Naming Conventions
      - Coordinates
        - ::dic_t (domain Integer Coordinate) -- A single integer
        - ::src_t (domain/range space Real Coordinate) -- Floating point
      - Domain Integer Coordinate Tuples
        - ::diti_t (domain int tuple Integer) -- Encoded as a packed integer
        - ::dita_t (domain int tuple Array)   -- As a dom_dim element std::array
        - ::ditv_t (domain int tuple Vector)  -- As a dom_dim element std::vector
      - Domain/range real Coordinate  tuples
        - ::rrta_t & ::drta_t (range/domain real tuple Array)  -- As a `dom_dim` or `rng_dim` element std::array
        - ::rrtv_t & ::drtv_t (range/domain real tuple Vector) -- As a `dom_dim` or `rng_dim` element std::vector
        - ::rrpt_t & ::drpt_t (range/domain real psudo-tuple)  -- As a ::src_t when `dim==1`, and an array otherwise
    Function arguments of type "`foo_t`" are frequently called "`foo`".

    Details
    =======

    @tparam max_level  The maximum depth of the tree -- use one minus a power of two for highest performance.
    @tparam spc_real_t The base floating type to use for both domain & range.
    @tparam dom_dim    Domain dimension.
    @tparam rng_dim    Range dimension. */
  template <int max_level, class spc_real_t, int dom_dim, int rng_dim>
  requires ((max_level>0)                                  &&
            (dom_dim>0)                                    &&
            (rng_dim>0)                                    &&
            (dom_dim*max_level<=CHAR_BIT*sizeof(uint64_t)) &&
            (std::is_floating_point<spc_real_t>::value))
  class MR_rect_tree {
    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Template Parameter Constants */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      constexpr static int domain_dimension = dom_dim;   //!< The value of the template parameter dom_dim
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      constexpr static int range_dimension  = rng_dim;   //!< The value of the template parameter rng_dim
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      constexpr static int maximum_level    = max_level; //!< The value of the template parameter max_level
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Template Parameter Types */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Externally exposed typedef for spc_real_t */
      typedef MR_rect_tree<max_level, spc_real_t, dom_dim, rng_dim> this_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Externally exposed typedef for spc_real_t */
      typedef spc_real_t src_t;
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Domain Real Coordinates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** std::vector for values in the domain space. */
      typedef std::vector<src_t> drtv_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** An std::array for points in the domain space. */
      typedef std::array<src_t, dom_dim> drta_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** For values in the domain space.
            - When dom_dim==1, this will be src_t
            - When dom_dim!=1, this will be drta_t (an std::array) */
      typedef typename std::conditional<std::cmp_equal(dom_dim, 1), src_t, drta_t>::type drpt_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** A nicely descriptive typedef for drpt_t. */
      typedef drpt_t real_domain_t;
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Range Real Coordinates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** std::vector for values in the range space. */
      typedef std::vector<src_t> rrtv_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** An std::array for points in the range space. */
      typedef std::array<src_t, rng_dim> rrta_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** For values in the range space.
            - When rng_dim==1, this will be src_t
            - When rng_dim!=1, this will be a rrta_t (an std::array) */
      typedef typename std::conditional<std::cmp_equal(rng_dim, 1), src_t, rrta_t>::type rrpt_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** A nicely descriptive typedef for rrpt_t */
      typedef rrpt_t real_range_t;
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Domain Integer Coordinates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** The number of bits used by a single component of an integer coordinate tuple.
          @warning dic_bits >= sizeof(dic_t), but it might not be equal. */
      constexpr static int dic_bits = max_level+1;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** The number of bits used by an entire integer coordinate tuple.
          @warning diti_bits >= sizeof(diti_t), but it might not be equal. */
      constexpr static int diti_bits = dic_bits * dom_dim;
      //@}

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Domain Integer Coordinates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /* Private to keep it undocumented in doxygen because the typedef definition is too large */
      typedef typename std::conditional<std::cmp_greater_equal(CHAR_BIT*sizeof(int8_t),  dic_bits), uint8_t,
              typename std::conditional<std::cmp_greater_equal(CHAR_BIT*sizeof(int16_t), dic_bits), uint16_t,
              typename std::conditional<std::cmp_greater_equal(CHAR_BIT*sizeof(int32_t), dic_bits), uint32_t,
                                                                                                    uint64_t
                                        >::type>::type>::type priv_dic_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /* Private to keep it undocumented in doxygen because the typedef definition is too large */
      typedef typename std::conditional<std::cmp_greater_equal(CHAR_BIT*sizeof(int8_t),  diti_bits), uint8_t,
              typename std::conditional<std::cmp_greater_equal(CHAR_BIT*sizeof(int16_t), diti_bits), uint16_t,
              typename std::conditional<std::cmp_greater_equal(CHAR_BIT*sizeof(int32_t), diti_bits), uint32_t,
                                                                                                     uint64_t
                                        >::type>::type>::type priv_diti_t;
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Domain Integer Coordinates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Unsigned integer type large large enough to hold an integer coordiante component */
      typedef priv_dic_t dic_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Unsigned integer type large large enough to hold an integer coordiante tuple */
      typedef priv_diti_t diti_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Domain point given as an std::array integer tuple. */
      typedef std::array<dic_t, 3> dita_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Domain point given as an std::vector integer tuple. */
      typedef std::vector<dic_t> ditv_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** A std::vector used to pass lists of diti_t types around.  */
      typedef std::vector<diti_t> diti_list_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      constexpr static dic_t dic_max = (static_cast<dic_t>(1) << max_level);     //!< Maximum allowd for a dic_t
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      constexpr static dic_t dic_ctr = (static_cast<dic_t>(1) << (max_level-1)); //!< Center value for a dic_t
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      constexpr static dic_t dic_min = 0;                                        //!< Minimum allowd for a dic_t
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Function Types */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Real input sample function */
      typedef std::function<rrpt_t(drpt_t)> drpt2rrpt_func_t;       // dr2rr (Sample Function)
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Integer input predicate function */
      typedef std::function<bool(diti_t)>   diti2bool_func_t;      //  di2b (Domain Index Predicate)
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Real input predicate function */
      typedef std::function<bool(drpt_t)>   drpt2bool_func_t;      //  dr2b (Domain Point Predicate)
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Real input, single real variable output function */
      typedef std::function<src_t(drpt_t)>  drpt2real_func_t;      //  dr2r (Domain Point SDF)
      //@}

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Domain Integer Coordinates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      constexpr static diti_t diti_ones = ~static_cast<diti_t>(0);                         // diti_t int with all ones
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      constexpr static diti_t diti_msk0 = static_cast<diti_t>(~(diti_ones << dic_bits));   // diti_t int with ones on 0th coord
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Data Members */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      drpt_t bbox_min;      //!< Holds the minimal point for the real domain range
      drpt_t bbox_max;      //!< Holds the maximal point for the real domain range
      drpt_t bbox_delta;    //!< The wdith of the real domain range
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      std::unordered_map<diti_t, rrpt_t> samples; //!< Holds the sampled data
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Constructors */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Set real coordinates to defaults.  see: set_bbox_default(). */
      MR_rect_tree()                                                            { set_bbox_default();                     }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Set real coordinate as specified.
          @param new_bbox_min Value to use for bounding box minimum point
          @param new_bbox_max Value to use for bounding box maximum point */
      MR_rect_tree(drpt_t new_bbox_min, drpt_t new_bbox_max)                    { set_bbox(new_bbox_min, new_bbox_max);   }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Construction Helpers */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Update the value of bbox_delta.
          Use this after modifying the value of bbox_min or bbox_max. */
      void update_bbox_delta() {
        if constexpr (dom_dim == 1) {
          bbox_delta = (bbox_max - bbox_min) / ((dic_t(1) << max_level));
          if (bbox_max <= bbox_min) {
            std::cerr << "ERROR(update_bbox_delta): bbox_min must be less than bbox_max!" << std::endl;
            exit(1);
          }
        } else {
          for(int i=0; i<dom_dim; i++) {
            bbox_delta[i] = (bbox_max[i] - bbox_min[i]) / ((dic_t(1) << max_level));
            if (bbox_max[i] <= bbox_min[i]) {
              std::cerr << "ERROR(update_bbox_delta): Corrisponding elements of bbox_min must be less than bbox_max!" << std::endl;
              exit(1);
            }
          }
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Set the bounding box
          @param new_bbox_min Value to use for bounding box minimum point
          @param new_bbox_max Value to use for bounding box maximum point */
      void set_bbox(drpt_t new_bbox_min, drpt_t new_bbox_max) {
        bbox_min = new_bbox_min;
        bbox_max = new_bbox_max;
        update_bbox_delta();
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Set the bounding box to the default: -1 for all minimum components and +1 for all maximum components. */
      void set_bbox_default() {
        if constexpr (dom_dim == 1) {
          bbox_min = static_cast<src_t>(-1.0);
          bbox_max = static_cast<src_t>( 1.0);
        } else {
          bbox_min.fill(static_cast<src_t>(-1.0));
          bbox_max.fill(static_cast<src_t>( 1.0));
        }
        update_bbox_delta();
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Set the bounding box minimum.
          @param new_bbox_min Value to use for bounding box minimum point*/
      void set_bbox_min(drpt_t new_bbox_min) { set_bbox(new_bbox_min, bbox_max); };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Set the bounding box max_level.
          @param new_bbox_max Value to use for bounding box maximum point*/
      void set_bbox_max(drpt_t new_bbox_max) { set_bbox(bbox_min, new_bbox_max); };
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Basic Class Info */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the bounding box minimum point */
      inline drpt_t get_bbox_min() const     { return (bbox_min); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the bounding box minimum point */
      inline drpt_t get_bbox_max() const     { return (bbox_max); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the bounding box minimum point */
      inline drpt_t get_bbox_delta() const   { return (bbox_delta); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the sample value for vertex.
          @param vertex Input vertex */
      inline rrpt_t get_sample(diti_t vertex) const { return samples.at(vertex); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the sample value for vertex as an rrta_t (an std::array)
          @param vertex Input vertex */
      inline rrta_t get_sample_rrta(diti_t vertex) const {
        rrta_t ret;
        if constexpr (rng_dim == 1)
          ret[0] = get_sample(vertex);
        else
          ret = get_sample(vertex);
        return ret;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Provide a constant forward iterator for the sample data.
          Sample data is stored as a pair with the first element being the packed integer domain coordinates and the second being the sampled data. */
      inline  std::unordered_map<diti_t, rrpt_t>::const_iterator cbegin_samples()  { return samples.cbegin(); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Provide a constant end iterator for the sample data.
          see: cbegin_samples(). */
      inline  std::unordered_map<diti_t, rrpt_t>::const_iterator   cend_samples()  { return samples.cend();   }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cell Oriented Coordinate Computation
          These functions compute theoretical values based on a given cell center coordinate, and have nothing to do with sample state.

          For example, the ccc_get_top_cell() function returns the coordinates that would be used to identify the tree's top cell; however, this function tells us
          nothing about if that cell exists (as been sampled) in the tree.

          In general these routines are optimized for performance, and do not preform any error checking.

            - Many functions require a valid cell coordinate (that is a coordinate that could be the center of a cell).  For example, if one of these
              functions is given 0, then erroneous results are likely because 0 represents the coordinate for the corner of a tree which can never the the
              center of a cell.

            - When possible, the fewest coordinate components are used.  For example ccc_cell_level() uses just the first coordinate in it's computation. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the top cell coordinates for the tree
          @warning Just returns the coordinates regardless of if the top cell exists (is sampled).
          @return Top cell coordinate for the tree. */
      inline diti_t ccc_get_top_cell() const { return cuc_set_all_crd(dic_ctr); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute cell level
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @warning Only uses the last dic in the cell to compute value.  Thus inconsistant components will not be detected.
          @param cell Input cell
          @return The level of the given cell. */
      inline dic_t ccc_cell_level(diti_t cell) const     { return static_cast<dic_t>(max_level-static_cast<diti_t>(1)-std::countr_zero(cell)); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute cell quarter width
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell Input cell
          @return Quarter width of the given cell. */
      inline dic_t ccc_cell_quarter_width(diti_t cell) const { return static_cast<dic_t>(ccc_cell_half_width(cell) >> static_cast<dic_t>(1)); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute cell half width
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell Input cell
          @return Half width of the given cell. */
       inline dic_t ccc_cell_half_width(diti_t cell) const { return static_cast<dic_t>(cell & (~cell+static_cast<diti_t>(1))); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute cell full width
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell Input cell
          @return Full width of the given cell. */
      inline dic_t ccc_cell_full_width(diti_t cell) const { return static_cast<dic_t>(ccc_cell_half_width(cell) << static_cast<dic_t>(1)); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Cell first corner
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell Input cell
          @return First corner of given cell */
      inline diti_t ccc_cell_get_corner_min(diti_t cell) const { return (cuc_dec_all_crd(cell, ccc_cell_half_width(cell))); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Cell last corner
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell Input cell
          @return Last corner of given cell */
      inline diti_t ccc_cell_get_corner_max(diti_t cell) const { return (cuc_inc_all_crd(cell, ccc_cell_half_width(cell))); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return a list of the corners of the given cell
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell Input cell */
      diti_list_t ccc_get_corners(diti_t cell) const { return cuc_two_cross(cell, ccc_cell_half_width(cell)); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return a list of the corners of the given cell
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell Input cell
          @param index     The index of the axis.  Must be in [0, dom_dim-1].  No error checking.
          @param direction The direction on the given index.  Must be 1 or -1.  No error checking. */
      diti_list_t ccc_get_corners(diti_t cell, int index, int direction) const { return cuc_two_cross(cell, ccc_cell_half_width(cell), index, direction); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return a list of potential neighbor cells of the specified cell
          Note the cells are not in canonical order!
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell Input cell */
      diti_list_t ccc_get_neighbors(diti_t cell) const { return cuc_axis_cross(cell, ccc_cell_full_width(cell)); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the potential neighbor cell along the given axis in the specified direction.
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          @param cell      Input cell
          @param index     The index of the axis.  Must be in [0, dom_dim-1].  No error checking.
          @param direction The direction on the given index.  Must be 1 or -1.  No error checking.
          @return A neighbor cell coordinate or 0 if no neighbor exists.  Note 0 is a valid coordinate, but not the center of any cell. */
      diti_t ccc_get_neighbor(diti_t cell, int index, int direction) const {
        dic_t tmp   = cuc_get_crd(cell, index);
        dic_t delta = ccc_cell_full_width(cell);
        if (direction == 1) {
          if ((dic_max-tmp) >= delta)
            return cuc_inc_crd(cell, index, delta);
        } else {
          if (tmp >= delta)
            return cuc_dec_crd(cell, index, delta);
        }
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return a list of child cells of the specified cell
          @warning No error checking -- cell must be a valid center coordinate. See: cell_good_cords()
          If cell can't have children because it is at max_level, then an empty vector is returned.
          @warning This isn't a check for existing, sampled children -- it simply returns the coordinates of potential children.
          @param cell Input cell */
      diti_list_t ccc_get_children(diti_t cell) const {
        if (ccc_cell_level(cell) < (max_level-1)) {
          return cuc_two_cross(cell, ccc_cell_quarter_width(cell));
        } else {
          return diti_list_t();
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return a list of child cells of the specified cell
          An empty vector is returned if no children are possible.
          @warning This isn't a check for existing, sampled children -- it simply returns the coordinates of potential children.
          @param cell Input cell
          @param index     The index of the axis.  Must be in [0, dom_dim-1].  No error checking.
          @param direction The direction on the given index.  Must be 1 or -1.  No error checking.
      */
      diti_list_t ccc_get_children(diti_t cell, int index, int direction) const {
        if (ccc_cell_level(cell) < (max_level-1)) {
          return cuc_two_cross(cell, ccc_cell_quarter_width(cell), index, direction);
        } else {
          return diti_list_t();
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return a list of the vertexes (corners and center) of the given cell
          @param cell Input cell */
      diti_list_t ccc_get_vertexes(diti_t cell) const {
        diti_list_t rv = ccc_get_corners(cell);
        rv.push_back(cell);
        return rv;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Low Level Integer Tuple Computation
          These functions encapsulate various coordinate computations.

          In general these routines are optimized for performance, and do not preform any error checking.  For example it is possible to underflow/overflow
          with inappropriate values. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Extract a component from a tuple
          @param diti The input tuple
          @param index Which component to extract
          @return The component at index position */
      inline dic_t  cuc_get_crd(diti_t diti, int index) const { return static_cast<dic_t>(diti_msk0 & static_cast<diti_t>(diti >> (dic_bits * index))); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Incriment a component in a a tuple
          @param diti The input tuple
          @param index Which component to incriment
          @param value Amout by which to increment the component
          @return New tuple with the specified component incrimented */
      inline diti_t cuc_inc_crd(diti_t diti, int index, dic_t value) const { return (diti + (static_cast<diti_t>(static_cast<diti_t>(value) << (dic_bits*index)))); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Decriment a component in a a tuple
          @param diti The input tuple
          @param index Which component to incriment
          @param value Amout by which to decrement the component
          @return New tuple with the specified component decrimented */
      inline diti_t cuc_dec_crd(diti_t diti, int index, dic_t value) const { return (diti - (static_cast<diti_t>(static_cast<diti_t>(value) << (dic_bits*index)))); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Decriment all components in a a tuple
          @param diti The input tuple
          @param value Amout by which to decrement each component
          @return New tuple with decrimented components */
      inline diti_t cuc_dec_all_crd(diti_t diti, dic_t value) const {
        if constexpr (dom_dim == 1) {
          return (diti - value);
        } else {
          diti_t rv = diti;
          diti_t tmp = value;
          for(int i=0; i<dom_dim; i++) {
            rv -= tmp;
            tmp = static_cast<diti_t>(tmp << dic_bits);
          }
          return rv;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Incriment all components in a a tuple
          @param diti The input tuple
          @param value Amout by which to increment each component
          @return New tuple with incrimented components */
      inline diti_t cuc_inc_all_crd(diti_t diti, dic_t value) const {
        if constexpr (dom_dim == 1) {
          return (diti + value);
        } else {
          diti_t rv = diti;
          diti_t tmp = value;
          for(int i=0; i<dom_dim; i++) {
            rv += tmp;
            tmp = static_cast<diti_t>(tmp << dic_bits);
          }
          return rv;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Set all components in a a tuple to a constant
          @param value The value for each component
          @return New tuple */
      inline diti_t cuc_set_all_crd(dic_t value) const {
        if constexpr (dom_dim == 1) {
          return value;
        } else {
          diti_t tmp = value;
          diti_t rv = 0;
          for(int i=0; i<dom_dim; i++)
            rv |= (tmp << (i * dic_bits));
          return rv;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute cross product points centered at diti and delta away
          Examples:
            - dom_dom==1 => {diti-delta, diti+delta}
            - dom_dom==2 => {diti+[-delta, -delta], diti+[-delta, +delta], diti+[+delta, -delta], diti+[+delta, +delta]}
            - dom_dom==3 => {diti+[-delta, -delta, -delta], ..., diti+[+delta, +delta, +delta]} -- this list will have 8 points
          @param diti Center coordinates for the cross product points
          @param delta The Distance for the cross product points
          @return Last of cross product points */
      diti_list_t cuc_two_cross(diti_t diti, dic_t delta) const {
        //  MJR TODO NOTE <2024-07-11T11:50:36-0500> cuc_two_cross: If diti is close to an corner, some result points may be out of range.
        diti_list_t rv;
        for(int i=0; i<(1 << dom_dim); i++) {
          diti_t tmp = diti;
          for(int j=0; j<dom_dim; j++) {
            if (i & (1 << j)) {
              tmp = cuc_inc_crd(tmp, j, delta);
            } else {
              tmp = cuc_dec_crd(tmp, j, delta);
            }
          }
          rv.push_back(tmp);
        }
        return rv;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute directional cross product points centered at diti and delta away
          @warning If diti is close to an corner, then some of the cross product points may be out of range!
          Examples:
            - dom_dom==1, index==0, direction==-1 => {diti-delta}
            - dom_dom==1, index==0, direction== 1 => {diti+delta}
            - dom_dom==2, index==0, direction==-1 => {diti+[-delta, -delta], diti+[-delta, +delta]                                                  }
            - dom_dom==2, index==0, direction== 1 => {                                                  diti+[+delta, -delta], diti+[+delta, +delta]}
            - dom_dom==2, index==1, direction==-1 => {diti+[-delta, -delta],                            diti+[+delta, -delta]                       }
            - dom_dom==2, index==1, direction== 1 => {                       diti+[-delta, +delta],                            diti+[+delta, +delta]}
          @param diti Center coordinates for the cross product points
          @param delta     The Distance for the cross product points.
          @param index     The index to hold constant.  Must be in [0, dom_dom-1].  No error checking.
          @param direction The direction on the given index.  Must be 1 or -1.  No error checking.
          @return List of cross product points */
      diti_list_t cuc_two_cross(diti_t diti, dic_t delta, int index, int direction) const {
        diti_list_t rv;
        if (direction != 1)
          direction = 0;
        for(int i=0; i<(1 << dom_dim); i++) {
          if (((i >> index) & 1) == direction) {
            diti_t tmp = diti;
            for(int j=0; j<dom_dim; j++) {
              if ((i >> j) & 1) {
                tmp = cuc_inc_crd(tmp, j, delta);
              } else {
                tmp = cuc_dec_crd(tmp, j, delta);
              }
            }
            rv.push_back(tmp);
          }
        }
        return rv;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute axis aligned cross points centered at diti and delta away
          Axis aligned cross points that fall outside of the domain will not be returned.  For example, if diti=0, then we get no points to the left.
          Examples:
            - dom_dom==1 => {diti-delta, diti+delta}
            - dom_dom==2 => {diti+[0, -delta], diti+[0, +delta], diti+[-delta, 0], diti+[+delta, 0]}
            - dom_dom==3 => {diti+[0, 0, -delta], diti+[0, 0, +delta],
                            diti+[0, -delta, 0], diti+[0, +delta, 0],
                            diti+[-delta, 0, 0], diti+[+delta, 0, 0]}
          @param diti Center coordinates for the cross
          @param delta The Distance for the cross points
          @return Last of cross product points */
      diti_list_t cuc_axis_cross(diti_t diti, dic_t delta) const {
        diti_list_t rv;
        for(int idx=0; idx<dom_dim; idx++) {
          for(int up=0; up<2; up++) {
            dic_t tmp = cuc_get_crd(diti, idx);
            if (up) {
              if ((dic_max-tmp) >= delta)
                rv.push_back(cuc_inc_crd(diti, idx, delta));
            } else {
              if (tmp >= delta)
                rv.push_back(cuc_dec_crd(diti, idx, delta));
            }
          }
        }
        return rv;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Easy way to treat rrpt_t & drpt_t types as indexable regardless of dom_dim & rng_dim. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /* Helper to index range values. */
      inline src_t rng_at(rrpt_t value, int index) const {
        if constexpr (rng_dim == 1) {
          return value;
        } else {
          return value[index];
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /* Helper to index domain values. */
      inline src_t dom_at(drpt_t value, int index) const {
        if constexpr (dom_dim == 1) {
          return value;
        } else {
          return value[index];
        }
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Packed Integer Tuple <-> std::array of Integer Tuple Conversions */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert an std::array of integer coordinate tuple to an packed integer coordinate tuple.
          @param dita The std::array of integer coordinate
          @return New packed integer tuple. */
      inline diti_t dita_to_diti(const dita_t& dita) const {
        if constexpr (dom_dim == 1) {
          return (dita[0]);
        } else {
          diti_t rv = 0;
          for(int i=0; i<dom_dim; i++)
            rv |= (dita[i] << (i * dic_bits));
          return rv;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a packed integer coordinate tuple to an std::array of integer coordinate tuple
          Note dita_t is *always* an std::array even when dom_dim==1.
          @param diti Packed integer coordinate tuple
          @return New std::array of integer coordinate tuple */
      inline dita_t diti_to_dita(diti_t diti) const {
        if constexpr (dom_dim == 1) {
          return (dita_t({diti}));
        } else {
          dita_t rv;
          for(int i=0; i<dom_dim; i++)
            rv[i] = (diti >> ((dom_dim-1-i) * dic_bits)) & diti_msk0;
          return rv;
        }
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Integer Tuple to Domain Space Tuple Conversion */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a packed integer coordinate tuple to the coorisponding coordinate in drpt_t
          Note this function might return a scalar (float, double, etc...) or a std::array!
          @param diti Input integer tuple
          @return domain space value. */
      inline drpt_t diti_to_drpt(diti_t diti) const {
        if constexpr (dom_dim == 1) {
          return (bbox_min + bbox_delta * diti);
        } else {
          drpt_t rv;
          rv.fill(static_cast<src_t>(0.0));
          for(int i=0; i<dom_dim; i++) {
            rv[i] = (bbox_min[i] + bbox_delta[i] * cuc_get_crd(diti, i));
          }
          return rv;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a packed integer coordinate tuple to the coorisponding std::array tuple of src_t
          Unlike diti_to_drpt, this function ALWAYS returns an std::array!
          @param diti Input integer tuple
          @return domain space value. */
      inline drta_t diti_to_drta(diti_t diti) const {
        drta_t rv;
        rv.fill(static_cast<src_t>(0.0));
        if constexpr (dom_dim == 1) {
          rv[0] = (bbox_min + bbox_delta * diti);
        } else {
          for(int i=0; i<dom_dim; i++) {
            rv[i] = (bbox_min[i] + bbox_delta[i] * cuc_get_crd(diti, i));
          }
        }
        return rv;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Function Sampleing

       Note that refine_grid() is a hybrid in that it can be used for both refinement and for sampling. */
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Sample a cell.
          @param cell Cell to sample
          @param func Function to use for samples */
      void sample_cell(diti_t cell, drpt2rrpt_func_t func) {
        if (sample_point_maybe(cell, func)) {
          for(auto const e: ccc_get_corners(cell)) {
            sample_point_maybe(e, func);
          }
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      void sample_cell(drpt2rrpt_func_t func) {
        sample_cell(ccc_get_top_cell(), func);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Sample a point if it has not already been sampled.
          @param diti Point at which to sample
          @param func Function to sample
          @return true if we sampled the point, and false otherwise. */
      inline bool sample_point_maybe(diti_t diti, drpt2rrpt_func_t func) {
        if ( !(vertex_exists(diti))) {
          drpt_t xvec = diti_to_drpt(diti);
          rrpt_t val = func(xvec);
          samples[diti] = val;
          return true;
        } else {
          return false;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Sample, or resample, a point.
          @param diti Point at which to sample
          @param func Function to sample */
      inline void sample_point(diti_t diti, drpt2rrpt_func_t func) {
        drpt_t xvec = diti_to_drpt(diti);
        rrpt_t val = func(xvec);
        samples[diti] = val;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Real Range Space Computation */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a point in the range space contains a NaN coordinate value
          @param val Value in  space */
      inline bool rrpt_is_nan(rrpt_t val) const {
        if constexpr (rng_dim == 1) {
          return std::isnan(val);
        } else {
          return (std::any_of(val.cbegin(), val.cend(), [](src_t v) { return (std::isnan(v)); }));
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Distance between two points in the range space using the infinity norm (max absolute value)
          @param val1 Value in range space
          @param val2 Value in range space */
      inline src_t rrpt_distance_inf(rrpt_t val1, rrpt_t val2) const {
        if constexpr (rng_dim == 1) {
          return std::abs(val1-val2);
        } else {
          src_t ret = 0;
          for(int i=0; i<rng_dim; ++i)
            ret = std::max(ret, std::abs(val1[i]-val2[i]));
          return ret;
        }
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Real Domain Space Computation */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Distance between two points in the domain space using the infinity norm (max absolute value)
          @param val1 Value in domain space
          @param val2 Value in domain space */
      inline src_t drpt_distance_inf(drpt_t val1, drpt_t val2) const {
        if constexpr (dom_dim == 1) {
          return std::abs(val1-val2);
        } else {
          src_t ret = 0;
          for(int i=0; i<dom_dim; ++i)
            ret = std::max(ret, std::abs(val1[i]-val2[i]));
          return ret;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the midpoint between two points in the domain space
          @param val1 Value in domain space
          @param val2 Value in domain space */
      inline drpt_t drpt_midpoint(drpt_t val1, drpt_t val2) const {
        if constexpr (dom_dim == 1) {
          return (val1+val2)/static_cast<src_t>(2.0);
        } else {
          drpt_t ret;
          for(int i=0; i<dom_dim; ++i)
            ret[i] = (val1[i] + val2[i])/static_cast<src_t>(2.0);
          return ret;
        }
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Tree/Cell Refinement

          Three general function smalpling stratigies are provided:
            - Manual: refine_once()
            - Uniformly sample a cell: refine_grid() & refine_recursive()
            - Recine recursively on a predicate: refine_recursive_cell_pred()
            - Refine leaves with a predicate:
              - Once: refine_leaves_once_if_cell_pred()
              - Repeatedly atomically: refine_leaves_atomically_if_cell_pred
              - Recursively: refine_leaves_recursive_cell_pred()
      */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Refine a cell if possible.
          @param cell Cell to refine -- no error checking!!
          @param func Function to use for samples
          @return 1 if cell was refined, and 0 otherwise -- i.e. the number of cells refined. */
      bool refine_once(diti_t cell, drpt2rrpt_func_t func) {
        auto children = ccc_get_children(cell);
        if (children.empty()) {
          return 0;
        } else {
          for(auto const c : children)   // ccc_get_children returns an empty list if the cell is too small to refine.
            sample_cell(c, func);
          return 1;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Sample a function uniformly across given cell to the given level.
          The given cell need not exist in the tree yet.

          @warning Will resample previously sampled points in the cell.

          @param cell        The cell to sample within
          @param level_delta The relative level at which to sample
                               level_delta=0 is equivalent to calling sample_cell(cell, func).
                               level_delta=1 is equivalent to calling sample_cell(cell, func) followed by refine_once(cell, func).
          @param func        Function to use for samples */
      void refine_grid(diti_t cell, int level_delta, drpt2rrpt_func_t func) {
        diti_t tmp = ccc_cell_get_corner_min(cell);
        if constexpr (dom_dim == 1) {
          dic_t  del = ccc_cell_half_width(cell) >> level_delta;
          for(int i=0; i<(1 << (1+level_delta))+1; i++) {
            sample_point(tmp, func);
            tmp += del;
          }
        } else if constexpr (dom_dim == 2) {
          dic_t  del = static_cast<dic_t>(ccc_cell_full_width(cell) >> level_delta);
          // Corners
          for(dic_t i=0; i<(1 << level_delta)+1; i++) {
            diti_t tmp2 = cuc_inc_crd(tmp, 0, i*del);
            for(dic_t j=0; j<(1 << level_delta)+1; j++) {
              diti_t tmp3  = cuc_inc_crd(tmp2, 1, j*del);
              sample_point(tmp3, func);
            }
          }
          // Centers
          tmp = cuc_inc_all_crd(tmp, del/2);
          for(dic_t i=0; i<(1 << level_delta); i++) {
            diti_t tmp2 = cuc_inc_crd(tmp, 0, i*del);
            for(dic_t j=0; j<(1 << level_delta); j++) {
              diti_t tmp3  = cuc_inc_crd(tmp2, 1, j*del);
              sample_point(tmp3, func);
            }
          }
        } else if constexpr (dom_dim == 3) {
          dic_t  del = static_cast<dic_t>(ccc_cell_full_width(cell) >> level_delta);
          // Corners
          for(dic_t i=0; i<(1 << level_delta)+1; i++) {
            diti_t tmp2 = cuc_inc_crd(tmp, 0, i*del);
            for(dic_t j=0; j<(1 << level_delta)+1; j++) {
              diti_t tmp3  = cuc_inc_crd(tmp2, 1, j*del);
              for(dic_t k=0; k<(1 << level_delta)+1; k++) {
                diti_t tmp4  = cuc_inc_crd(tmp3, 2, k*del);
                sample_point(tmp4, func);
              }
            }
          }
          // Centers
          tmp = cuc_inc_all_crd(tmp, del/2);
          for(dic_t i=0; i<(1 << level_delta); i++) {
            diti_t tmp2 = cuc_inc_crd(tmp, 0, i*del);
            for(dic_t j=0; j<(1 << level_delta); j++) {
              diti_t tmp3  = cuc_inc_crd(tmp2, 1, j*del);
              for(dic_t k=0; k<(1 << level_delta); k++) {
                diti_t tmp4  = cuc_inc_crd(tmp3, 2, k*del);
                sample_point(tmp4, func);
              }
            }
          }
        } else {
          std::cout << "ERROR: refine_grid with dom_dim>3 not supported (use refine_recursive)!" << std::endl;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      void refine_grid(int level_delta, drpt2rrpt_func_t func) {
        refine_grid(ccc_get_top_cell(), level_delta, func);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Refine a cell until refined cells reach specified level

          @note Unlike most of the other refine functions, this one will sample the given cell. i.e. the provided cell need not be already sampled.  This
                difference is because this function is frequently called to initially sample a tree.

          @param cell Cell to refine
          @param level Maximum level of refined cells.  -1 means refine to the limit.
          @param func Function to use for samples */
      void refine_recursive(diti_t cell, int level, drpt2rrpt_func_t func) {
        sample_cell(cell, func);
        if ((level < 0) || (ccc_cell_level(cell) < level)) {
          for(auto const c : ccc_get_children(cell)) {
            sample_cell(c, func);
            refine_recursive(c, level, func);
          }
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      void refine_recursive(int level, drpt2rrpt_func_t func) {
        refine_recursive(ccc_get_top_cell(), level, func);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Refine a cells matching predicate until refined cells reach specified level

          @warning cell must exist (i.e. must be sampled)

          @note Nothing happens if pred isn't true for cell.  i.e. this is a worker function to refine a cell with a known condition until the condition is
          no longer true, or a maximum level has been reached.  This is normally the wrong function to call with the ccc_get_top_cell().  If you wish to
          refine recursively all cells in a tree that match a predicate, then use refine_leaves_recursive_cell_pred().

          @param cell Cell to refine
          @param level Maximum level of refinded cells.  -1 means refine to the limit.
          @param func Function to use for samples
          @param pred Predicate function. */
      void refine_recursive_cell_pred(diti_t cell, int level, drpt2rrpt_func_t func, diti2bool_func_t pred) {
        if ((level < 0) || (ccc_cell_level(cell) < level)) {
          if (pred(cell)) {
            refine_once(cell, func);
            for(auto const c : ccc_get_children(cell)) {
              refine_recursive_cell_pred(c, level, func, pred);
            }
          }
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Refine a cells matching predicate until refined cells reach specified level
          @warning cell must exist (i.e. must be sampled)
          @param cell Cell to refine
          @param level Maximum level of refinded cells.  -1 means refine to the limit.
          @param func Function to use for samples
          @param pred Predicate function. */
      void refine_leaves_recursive_cell_pred(diti_t cell, int level, drpt2rrpt_func_t func, diti2bool_func_t pred) {
        for(auto c: get_leaf_cells(cell))
          refine_recursive_cell_pred(c, level, func, pred);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      void refine_leaves_recursive_cell_pred(int level, drpt2rrpt_func_t func, diti2bool_func_t pred) {
        refine_leaves_recursive_cell_pred(ccc_get_top_cell(), level, func, pred);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Refine each leaf cell if pred returns true and the cell level is less than the given level value.

          All leaf cells are tested with the predicate first.  Then the list of cells requiring refinement are refined.  This guarantees that no leaf cell
          will be refined more than once even if the refinement process would impact the value of the predicate.  This can produce qualitatively different
          behavior than refine_leaves_recursive_cell_pred() when used with predicates like cell_is_unbalanced().

          @param cell  Input cell. Must be a valid cell. -- no error checking.
          @param level Maximum level of refinded cells.  -1 means refine to the limit.
          @param func  Function to sample
          @param pred  Predicate function. */
      int  refine_leaves_once_if_cell_pred(diti_t cell, int level, drpt2rrpt_func_t func, diti2bool_func_t pred) {
        diti_list_t cells_to_check  = get_leaf_cells(cell);
        diti_list_t cells_to_refine;
        for(auto c: cells_to_check)
          if (pred(c))
            if ((level < 0) || (ccc_cell_level(c) < level))
              cells_to_refine.push_back(c);
        int refined_count = 0;
        for(auto c: cells_to_refine)
          refined_count += refine_once(c, func);
        return refined_count;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      int  refine_leaves_once_if_cell_pred(int level, drpt2rrpt_func_t func, diti2bool_func_t pred) {
        return refine_leaves_once_if_cell_pred(ccc_get_top_cell(), level, func, pred);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Repeatedly refine leaves of the given cell for which pred is true.

          Apply refine_leaves_once_if_cell_pred() repeatedly untill no new cells are refined.

          @param cell  Input cell. Must be a valid cell. -- no error checking.
          @param level Maximum level of refinded cells.  -1 means refine to the limit.
          @param func  Function to sample
          @param pred  Predicate function. */
      int refine_leaves_atomically_if_cell_pred(diti_t cell, int level, drpt2rrpt_func_t func, diti2bool_func_t pred) {
        int refined_count = 0;
        while (0 < (refined_count = refine_leaves_once_if_cell_pred(cell, level, func, pred)))
          ;
        return refined_count;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      int refine_leaves_atomically_if_cell_pred(int level, drpt2rrpt_func_t func, diti2bool_func_t pred) {
        return refine_leaves_atomically_if_cell_pred(ccc_get_top_cell(), level, func, pred);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Refine a cells with NaNs until refined cells reach specified level

          This is a convenience function combining refine_recursive_cell_pred(), cell_vertex_is_nan(), & ccc_get_top_cell().

          @param level Maximum level of refinded cells.  -1 means refine to the limit.
          @param func Function to use for samples */
      void refine_recursive_if_cell_vertex_is_nan(int level, drpt2rrpt_func_t func) {
        refine_leaves_recursive_cell_pred(level, func, std::bind_front(&MR_rect_tree::cell_vertex_is_nan, this));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Refine leaf cells if they are unbalanced at the given level.

          This is a convenience function combining refine_leaves_atomically_if_cell_pred(), cell_is_unbalanced(), & ccc_get_top_cell().

          The primary use case is to demonstrate the step wise balancing of a tree.

          @param level_delta The Level.
          @param func        Function to sample */
      int refine_leaves_once_if_unbalanced(int level_delta, drpt2rrpt_func_t func) {
        return refine_leaves_once_if_cell_pred(ccc_get_top_cell(), -1, func, std::bind_front(&cell_is_unbalanced, this, level_delta));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Balance the cell to the given level.

          This is a convenience function combining refine_leaves_atomically_if_cell_pred(), cell_is_unbalanced(), & ccc_get_top_cell().

          @param level_delta  The Level.
          @param func         Function to sample */
      void balance_tree(int level_delta, drpt2rrpt_func_t func) {
        refine_leaves_atomically_if_cell_pred(ccc_get_top_cell(), -1, func, std::bind_front(&this_t::cell_is_unbalanced, this, level_delta));
      }
      //@}


      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cell Predicates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if cell coordinates are in range to be a cell center
          @param cell Input cell
          @return True if cell is in bounds, false otherwise */
      inline bool cell_good_cords(diti_t cell) const {
        if constexpr (dom_dim == 1) {
          return ((0 < cell) && (dic_max > cell));
        } else {
          diti_t working_tuple = cell;
          for(int i=0; i<dom_dim; i++) {
            dic_t cur_dic = diti_msk0 & working_tuple;
            if ((cur_dic <= 0) || (cur_dic >= dic_max)) {
              return false;
            }
            working_tuple = working_tuple >> dic_bits;
          }
          return true;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell has been sampled.
          @warning Simply checks that cell has been sampled -- identical to vertex_exists().
          @param cell Input cell*/
      inline bool cell_exists(diti_t cell) const {
        return (vertex_exists(cell));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell has been sampled.
          @warning Simply checks that cell has been sampled -- identical to vertex_exists().
          @param cell Input cell*/
      inline bool cell_is_sampled(diti_t cell) const {
        diti_list_t verts = ccc_get_vertexes(cell);
        return (std::all_of(verts.cbegin(), verts.cend(), [this](diti_t i) { return (vertex_exists(i)); }));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell has a vertex with a NaN value for a sample
          @param cell Input cell */
      inline bool cell_vertex_is_nan(diti_t cell) {
        diti_list_t verts = ccc_get_vertexes(cell);
        return (std::any_of(verts.cbegin(), verts.cend(), [this](diti_t i) { return (vertex_is_nan(i)); }));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell has an corner with a NaN value for a sample
          @param cell Input cell */
      inline bool cell_corner_is_nan(diti_t cell) {
        diti_list_t corners = ccc_get_corners(cell);
        return (std::any_of(corners.cbegin(), corners.cend(), [this](diti_t i) { return (vertex_is_nan(i)); }));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell has children
          @warning Assumes the tree is well formed, and only checks that the "lower left" child's center been sampled!
          @param cell Input cell */
      inline bool cell_has_child(diti_t cell) const {
        return (cell_can_have_children(cell) && cell_exists(cuc_dec_all_crd(cell, ccc_cell_quarter_width(cell))));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell has no children
          @warning Assumes the tree is well formed, and only checks that the "lower left" child's center been sampled!
          @param cell Input cell */
      inline bool cell_has_no_child(diti_t cell) const {
        return ( !(cell_has_child(cell)));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell can have children (it not of minimal size)
          @param cell Input cell */
      inline bool cell_can_have_children(diti_t cell) const {
        return (ccc_cell_level(cell) < (max_level - 1));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if the specified neighbor cell exixts (is smapled).
          @param index     The index of the axis.  Must be in [0, dom_dim-1].  No error checking.
          @param direction The direction on the given index.  Must be 1 or -1.  No error checking.
          @param cell Input cell */
      inline bool cell_has_neighbor(diti_t cell, int index, int direction) {
        diti_t tmp = ccc_get_neighbor(cell, index, direction);
        if (tmp)
          return cell_is_sampled(cell);
        else
          return false;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell crosses, or is on, a signed distance function boundry.
          "Crossing" is defined as: The center is zero or any corner has a diffrent sign from the center.
          @warning Incorrect result if ALL vertexes of the cell are zero.
          @param cell Input Cell
          @param sdf Signed distance function
          @return true if the cell crosses, or is on, the signed distance function boundry. */
      inline bool cell_cross_sdf(diti_t cell, drpt2real_func_t sdf) {
        /* The algorithm below directly expresses the RHS of the following iff which is equivalent to the LHS (and the statement in the documentation).
           @f[
           (\mathrm{sgn}(\vec{\mathbf{c}})=0)\lor(\exists \vec{\mathbf{v}}\in E(\vec{\mathbf{c}})\,\mathrm{st}\,\mathrm{sgn}(\vec{\mathbf{c}})\ne\mathrm{sgn}(\vec{\mathbf{v}}))
           \Longleftrightarrow
           (\exists \vec{\mathbf{v}}\in E(\vec{\mathbf{c}})\,\mathrm{st}\,\mathrm{sgn}(\vec{\mathbf{v}})=0)
           \lor
           (\exists \vec{\mathbf{v_1}},\vec{\mathbf{v_2}}\in E(\vec{\mathbf{c}})\,\mathrm{st}\,\mathrm{sgn}(\vec{\mathbf{v1}})\ne\mathrm{sgn}(\vec{\mathbf{v2}}))
           @f]
        */
        int center_sign = mjr::math::sfun::sgn(sdf(diti_to_drpt(cell)));
        if (center_sign == 0)
          return true;
        for(diti_t& v: ccc_get_corners(cell))
          if (center_sign != mjr::math::sfun::sgn(sdf(diti_to_drpt(v))))
            return true;
        return false;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell contains, or is close to, a point in the domain.
          @param domain_point The point in the domain
          @param epsilon      How close the point must be
          @param cell Input Cell
          @return true if a cell contains, or is close to, a domain_point. */
      inline bool cell_near_domain_point(drpt_t domain_point, src_t epsilon, diti_t cell) {
        drpt_t min_drpt = diti_to_drpt(ccc_cell_get_corner_min(cell));
        for(int i=0; i<dom_dim; i++)
          if (dom_at(min_drpt, i)-epsilon > dom_at(domain_point, i))
            return false;
        drpt_t max_drpt = diti_to_drpt(ccc_cell_get_corner_max(cell));
        for(int i=0; i<dom_dim; i++)
          if (dom_at(max_drpt, i)+epsilon < dom_at(domain_point, i))
            return false;
        return true;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell crosses the given domain component value
          @param cell Input Cell
          @param domain_index The index of the domain component we are testing
          @param domain_level The level, or value, of the domain component we are testing
          @param epsilon      Used to fuzz floating point comparisons
          @return true if the cell crosses the domain level. */
      inline bool cell_near_domain_level(diti_t cell, int domain_index, src_t domain_level, src_t epsilon) {
        return ( (dom_at(diti_to_drpt(ccc_cell_get_corner_min(cell)), domain_index) < domain_level+epsilon) &&
                 (dom_at(diti_to_drpt(ccc_cell_get_corner_max(cell)), domain_index) > domain_level-epsilon) );
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell is below the given domain component value
          @param cell Input Cell
          @param domain_index The index of the domain component we are testing
          @param domain_level The level, or value, of the domain component we are testing
          @return true if the cell below the domain level. */
      inline bool cell_below_domain_level(diti_t cell, int domain_index, src_t domain_level) {
        return (dom_at(diti_to_drpt(ccc_cell_get_corner_max(cell)), domain_index) < domain_level);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell above the given domain component value
          @param cell Input Cell
          @param domain_index The index of the domain component we are testing
          @param domain_level The level, or value, of the domain component we are testing
          @return true if the cell above the domain level. */
      inline bool cell_above_domain_level(diti_t cell, int domain_index, src_t domain_level) {
        return ((dom_at(diti_to_drpt(ccc_cell_get_corner_min(cell)), domain_index) > domain_level));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell crosses the given range component value.
          See ::cell_cross_sdf for algorithm notes.
          @param cell Input Cell
          @param range_index The index of the range component we are testing
          @param range_level The level, or value, of the range component we are testing
          @return true if the cell crosses the range level. */
      inline bool cell_cross_range_level(diti_t cell, int range_index, src_t range_level) {
        int center_sign = mjr::math::sfun::sgn(rng_at(samples[cell], range_index)-range_level);
        if (center_sign == 0)
          return true;
        for(diti_t& v: ccc_get_corners(cell))
          if (center_sign != mjr::math::sfun::sgn(rng_at(samples[v], range_index)-range_level))
            return true;
        return false;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell is below given range component value
          @param cell Input Cell
          @param range_index The index of the range component we are testing
          @param range_level The level, or value, of the range component we are testing
          @return true if the cell is below the range level. */
      inline bool cell_below_range_level(diti_t cell, int range_index, src_t range_level) {
        diti_list_t verts = ccc_get_vertexes(cell);
        return std::all_of(verts.cbegin(), verts.cend(), [this, range_index, range_level](int i) { return (rng_at(samples[i], range_index) < range_level); });
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell is above given range component value
          @param cell Input Cell
          @param range_index The index of the range component we are testing
          @param range_level The level, or value, of the range component we are testing
          @return true if the cell is above the range level. */
      inline bool cell_above_range_level(diti_t cell, int range_index, src_t range_level) {
        diti_list_t verts = ccc_get_vertexes(cell);
        return std::all_of(verts.cbegin(), verts.cend(), [this, range_index, range_level](int i) { return (rng_at(samples[i], range_index) > range_level); });
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a cell is unbalanced at the given level
          @param cell        Input Cell
          @param level_delta Signed distance function
          @return true if the cell is unbalanced at the given level. */
      bool cell_is_unbalanced(int level_delta, diti_t cell) {
        //  MJR TODO NOTE <2024-07-11T16:04:10-0500> cell_is_unbalanced: Optimize run time (nix the use of get_smallest_neighbor_level)
        int max_level_neighbor     = get_smallest_neighbor_level(cell);
        int current_ccc_cell_level = ccc_cell_level(cell);
        if (max_level_neighbor > current_ccc_cell_level+level_delta)
          return true;
        else
          return false;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Vertex Predicates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Test if a vertex is NaN or, when it is an std::array, if it contains a NaN element
          @param vertex Input vertex */
      inline bool vertex_is_nan(diti_t vertex) {
        if constexpr (rng_dim == 1) {
          return std::isnan(samples[vertex]);
        } else {
          return (std::any_of(samples[vertex].cbegin(), samples[vertex].cend(), [this](src_t v) { return (std::isnan(v)); }));
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Point has been sampled.
          @param vertex Input vertex */
      inline bool vertex_exists(diti_t vertex) const {
        return (samples.contains(vertex));
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Extract Cells */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Extract a list of all leaf cells starting from the given cell
          @param cell Starting cell */
      diti_list_t get_leaf_cells(diti_t cell) const {
        diti_list_t rv;
        if (cell_has_child(cell)) {
          for(auto const c : ccc_get_children(cell)) {
            auto mv = get_leaf_cells(c);
            std::move(mv.begin(), mv.end(), std::back_inserter(rv));
          }
        } else {
          rv.push_back(cell);
        }
        return rv;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      diti_list_t get_leaf_cells() const {
        return get_leaf_cells(ccc_get_top_cell());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Extract a list of all leaf cells starting from the given cell that match the given predicate
          @warning The given cell need not match the predicate.
          @param cell Starting cell
          @param pred Predicate function. */
      diti_list_t get_leaf_cells_pred(diti_t cell, diti2bool_func_t pred) const {
        diti_list_t cells_to_return;
        for(auto c: get_leaf_cells(cell))
          if (pred(c))
            cells_to_return.push_back(c);
        return cells_to_return;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Extract a list of all leaf cells starting from the given cell
          @param cell      Input cell. Must be a valid cell. -- no error checking.
          @param index     The index of the axis.  Must be in [0, dom_dim-1].  No error checking.
          @param direction The direction on the given index.  Must be 1 or -1.  No error checking. */
      diti_list_t get_leaf_cells(diti_t cell, int index, int direction) const {
        diti_list_t rv;
        if (cell_has_child(cell)) {
          for(auto const c : ccc_get_children(cell, index, direction)) {
            auto mv = get_leaf_cells(c, index, direction);
            std::move(mv.begin(), mv.end(), std::back_inserter(rv));
          }
        } else {
          rv.push_back(cell);
        }
        return rv;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Count the number of leaf cells starting from the given cell
          @param cell Starting cell */
      int count_leaf_cells(diti_t cell) const {
        int rv=0;
        if (cell_has_child(cell)) {
          for(auto const c : ccc_get_children(cell)) {
            int mv = count_leaf_cells(c);
            rv+=mv;
          }
        } else {
          rv++;
        }
        return rv;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return a list of sampled neighbors to the given cell.
          @warning Unlike ccc_get_neighbor(), this member returns a *list* of *existing* (sampled) cell centers.
          @param cell      Input cell. Must be a valid cell. -- no error checking.
          @param index     The index of the axis.  Must be in [0, dom_dim-1].  No error checking.
          @param direction The direction on the given index.  Must be 1 or -1.  No error checking.
          @return  List of existing neighbors. */
      diti_list_t get_existing_neighbor(diti_t cell, int index, int direction) const {
        diti_t tmp = ccc_get_neighbor(cell, index, direction);
        if ( (tmp!=0) && cell_exists(tmp)) {
          if (cell_has_child(tmp)) {
              return get_leaf_cells(tmp, index, -direction);
            } else {
              return diti_list_t({tmp});
            }
        }
        return diti_list_t();
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the level of the smallest, existing neighbor.
          Returns -1 if no neighbors exist.
          @param cell      Input cell. Must be a valid cell. -- no error checking. */
      int get_smallest_neighbor_level(diti_t cell) const {
        int maximum_level = -1;
        int start_level = ccc_cell_level(cell);
        for(int aix=0; aix<dom_dim; aix++) {
          for(int dir=-1; dir<2; dir+=2) {
            diti_t nbr0 = ccc_get_neighbor(cell, aix, dir);
            if ( (nbr0!=0) && cell_exists(nbr0)) {
              if (cell_has_child(nbr0)) {
                diti_list_t nbr_leafs = get_leaf_cells(nbr0, aix, -dir);
                for(auto nbr_leaf: nbr_leafs) {
                  int nbr_leaf_lvl = ccc_cell_level(nbr_leaf);
                  if (nbr_leaf_lvl > maximum_level)
                    maximum_level = nbr_leaf_lvl;
                }
              } else {
                if (start_level > maximum_level)
                  maximum_level = start_level;
              }
            }
          }
        }
        return maximum_level;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Debug */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a diti_t to a string representation
          @param diti      The diti_t to convert
          @param do_hex    Print hex if true
          @param include_domain When true, also include the real coordinates in the string representation. */
      std::string diti_to_string(diti_t diti, bool include_domain, bool do_hex) const {
        std::ostringstream convert;
        int pwid = (do_hex ? (max_level+1)/4+1 : (max_level+1)/2+1);
        dita_t tmp1 = diti_to_dita(diti);
        for(int i=0; i<dom_dim; i++)
          convert << std::setw(pwid) << std::setfill('0') << (do_hex ? std::hex : std::dec ) << (uint32_t)tmp1[dom_dim-1-i] << std::dec << std::setw(0) << std::setfill('\0') << " ";
        if (include_domain) {
          drta_t tmp3 = diti_to_drta(diti);
          convert << "[ ";
          for (const auto tmp4 : tmp3)
            convert << std::setprecision(5) << static_cast<src_t>(tmp4) << " ";
          convert << "] ";
        }
          return(convert.str());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      std::string diti_to_string(diti_t diti, bool include_domain) const { return diti_to_string(diti, include_domain, true); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      std::string diti_to_string(diti_t diti)                      const { return diti_to_string(diti, true,          true); }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a drpt_t to a string representation */
      std::string drpt_to_string(drpt_t x) const {
        std::ostringstream convert;
        if constexpr (dom_dim == 1) {
          convert << "[ " << x << " ]";
        } else {
          convert << "[ ";
          for(auto c: x)
            convert << std::setprecision(5) << static_cast<src_t>(c) << " ";
          convert << "]";
        }
        return(convert.str());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert rrpt value to a string representation */
      std::string rrpt_to_string(rrpt_t x) const {
        std::ostringstream convert;
        if constexpr (rng_dim == 1) {
          convert << "[ " << x << " ]";
        } else {
          convert << "[ ";
          for(auto c: x)
            convert << std::setprecision(5) << static_cast<src_t>(c) << " ";
          convert << "]";
        }
          return(convert.str());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump tree to STDOUT
          @param max_num_print Maximum number of samples to print.  Use 0 to print all samples. */
      void dump_tree(int max_num_print) const {
        std::cout << "Meta Data" << std::endl;
        std::cout << "  get_bbox_min ... " << drpt_to_string(get_bbox_min()) << std::endl;
        std::cout << "  get_bbox_max ... " << drpt_to_string(get_bbox_max()) << std::endl;
        std::cout << "  dom_dim ........ " << dom_dim << std::endl;
        std::cout << "  rng_dim ........ " << rng_dim << std::endl;
        std::cout << "  max_level ...... " << max_level << std::endl;
        std::cout << "  size icrd Cmp .. " << sizeof(dic_t)  << std::endl;
        std::cout << "  size icrd Tup .. " << sizeof(diti_t) << std::endl;
        std::cout << "  Samples ........ " << samples.size() << std::endl;
        std::cout << "  Leaf Cells ..... " << count_leaf_cells(ccc_get_top_cell()) << std::endl;
        std::cout << "Samples" << std::endl;
        int num_printed = 0;
        for (const auto& kvp : samples) {
          std::cout << "  c=" << diti_to_string(kvp.first, true);
          std::cout << " v=" << rrpt_to_string(kvp.second) << std::endl;
          num_printed++;
          if ((max_num_print > 0) && (num_printed >= max_num_print)) {
            std::cout << "Maximum number of samples reached.  Halting tree dump." << std::endl;
            break;
          }
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump tree points to file -- one tuple per line with domain coordinates followed by range values. */
      int dump_tree_datafile(std::string file_name) const {
        std::ofstream out_stream;
        out_stream.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);
        if (out_stream.is_open()) {
          out_stream.imbue(std::locale::classic());
        } else {
          std::cout << "ERROR(write_xml_vtk): Could not open file!" << std::endl;
          return 1;
        }
        for (const auto& kvp : samples) {
          for(int i=0; i<dom_dim; i++)
            out_stream << std::setprecision(5) << dom_at(diti_to_drpt(kvp.first), i) << " ";
          for(int i=0; i<rng_dim; i++)
            out_stream << std::setprecision(5) << rng_at(kvp.second, i) << " ";
          out_stream << std::endl;
        }
        out_stream.close();
        return 0;
      }

      //@}

  };

    //--------------------------------------------------------------------------------------------------------------------------------------------------------
    /* 7-bit per coordinate */
    typedef mjr::MR_rect_tree<7, double, 1, 1> tree7b1d1rT;
    typedef mjr::MR_rect_tree<7, double, 2, 1> tree7b2d1rT;
    typedef mjr::MR_rect_tree<7, double, 3, 1> tree7b3d1rT;
    typedef mjr::MR_rect_tree<7, double, 4, 1> tree7b4d1rT;
    typedef mjr::MR_rect_tree<7, double, 5, 1> tree7b5d1rT; // For testing

    typedef mjr::MR_rect_tree<7, double, 1, 2> tree7b1d2rT;
    typedef mjr::MR_rect_tree<7, double, 2, 2> tree7b2d2rT;
    typedef mjr::MR_rect_tree<7, double, 3, 2> tree7b3d2rT;
    typedef mjr::MR_rect_tree<7, double, 4, 2> tree7b4d2rT;

    typedef mjr::MR_rect_tree<7, double, 1, 3> tree7b1d3rT;
    typedef mjr::MR_rect_tree<7, double, 2, 3> tree7b2d3rT;
    typedef mjr::MR_rect_tree<7, double, 3, 3> tree7b3d3rT;
    typedef mjr::MR_rect_tree<7, double, 4, 3> tree7b4d3rT;

    typedef mjr::MR_rect_tree<7, double, 1, 4> tree7b1d4rT;
    typedef mjr::MR_rect_tree<7, double, 2, 4> tree7b2d4rT;
    typedef mjr::MR_rect_tree<7, double, 3, 4> tree7b3d4rT;
    typedef mjr::MR_rect_tree<7, double, 4, 4> tree7b4d4rT;


    //--------------------------------------------------------------------------------------------------------------------------------------------------------
    /* 15-bit per coordinate */

    typedef mjr::MR_rect_tree<15, double, 1, 1> tree15b1d1rT;  // Curve defined by [x, f(x)]
    typedef mjr::MR_rect_tree<15, double, 2, 1> tree15b2d1rT;  // Surface defined by [x, y, f(x, y)]
                                                               // Monochrome image defined by [x, y] and colored by f(x, y)
                                                               // 2D scalar field
    typedef mjr::MR_rect_tree<15, double, 3, 1> tree15b3d1rT;  // 3D scalar field
    typedef mjr::MR_rect_tree<15, double, 4, 1> tree15b4d1rT;

    typedef mjr::MR_rect_tree<15, double, 1, 2> tree15b1d2rT;  // Parametric Plane Curve defined by [x(t), y(t)]
    typedef mjr::MR_rect_tree<15, double, 2, 2> tree15b2d2rT;  // Complex function magnutude surface defined by [Re(z), Im(z), mag(z)]
    typedef mjr::MR_rect_tree<15, double, 3, 2> tree15b3d2rT;
    typedef mjr::MR_rect_tree<15, double, 4, 2> tree15b4d2rT;

    typedef mjr::MR_rect_tree<15, double, 1, 3> tree15b1d3rT;  // Parametric Space Curve defined by [x(t), y(t), z(t)]
    typedef mjr::MR_rect_tree<15, double, 2, 3> tree15b2d3rT;  // Parametric Surface defined by [x(u, v), y(u, v), z(u, v)]
                                                               // RGB image defined by [x, y] and colored by [f_1(x, y), f_2(x, y), f_3(x, y)]
    typedef mjr::MR_rect_tree<15, double, 3, 3> tree15b3d3rT;
    typedef mjr::MR_rect_tree<15, double, 4, 3> tree15b4d3rT;

    typedef mjr::MR_rect_tree<15, double, 1, 4> tree15b1d4rT;  //  curve defined by [x, f(x)] with point color
                                                               // Parametric Plane Curve defined by [x(t), y(t)] with velocity vector
                                                               // 3D scalar field with piont color
    typedef mjr::MR_rect_tree<15, double, 2, 4> tree15b2d4rT;  // Complex function magnutude surface defined by [Re(z), Im(z), mag(z)] with point color XOR surface normals
                                                               // Surface defined by [x, y, f(x, y)] with surface normal XOR point color
                                                               // 2D scalar field with point color
    typedef mjr::MR_rect_tree<15, double, 3, 4> tree15b3d4rT;
    typedef mjr::MR_rect_tree<15, double, 4, 4> tree15b4d4rT;

    typedef mjr::MR_rect_tree<15, double, 1, 5> tree15b1d5rT;  // Parametric Plane Curve defined by [x(t), y(t)] with point color
    typedef mjr::MR_rect_tree<15, double, 2, 5> tree15b2d5rT;
    typedef mjr::MR_rect_tree<15, double, 3, 5> tree15b3d5rT;
    typedef mjr::MR_rect_tree<15, double, 4, 5> tree15b4d5rT;

    typedef mjr::MR_rect_tree<15, double, 1, 6> tree15b1d6rT;  // Parametric Space Curve defined by [x(t), y(t), z(t)] with velocity vector XOR point color
    typedef mjr::MR_rect_tree<15, double, 2, 6> tree15b2d6rT;  // Parametric Surface defined by [x(u, v), y(u, v), z(u, v)] with surface normal XOR point color
    typedef mjr::MR_rect_tree<15, double, 3, 6> tree15b3d6rT;
    typedef mjr::MR_rect_tree<15, double, 4, 6> tree15b4d6rT;

    typedef mjr::MR_rect_tree<15, double, 1, 7> tree15b1d7rT;  // Parametric Plane Curve defined by [x(t), y(t)] with velocity vector AND point color
    typedef mjr::MR_rect_tree<15, double, 2, 7> tree15b2d7rT;  // Surface defined by [x, y, f(x, y)] with surface normal AND point color
                                                               // Complex function magnutude surface defined by [Re(z), Im(z), mag(z)] with point color AND surface normals
    typedef mjr::MR_rect_tree<15, double, 3, 7> tree15b3d7rT;
    typedef mjr::MR_rect_tree<15, double, 4, 7> tree15b4d7rT;

    typedef mjr::MR_rect_tree<15, double, 1, 8> tree15b1d8rT;
    typedef mjr::MR_rect_tree<15, double, 2, 8> tree15b2d8rT;
    typedef mjr::MR_rect_tree<15, double, 3, 8> tree15b3d8rT;
    typedef mjr::MR_rect_tree<15, double, 4, 8> tree15b4d8rT;

    typedef mjr::MR_rect_tree<15, double, 1, 9> tree15b1d9rT;  // Parametric Space Curve defined by [x(t), y(t), z(t)] with velocity vector AND point color
    typedef mjr::MR_rect_tree<15, double, 2, 9> tree15b2d9rT;  // Parametric Surface defined by [x(u, v), y(u, v), z(u, v)] with surface normal AND point color
    typedef mjr::MR_rect_tree<15, double, 3, 9> tree15b3d9rT;
    typedef mjr::MR_rect_tree<15, double, 4, 9> tree15b4d9rT;

    typedef mjr::MR_rect_tree<15, double, 1, 10> tree15b1d10rT;
    typedef mjr::MR_rect_tree<15, double, 2, 10> tree15b2d10rT;
    typedef mjr::MR_rect_tree<15, double, 3, 10> tree15b3d10rT;
    typedef mjr::MR_rect_tree<15, double, 4, 10> tree15b4d10rT;

    typedef mjr::MR_rect_tree<15, double, 1, 11> tree15b1d11rT;
    typedef mjr::MR_rect_tree<15, double, 2, 11> tree15b2d11rT;
    typedef mjr::MR_rect_tree<15, double, 3, 11> tree15b3d11rT;
    typedef mjr::MR_rect_tree<15, double, 4, 11> tree15b4d11rT;

    typedef mjr::MR_rect_tree<15, double, 1, 12> tree15b1d12rT;
    typedef mjr::MR_rect_tree<15, double, 2, 12> tree15b2d12rT;
    typedef mjr::MR_rect_tree<15, double, 3, 12> tree15b3d12rT;
    typedef mjr::MR_rect_tree<15, double, 4, 12> tree15b4d12rT;

    typedef mjr::MR_rect_tree<15, double, 1, 13> tree15b1d13rT;
    typedef mjr::MR_rect_tree<15, double, 2, 13> tree15b2d13rT;
    typedef mjr::MR_rect_tree<15, double, 3, 13> tree15b3d13rT;
    typedef mjr::MR_rect_tree<15, double, 4, 13> tree15b4d13rT;

    typedef mjr::MR_rect_tree<15, double, 1, 14> tree15b1d14rT;
    typedef mjr::MR_rect_tree<15, double, 2, 14> tree15b2d14rT;
    typedef mjr::MR_rect_tree<15, double, 3, 14> tree15b3d14rT;
    typedef mjr::MR_rect_tree<15, double, 4, 14> tree15b4d14rT;

    typedef mjr::MR_rect_tree<15, double, 1, 15> tree15b1d15rT;
    typedef mjr::MR_rect_tree<15, double, 2, 15> tree15b2d15rT;
    typedef mjr::MR_rect_tree<15, double, 3, 15> tree15b3d15rT;
    typedef mjr::MR_rect_tree<15, double, 4, 15> tree15b4d15rT;

    //--------------------------------------------------------------------------------------------------------------------------------------------------------
    /* 31-bit per coordinate */
    typedef mjr::MR_rect_tree<31, double, 1, 1> tree31b1d1rT;
    typedef mjr::MR_rect_tree<31, double, 2, 1> tree31b2d1rT;

    typedef mjr::MR_rect_tree<31, double, 1, 2> tree31b1d2rT;
    typedef mjr::MR_rect_tree<31, double, 2, 2> tree31b2d2rT;

    typedef mjr::MR_rect_tree<31, double, 1, 3> tree31b1d3rT;
    typedef mjr::MR_rect_tree<31, double, 2, 3> tree31b2d3rT;

    typedef mjr::MR_rect_tree<31, double, 1, 4> tree31b1d4rT;
    typedef mjr::MR_rect_tree<31, double, 2, 4> tree31b2d4rT;

    //--------------------------------------------------------------------------------------------------------------------------------------------------------
    /* 63-bit per coordinate */

    typedef mjr::MR_rect_tree<63, double, 1, 1> tree63b1d1rT;
    typedef mjr::MR_rect_tree<63, double, 1, 2> tree63b1d2rT;
    typedef mjr::MR_rect_tree<63, double, 1, 3> tree63b1d3rT;
    typedef mjr::MR_rect_tree<63, double, 1, 4> tree63b1d4rT;
}

#define MJR_INCLUDE_MR_rect_tree
#endif
