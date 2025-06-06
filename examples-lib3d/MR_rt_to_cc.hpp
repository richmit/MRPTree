// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      MR_rt_to_cc.hpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Implimentation for the MR_rt_to_cc class.@EOL
 @keywords  tree cell complex
 @std       C++23
 @see       MR_rect_tree.hpp, MR_cell_cplx.hpp
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
#ifndef MJR_INCLUDE_MR_rt_to_cc

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <tuple>                                                         /* STL tuples              C++11    */
#include <vector>                                                        /* STL vector              C++11    */
#include <string>                                                        /* C++ strings             C++11    */
#include <variant>                                                       /* C++ variant type        C++17    */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put everything in the mjr namespace
namespace mjr {
  /** @brief Tessellates a MR_rect_tree object, and places the result into an MR_cell_cplx object.

      From a structural perspective this class simply a templated collection of types and static methods all designed to work with pairs of MR_rect_tree and
      MR_cell_cplx objects.  From a functional, or pattern, perspective this class might classified as a MR_cell_cplx pseudo-constructor or proto-factory. I
      think of it as a bridge between MR_rect_tree and MR_cell_cplx objects.  A collection of helper types that ease working with pairs of MR_rect_tree and
      MR_cell_cplx objects containing data derived from the same sampled function.

      @tparam rt_t The type of supported MR_rect_tree objects
      @tparam cc_t The type of supported MR_cell_cplx objects */
  template <class rt_t, class cc_t>
  requires (std::is_same<typename rt_t::src_t, typename cc_t::uft_t>::value)
  class MR_rt_to_cc {

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Types imported from MR_rect_tree and MR_cell_cplx. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      typedef typename cc_t::node_data_idx_lst_t cc_node_data_idx_lst_t;
      typedef typename cc_t::node_data_t         cc_node_data_t;
      typedef typename cc_t::node_idx_list_t     cc_node_idx_list_t;
      typedef typename cc_t::node_idx_t          cc_node_idx_t;
      typedef typename cc_t::cell_verts_t              cc_cell_verts_t;
      typedef typename cc_t::uft_t               cc_uft_t;
      typedef typename rt_t::diti_list_t         rt_diti_list_t;
      typedef typename rt_t::diti_t              rt_diti_t;
      typedef typename rt_t::drpt_t              rt_drpt_t;
      typedef typename rt_t::drpt2real_func_t    rt_drpt2real_func_t;
      typedef typename rt_t::rrpt_t              rt_rrpt_t;
      typedef typename rt_t::drpt2rrpt_func_t    rt_drpt2rrpt_func_t;
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Describe point source. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Specify a source space for a data index. */
      enum class val_src_spc_t { FDOMAIN,    //!< Sample function domain space.
                                 FRANGE,     //!< Sample function range space.
                                 CONSTANT   //!< A pseudo-source that returns a constant.
                               };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold an integer or float. */
      typedef std::variant<int, cc_uft_t> iorf_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type used to hold a description of how to extract a scalar value from a tree object */
      typedef std::tuple<val_src_spc_t, iorf_t> val_src_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** A list of val_src_t objects.  */
      typedef std::vector<val_src_t> val_src_lst_t;
      //@}

    private:
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Utility Functions. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Create approprate inputs for ccplx.create_dataset_to_point_mapping() from tree data.  
          @param rtree  The MR_rect_tree with source data
          @param ccplx  The MR_cell_cplx to populate wiht a pont mapping
          @param rt_dil Description of data sources */
      inline static void create_dataset_to_point_mapping(const rt_t& rtree, cc_t& ccplx, const val_src_lst_t& rt_dil) {
        cc_node_data_idx_lst_t cc_data_idx_lst(3);
        for(int i=0; i<3; ++i)
          if(get<0>(rt_dil[i]) == val_src_spc_t::FDOMAIN)
            cc_data_idx_lst[i] = get<int>(get<1>(rt_dil[i]));
          else if(get<0>(rt_dil[i]) == val_src_spc_t::FRANGE)
            cc_data_idx_lst[i] = get<int>(get<1>(rt_dil[i])) + rtree.domain_dimension;
          else if(get<0>(rt_dil[i]) == val_src_spc_t::CONSTANT)
            cc_data_idx_lst[i] = get<cc_uft_t>(get<1>(rt_dil[i]));
        ccplx.create_dataset_to_point_mapping(cc_data_idx_lst);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Given rt coordinates, extract point/scalar/vector data, and add point/data to cc
          @param ccplx  The MR_cell_cplx to populate with geometry
          @param rtree  The MR_rect_tree with source data
          @param diti   The point coordinate in rtree */
      inline static cc_node_idx_t add_node(cc_t& ccplx, const rt_t& rtree, rt_diti_t diti) {
        return add_node(ccplx, rtree.diti_to_drpt(diti), rtree.get_sample(diti));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Given rt coordinates, extract point/scalar/vector data, and add point/data to cc
          @param ccplx    The MR_cell_cplx to populate with geometry
          @param dom_pnt  Domain point
          @param rng_pnt  Range point */
      inline static cc_node_idx_t add_node(cc_t& ccplx, rt_drpt_t dom_pnt, rt_rrpt_t rng_pnt) {
        return ccplx.add_node(rt_pnt_to_cc_pnt(dom_pnt, rng_pnt));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Given rt_t domain & range points, produce a cc_t point data vector.
          @param dom_pnt  Domain point
          @param rng_pnt  Range point */
      inline static cc_node_data_t rt_pnt_to_cc_pnt(rt_drpt_t dom_pnt, rt_rrpt_t rng_pnt)  {
        cc_node_data_t pd;
        if constexpr (rt_t::domain_dimension == 1)
          pd.push_back(dom_pnt);
        else
          for(auto v: dom_pnt)
            pd.push_back(v);
        if constexpr (rt_t::range_dimension == 1)
          pd.push_back(rng_pnt);
        else
          for(auto v: rng_pnt)
            pd.push_back(v);
        return pd;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Transform a MR_cell_cplx::node_data_t value into MR_rect_tree::drpt_t.  
          @param pd The point data to convert.*/
      inline static rt_drpt_t node_data_to_drpt(const cc_node_data_t& pd) {
        rt_drpt_t ret;
        if constexpr (rt_t::domain_dimension == 1) {
          ret = pd[0];
        } else {
          for(int i=0; i<rt_t::domain_dimension; ++i)
            ret[i] = pd[i];
        }
        return ret;
      }
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Utility Functions. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a MR_rect_tree range index into an index for a point data array
          @param tree_range_index value to convert */
      inline static int rt_rng_idx_to_pd_idx(int tree_range_index) {
        return (tree_range_index + rt_t::domain_dimension);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Given an edge with one good point and one NaN point, find the longest segment from the good point toward the NaN point.

          Note we normally use this function when we detect a NaN in a geometric point (i.e. the things with a node_idx_t).  This solver solves until
          the return of func has no NaNs.  Those two criteria might not be the same thing, but it's OK.

          @param ccplx                   The MR_cell_cplx to populate with geometry
          @param rtree                   The MR_rect_tree with source data
          @param good_point_ccplx_index  Good point index in the ccplx object
          @param good_point_rtree_index  Good point index in the rtree object
          @param sick_point_rtree_index  Bad point index in the rtree object
          @param func                    The function to use for the solver
          @param solver_epsilon          Used as a distance threshold between sick point and solved endpoint in the tree domain space. */
      static cc_node_idx_t nan_edge_solver(cc_t&               ccplx,
                                          const rt_t&         rtree,
                                          cc_node_idx_t        good_point_ccplx_index,
                                          rt_diti_t           good_point_rtree_index,
                                          rt_diti_t           sick_point_rtree_index,
                                          rt_drpt2rrpt_func_t func,
                                          cc_uft_t            solver_epsilon=cc_t::epsilon/100
                                         ) {
        // Solver cache.  Clear it if we have a different rtree object from last time.
        static std::unordered_map<rt_diti_t, std::unordered_map<rt_diti_t, cc_node_idx_t>> nan_solver_cache;
        static const rt_t* rtree_cache = nullptr;
        if (rtree_cache != &rtree) {
          nan_solver_cache.clear();
          rtree_cache = &rtree;
        }
        // Check to see if we solved this one before
        if (nan_solver_cache.contains(sick_point_rtree_index))
          if (nan_solver_cache[sick_point_rtree_index].contains(good_point_rtree_index))
            return  nan_solver_cache[sick_point_rtree_index][good_point_rtree_index];
        // Apparently we need to solve this one as it's not in the case
        rt_drpt_t good_point_drpt = rtree.diti_to_drpt(good_point_rtree_index);
        rt_drpt_t sick_point_drpt = rtree.diti_to_drpt(sick_point_rtree_index);
        rt_rrpt_t good_point_rrpt = rtree.get_sample(good_point_rtree_index);
        rt_drpt_t init_point_drpt = good_point_drpt;
        while ( (rtree.drpt_distance_inf(good_point_drpt, sick_point_drpt) > solver_epsilon) ) {
          rt_drpt_t md_point_drpt = rtree.drpt_midpoint(good_point_drpt, sick_point_drpt);
          rt_rrpt_t y = func(md_point_drpt);
          if (rtree.rrpt_is_nan(y)) {
            sick_point_drpt = md_point_drpt;
          } else {
            good_point_drpt = md_point_drpt;
            good_point_rrpt = y;
          }
        }
        // Figure out what to return, add it to the cache, and return.
        cc_node_idx_t ret;
        if (rtree.drpt_distance_inf(good_point_drpt, init_point_drpt) < (ccplx.epsilon)) // Use ccplx here!!!
          ret = good_point_ccplx_index;
        else
          ret = add_node(ccplx, good_point_drpt, good_point_rrpt);
        nan_solver_cache[sick_point_rtree_index][good_point_rtree_index] = ret;
        return ret;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Poly data construction */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Populate attached MR_cell_cplx object from data in attached MR_rect_tree object.

          construct_geometry_fans(), unlike the other geometry construction methods, is capable of "healing" some broken edges -- edges with one good point
          and one NaN point.  It uses the origional sampling function to solve along the edge to produce the longest non-NaN edge possible.  Then it
          constructs cells using the new piont(s).  If func is nullptr, then no edge healing is preformed.  This feature works for segments & triangles only.

          @verbatim
          | Geom       | Dom Dim | Out Dim | Result             |
          |------------+---------+---------+--------------------|
          | FANS       |       2 |       1 | Triangle Edges     |
          | FANS       |       2 |       2 | Triangles          |
          | FANS       |       3 |       1 | Pyramid Edges      |
          | FANS       |       3 |       2 | Pyramid Faces      |
          | FANS       |       3 |       3 | Solid Pyramids     |
          @endverbatim

          @param ccplx             The MR_cell_cplx to populate with geometry
          @param rtree             The MR_rect_tree with source data
          @param cells             List of cells to output from rtree
          @param output_dimension  Parts of cells to output
          @param point_src         Point sources
          @param func              The function was used to sample the tree */
      static int construct_geometry_fans(cc_t&               ccplx,
                                         const rt_t&         rtree,
                                         rt_diti_list_t      cells,
                                         int                 output_dimension,
                                         val_src_lst_t       point_src,
                                         rt_drpt2rrpt_func_t func = nullptr
                                        ) {
        create_dataset_to_point_mapping(rtree, ccplx, point_src);
        if (rtree.domain_dimension == 1) {
          for(auto& cell: cells) {
            cc_node_idx_t ctr_pnti = add_node(ccplx, rtree, cell);
            rt_diti_list_t corners = rtree.ccc_get_corners(cell);
            cc_node_idx_t cn0_pnti = add_node(ccplx, rtree, corners[0]);
            cc_node_idx_t cn1_pnti = add_node(ccplx, rtree, corners[1]);
            if (func) { // We have a func, so we can "heal" broken edges.
              if (ctr_pnti < 0) { // Center: Broken. Left:
                if(cn0_pnti >= 0) { // Center: Broken.  Left: Good.
                  cc_node_idx_t np = nan_edge_solver(ccplx, rtree, cn0_pnti, corners[0], cell, func);
                  ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {cn0_pnti, np}, output_dimension);
                }
                if(cn1_pnti >= 0) { // Center: Broken.  Right: Good.
                  cc_node_idx_t np = nan_edge_solver(ccplx, rtree, cn1_pnti, corners[1], cell, func);
                  ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {np, cn1_pnti}, output_dimension);
                }
              } else {             // Center: Good.
                if(cn0_pnti < 0) { // Center: Good.  Left: Broken.
                  cc_node_idx_t np = nan_edge_solver(ccplx, rtree, ctr_pnti, cell, corners[0], func);
                  ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {np, ctr_pnti}, output_dimension);
                } else {           // Center: Good.  Left: Good.
                  ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {cn0_pnti, ctr_pnti}, output_dimension);
                }
                if(cn1_pnti < 0) { // Center: Good.  Right: Broken.
                  cc_node_idx_t np = nan_edge_solver(ccplx, rtree, ctr_pnti, cell, corners[1], func);
                  ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {ctr_pnti, np}, output_dimension);
                } else {           // Center: Good.  Left: Good.
                  ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {ctr_pnti, cn1_pnti}, output_dimension);
                }
              }
            } else {
              ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {cn0_pnti, ctr_pnti}, output_dimension);
              ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {ctr_pnti, cn1_pnti}, output_dimension);
            }
          }
        } else if (rtree.domain_dimension == 2) {
          for(auto& cell: cells) {
            if (func) { // We have a func, so we can "heal" broken edges.
              for(int i=0; i<2; i++) {
                for(int j=-1; j<2; j+=2) {
                  std::vector<rt_diti_list_t> triangles;
                  rt_diti_list_t nbrs = rtree.get_existing_neighbor(cell, i, j);
                  if (nbrs.size() > 1) {
                    for(auto n: nbrs) {
                      rt_diti_list_t corners = rtree.ccc_get_corners(n, i, -j);
                      if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                        triangles.push_back({corners[1], corners[0], cell});
                      else
                        triangles.push_back({corners[0], corners[1], cell});
                    }
                  } else {
                    rt_diti_list_t corners = rtree.ccc_get_corners(cell, i, j);
                    if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                      triangles.push_back({corners[1], corners[0], cell});
                    else
                      triangles.push_back({corners[0], corners[1], cell});
                  }
                  for(auto triangle: triangles) {
                    std::array<cc_node_idx_t, 3> tpnts {add_node(ccplx, rtree, triangle[0]),
                                                       add_node(ccplx, rtree, triangle[1]),
                                                       add_node(ccplx, rtree, triangle[2])};
                    int num_bad = static_cast<int>(std::count_if(tpnts.begin(), tpnts.end(), [](cc_node_idx_t i) { return i<0; }));
                    if (num_bad == 0) {
                      ccplx.add_cell(cc_t::cell_kind_t::TRIANGLE, {tpnts[0], tpnts[1], tpnts[2]}, output_dimension);
                    } else if ((num_bad == 1) || (num_bad == 2)) {
                      // Rotate points so we only have two cases to think about...
                      std::array<int, 3> p {0, 1, 2};
                      if ( ((tpnts[1] < 0) && (num_bad == 1)) || ((tpnts[1] >= 0) && (num_bad == 2)) )
                        p = {1, 2, 0};
                      else if ( ((tpnts[2] < 0) && (num_bad == 1)) || ((tpnts[2] >= 0) && (num_bad == 2)) )
                        p = {2, 0, 1};
                      // Solve for edge 0-1 & 0-2
                      if (num_bad == 1) {
                        cc_node_idx_t np1 = nan_edge_solver(ccplx, rtree, tpnts[p[1]], triangle[p[1]], triangle[p[0]], func);
                        cc_node_idx_t np2 = nan_edge_solver(ccplx, rtree, tpnts[p[2]], triangle[p[2]], triangle[p[0]], func);
                        ccplx.add_cell(cc_t::cell_kind_t::TRIANGLE, {np1, tpnts[p[1]], tpnts[p[2]]}, output_dimension);
                        ccplx.add_cell(cc_t::cell_kind_t::TRIANGLE, {tpnts[p[2]], np2, np1}, output_dimension);
                      } else {
                        cc_node_idx_t np1 = nan_edge_solver(ccplx, rtree, tpnts[p[0]], triangle[p[0]], triangle[p[1]], func);
                        cc_node_idx_t np2 = nan_edge_solver(ccplx, rtree, tpnts[p[0]], triangle[p[0]], triangle[p[2]], func);
                        ccplx.add_cell(cc_t::cell_kind_t::TRIANGLE, {tpnts[p[0]], np1, np2}, output_dimension);
                      }
                    }
                  }
                }
              }
            } else { // We don't have a func, so we can can't "heal" broken edges.  This is much faster. ;)
              cc_node_idx_t ctr_pnti = add_node(ccplx, rtree, cell);
              if (ctr_pnti >= 0) { // Center point was good, let's try and make some triangles...
                for(int i=0; i<2; i++) {
                  for(int j=-1; j<2; j+=2) {
                    rt_diti_list_t nbrs = rtree.get_existing_neighbor(cell, i, j);
                    if (nbrs.size() > 1) {
                      for(auto n: nbrs) {
                        rt_diti_list_t corners = rtree.ccc_get_corners(n, i, -j);
                        cc_node_idx_t cn0_pnti = add_node(ccplx, rtree, corners[0]);
                        cc_node_idx_t cn1_pnti = add_node(ccplx, rtree, corners[1]);
                        if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                          std::swap(cn0_pnti, cn1_pnti);
                        ccplx.add_cell(cc_t::cell_kind_t::TRIANGLE, {cn0_pnti, cn1_pnti, ctr_pnti}, output_dimension);
                      }
                    } else {
                      rt_diti_list_t corners = rtree.ccc_get_corners(cell, i, j);
                      cc_node_idx_t cn0_pnti = add_node(ccplx, rtree, corners[0]);
                      cc_node_idx_t cn1_pnti = add_node(ccplx, rtree, corners[1]);
                      if( ((i == 0) && (j == -1)) || ((i == 1) && (j == 1)) )
                        std::swap(cn0_pnti, cn1_pnti);
                      ccplx.add_cell(cc_t::cell_kind_t::TRIANGLE, {cn0_pnti, cn1_pnti, ctr_pnti}, output_dimension);
                    }
                  }
                }
              }
            }
          }
        } else if (rtree.domain_dimension == 3) {
          for(auto& cell: cells) {
            cc_node_idx_list_t new_cell(5);
            new_cell[4] = add_node(ccplx, rtree, cell);
            std::array<int, 5> p {0, 1, 3, 2, 4};
            if (new_cell[4] >= 0) { // Center point was good, let's try and make some pyramids...
              for(int dim=0; dim<3; dim++) {
                for(int dir=-1; dir<2; dir+=2) {
                  rt_diti_list_t nbrs = rtree.get_existing_neighbor(cell, dim, dir);
                  if (nbrs.size() > 1) {
                    for(auto n: nbrs) {
                      rt_diti_list_t corners = rtree.ccc_get_corners(n, dim, -dir);
                      for(int k=0; k<4; ++k)
                        new_cell[p[k]] = add_node(ccplx, rtree, corners[k]);
                      ccplx.add_cell(cc_t::cell_kind_t::PYRAMID, new_cell, output_dimension);
                    }
                  } else {
                    rt_diti_list_t corners = rtree.ccc_get_corners(cell, dim, dir);
                    for(int k=0; k<4; ++k)
                      new_cell[p[k]] = add_node(ccplx, rtree, corners[k]);
                    ccplx.add_cell(cc_t::cell_kind_t::PYRAMID, new_cell, output_dimension);
                  }
                }
              }
            }
          }
        } else { // if (rtree.domain_dimension > 3) {
          std::cout << "ERROR: construct_geometry_fans: output_dimension>3 not supported for output_dimension>0!" << std::endl;
          return 1;
        }
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      static int construct_geometry_fans(cc_t&               ccplx,
                                         const rt_t&         rtree,
                                         int                 output_dimension,
                                         val_src_lst_t       point_src,
                                         rt_drpt2rrpt_func_t func = nullptr
                                        ) {
        return construct_geometry_fans(ccplx, rtree, rtree.get_leaf_cells(), output_dimension, point_src, func);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Populate attached MR_cell_cplx object from data in attached MR_rect_tree object.

          Only 0D vertex cells are produced.  While similar results may be obtained by setting the output_dimension to zero and calling
          construct_geometry_fan() or construct_geometry_rects(), this method is much faster.  This method also provides the option of only outputting centers without
          corners.

          @verbatim
          | Geom   | output_centers | output_corners | Result        | output_dimension=0 equivalent |
          |--------+----------------+----------------+---------------+-------------------------------|
          | POINTS | true           | true           | cell vertexes | construct_geometry_fans       |
          | POINTS | true           | false          | cell centers  |                               |
          | POINTS | false          | true           | cell corners  | construct_geometry_rects      |
          | POINTS | false          | false          | No Points     |                               |
          @endverbatim

          @param ccplx           The MR_cell_cplx to populate with geometry
          @param rtree           The MR_rect_tree with source data
          @param cells           List of tree cells from which to construct geometry
          @param point_src       Point sources
          @param output_centers  Create vertexes for cell  centers
          @param output_corners  Create vertexes for cell corners*/
      static int construct_geometry_points(cc_t&          ccplx,
                                           const rt_t&    rtree,
                                           rt_diti_list_t cells,
                                           val_src_lst_t  point_src,
                                           bool           output_centers,
                                           bool           output_corners
                                          ) {
        create_dataset_to_point_mapping(rtree, ccplx, point_src);
        if (output_centers && output_corners) {
          for(auto& cell: cells)
            for(auto& vert: rtree.ccc_get_vertexes(cell))
              ccplx.add_cell(cc_t::cell_kind_t::POINT, {add_node(ccplx, rtree, vert)});
        } else if (output_centers) {
          for(auto& cell: cells)
            ccplx.add_cell(cc_t::cell_kind_t::POINT, {add_node(ccplx, rtree, cell)});
        } else if (output_corners) {
          for(auto& cell: cells)
            for(auto& vert: rtree.ccc_get_corners(cell))
              ccplx.add_cell(cc_t::cell_kind_t::POINT, {add_node(ccplx, rtree, vert)});
        } else {
          std::cout << "WARNING: construct_geometry_points: Both output_centers & output_corners are FALSE.  No geometry created!" << std::endl;
          return 1;
        }
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      static int construct_geometry_points(cc_t&         ccplx,
                                           const rt_t&   rtree,
                                           val_src_lst_t point_src,
                                           bool          output_centers,
                                           bool          output_corners
                                          ) {
        return construct_geometry_points(ccplx, rtree, rtree.get_leaf_cells(), point_src, output_centers, output_corners);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Populate a MR_cell_cplx object from data in a MR_rect_tree object.

          The resulting geometric structure in the MR_cell_cplx object will consist of 'rectangular' cell types (points, segments, rectangles, & hexahedra).

          @verbatim
          | Geom       | Dom Dim | Out Dim | Result             |
          |------------+---------+---------+--------------------|
          | RECTANGLES |     1-3 |       0 | Cell Corner Points |
          | RECTANGLES |     2-3 |       1 | Cell Edges         |
          | RECTANGLES |       2 |       2 | 2D Rectangles      |
          | RECTANGLES |       3 |       2 | Cell Faces         |
          | RECTANGLES |       3 |       3 | Solid Hexahedra    |
          @endverbatim

          @param ccplx                The MR_cell_cplx to populate with geometry
          @param rtree                The MR_rect_tree with source data
          @param cells                List of tree cells from which to construct geometry
          @param output_dimension     Parts of cells to output
          @param point_src            Point sources
          @param degenerate_fallback  If the rectangle is degenerate, try and make a triangle. (only works for cc_t::cell_kind_t::QUAD) */
      static int construct_geometry_rects(cc_t&          ccplx,
                                          const rt_t&    rtree,
                                          rt_diti_list_t cells,
                                          int            output_dimension,
                                          val_src_lst_t  point_src,
                                          bool           degenerate_fallback = true
                                         ) {
        create_dataset_to_point_mapping(rtree, ccplx, point_src);
        for(auto& cell: cells) {
          std::vector<cc_node_idx_t> cnr_pti;
          rt_diti_list_t corners = rtree.ccc_get_corners(cell);
          for(auto& corner: corners) {
            cc_node_idx_t pnti = add_node(ccplx, rtree, corner);
            cnr_pti.push_back(pnti);
          }
          if (rtree.domain_dimension == 1) {
            ccplx.add_cell(cc_t::cell_kind_t::SEGMENT, {cnr_pti[0], cnr_pti[1]}, output_dimension);
          } else if (rtree.domain_dimension == 2) {
            const std::array<int, 4> p = {0, 1, 3, 2};
            bool try_harder = !(ccplx.add_cell(cc_t::cell_kind_t::QUAD, {cnr_pti[0], cnr_pti[1], cnr_pti[3], cnr_pti[2]}, output_dimension));
            if ( degenerate_fallback && try_harder) { // Try for a triangle if we have a NaN point or an adjacent pair of duplicate points
              for(int i=0; i<4; i++) {
                if ((cnr_pti[p[i]] < 0) || (cnr_pti[p[(i+0)%4]] == cnr_pti[p[(i+1)%4]])) {
                  ccplx.add_cell(cc_t::cell_kind_t::TRIANGLE, {cnr_pti[p[(i+1)%4]], cnr_pti[p[(i+2)%4]], cnr_pti[p[(i+3)%4]]}, output_dimension);
                  break;
                }
              }
            }
          } else { // if(rtree.domain_dimension == 3) {
            ccplx.add_cell(cc_t::cell_kind_t::HEXAHEDRON,
                           {cnr_pti[0], cnr_pti[1], cnr_pti[3], cnr_pti[2],
                            cnr_pti[4], cnr_pti[5], cnr_pti[7], cnr_pti[6]},
                           output_dimension);
          }
        }
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload */
      static int construct_geometry_rects(cc_t&         ccplx,
                                          const rt_t&   rtree,
                                          int           output_dimension,
                                          val_src_lst_t point_src,
                                          bool          degenerate_fallback = true
                                         ) {
        return construct_geometry_rects(ccplx, rtree, rtree.get_leaf_cells(), output_dimension, point_src, degenerate_fallback);
      }
      //@}


      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Function Adapters */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Adapt a MR_rect_tree::drpt2rrpt_func_t (sample function) to a MR_cell_cplx::p2data_func_t (Point Data Transform).
          @param func The function to adapt
          @param pd   Point data to be passed to func. */
      inline static cc_node_data_t tsampf_to_cdatf(rt_drpt2rrpt_func_t func,
                                                  cc_node_data_t       pd) {
        rt_drpt_t xpt = node_data_to_drpt(pd);
        return rt_pnt_to_cc_pnt(xpt, func(xpt));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Adapt a MR_rect_tree::drpt2real_func_t (Domain Point SDF) to MR_cell_cplx::p2real_func_t (Point Data SDF).
          @param func        The function to adapt
          @param pd          Point data to be passed to func. */
      inline static cc_uft_t tsdf_to_csdf(rt_drpt2real_func_t func,
                                          cc_node_data_t       pd) {
        return static_cast<cc_uft_t>(func(node_data_to_drpt(pd)));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Evalate an MR_cell_cplx::p2real_func_t (Point Data SDF) function from a MR_rect_tree::drpt2rrpt_func_t (Sample Function) and level data.
          @param range_index Index into range of origional sample function
          @param level       Level to check range element aginst
          @param func        The function to adapt
          @param pd          Point data to be passed to func. */
      inline static cc_uft_t tsampf_to_clcdf(int                 range_index,
                                             cc_uft_t            level,
                                             rt_drpt2rrpt_func_t func,
                                             cc_node_data_t       pd) {
        if constexpr (rt_t::domain_dimension == 1) {
          return static_cast<cc_uft_t>(func(node_data_to_drpt(pd))) - level;
        } else {
          return static_cast<cc_uft_t>(func(node_data_to_drpt(pd))[range_index]) - level;
        }
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Mathematical Tools */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Drop cells from a MR_cell_cplx object using an MR_rect_tree domain SDF function.
       */
      inline static int cull_cc_cells_on_domain_sdf_boundry(cc_t&               ccplx,
                                                            rt_drpt2real_func_t sdf_func) {
        return ccplx.cull_cells([&ccplx, &sdf_func](cc_cell_verts_t c) { return ccplx.cell_near_sdf_boundry(c, [&sdf_func](cc_node_data_t pd) { return (tsdf_to_csdf(sdf_func, pd)); }); });
      }


      //@}
  };
}

#define MJR_INCLUDE_MR_rt_to_cc
#endif
