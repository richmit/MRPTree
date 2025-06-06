// -*- Mode:C++; Coding:us-ascii-unix; fill-column:158 -*-
/*******************************************************************************************************************************************************.H.S.**/
/**
 @file      MR_cell_cplx.hpp
 @author    Mitch Richling http://www.mitchr.me/
 @date      2024-07-13
 @brief     Implementation for the MR_cell_cplx class.@EOL
 @keywords  VTK polydata PLY file MR_rect_tree polygon triangulation cell complex tessellation
 @std       C++23
 @see       MR_rect_tree.hpp, MR_rt_to_cc.hpp
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
#ifndef MJR_INCLUDE_MR_cell_cplx

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <algorithm>                                                     /* STL algorithm           C++11    */
#include <array>                                                         /* array template          C++11    */
#include <cmath>                                                         /* std:: C math.h          C++11    */
#include <fstream>                                                       /* C++ fstream             C++98    */
#include <functional>                                                    /* STL funcs               C++98    */
#include <iomanip>                                                       /* C++ stream formatting   C++11    */
#include <iostream>                                                      /* C++ iostream            C++11    */
#include <map>                                                           /* STL map                 C++11    */
#include <set>                                                           /* STL set                 C++98    */
#include <sstream>                                                       /* C++ string stream       C++      */
#include <string>                                                        /* C++ strings             C++11    */
#include <unordered_map>                                                 /* STL hash map            C++11    */
#include <variant>                                                       /* C++ variant type        C++17    */
#include <vector>                                                        /* STL vector              C++11    */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "MRMathCPP.hpp"                                                   /* My Simple Math Utilities          */

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Put everything in the mjr namespace
namespace mjr {
  /** @brief Template class used to store and transform cell complexes (mesh/triangluation/etc...) and data sets generated from  MR_rect_tree sample data.

      The primary use case is to store and manipulate geometric (i.e. mesh) data derived from MR_rect_tree objects.  On the mesh manipulation front the focus
      is on computations that require knowledge related to the function being approximated.

      This code is a quick-n-dirty hack job.  I simply wanted to quickly get to a point where I could visually test MR_rect_tree objects.  Logically this
      class should be split into several classes (point, point list, cell, cell list, 3D vectors, 3D analytic geometry, etc...).  It's a bit of a mess.

      Design Constraints:
        - Geometry
          - Only four types of cells are supported: points, segments, triangles, quads, pyramids, & hexahedrons.
          - No tetrahedrons!  This may seem bazaar; however, they are simply not a natural product of MR_rect_tree tessellation.
          - Unlike MR_rect_tree, this class can't operate on data of dimension higher than 3!
        - Data
          - Only point data is supported (i.e. no cell data)
          - Data elements are scalar and/or vector.  No tensors, fields, etc...
          - uft_t types are used for atomic data components (scalars are uft_t, and components of vectors are uft_t).
        - Computations/Manipulation
          - Supported is focused on manipulation that benefits from knowledge of the approximated function
          - The idea is that one can use one of the many mesh packages for standard types of mesh computations
        - I/O
          - Write only.
          - VTK files
            - Only Unstructured_Grid files are supported
            - Both legacy & XML files are supported, but they are ASCII only.
            - XML files are serial and self contained
            - Note that NaN's in ASCII legacy files are not properly handled by many VTK applications
        - Error Checking
          - It's pretty limited
          - Memory allocation -- you run out, the thing crashes
          - Cells that are a part of other cells may be added (i.e. a segment that is part of an existing triangle)
          - Geometry checks are quite basic, and I don't even bother with things like convexity. ;)

      Cells supported:

       - Point
         - Zero-dimensional cell.
         - Defined by a single point.

       - Segment
         - One-dimensional cell.
         - Defined an ordered list of two points.
         - The direction along the line is from the first point to the second point.

       - Triangle
         - Two-dimensional cell.
         - Defined by an ordered list of three points.
         - The points are ordered counterclockwise, defining a surface normal using the right-hand rule.

       - Quadrilateral
         - Two-dimensional cell.
         - Defined by an ordered list of four points.
         - The points are ordered counterclockwise, defining a surface normal using the right-hand rule.
         - The points must be coplainer
         - The quadrilateral is convex and its edges must not intersect.

       - Hexahedron
         - Three-dimensional cell consisting of six quadrilateral faces, twelve edges, and eight vertices.
         - The hexahedron is defined by an ordered list of eight points.
         - The faces and edges must not intersect any other faces and edges, and the hexahedron must be convex.

       - Pyramid
         - Three-dimensional cell consisting of one quadrilateral face, four triangular faces, eight edges, and five vertices.
         - The pyramid is defined by an ordered list of five points.
         - The four points defining the quadrilateral base plane must be convex
         - The fifth point, the apex, must not be co-planar with the base points.

         \verbatim
                                                                              HEXAHEDRON                            4
                                                                                                       PYRAMID     .^.
                                                                              3---------2                       .  / \ .
                                        TRIANGLE          QUAD               /|        /|                     .   /   \  .
                                                                            / |       / |                  .     /     \   .
          POINT    SEGMENT                 2            3---------2        7---------6  |                 1...../.......\.....0  Back
                                          / \           |         |        |  |      |  |                 |    /         \    |
          0        0-------1             /   \          |         |        |  0------|--1 Back            |   /           \   |
                                        /     \         |         |        | /       | /                  |  /             \  |
                                       0-------1        |         |        |/        |/                   | /               \ |
                                                        0---------1        4---------5 Front              |/                 \|
                                                                                                          2 ----------------- 3  Front
         \endverbatim

      A number of quality checks may be performed on points and cells before they are added to the object.  These checks can slow down execution of add_cell()
      by an order of magnitude, and almost double the RAM required for this class.  These checks are entirely optional.

        - chk_point_unique:
          - Many 3D file formats store the master list of potential points, and then use an integer index into this list when defining geometric objects.
          - Most visualization software is pretty tolerant of having duplicate points on this list when just doing rendering.
          - Duplicate points can break many software packages pretty badly when it comes to computation.
          - Depending on how points are added to this object, this check can avoid quite a bit of wasted RAM and produce *much* smaller output files.
          - I normally leave this check turned on.
        - chk_cell_unique:
          - Most packages can render objects with duplicate cells, but you might see graphical artifacts or glitching.
          - When using a mesh with duplicate cells, computations can sometimes fail or take a *very* long time
          - I normally leave this check turned on.
        - chk_cell_vertexes:
          - This check makes sure that cells have the correct number of vertexes and that those vertexes are valid (on the master point list)
          - This check also makes sure cells have no duplicate vertexes.
          - In general this check is pretty cheap, and catches a great many common mesh problems.
          - This is the only check required for cell_kind_t::POINT and cell_kind_t::SEG cells
        - chk_cell_dimension:
          - This check makes sure that PYRAMIDS & HEXAHEDRONS are not co-planar and that TRIANGLES & QUADS are not co-linear.
          - Combined with chk_cell_vertexes this check is all that is required for cell_kind_t::TRIANGLE cells
        - chk_cell_edges:
          - This check makes sure that no cell edges intersect in disallowed ways
          - This check when combined with the above checks, usually produce cells good enough for most software packages.
            - cell_kind_t::QUAD cells might be concave or non-plainer.
            - cell_kind_t::HEXAHEDRON cells might be concave or degenerate (bad edge-face intersections)
            - cell_kind_t::PYRAMID cells might be concave.

      Levels Of Cell Goodness

       \verbatim
       +-------------+-------------------+--------------------+--------------------+-------------------------------+-------------------------------+
       |             | chk_point_unique  | chk_point_unique   | chk_point_unique   | chk_point_unique              | chk_point_unique              |
       |             | chk_cell_unique   | chk_cell_unique    | chk_cell_unique    | chk_cell_unique               | chk_cell_unique               |
       |             | chk_cell_vertexes | chk_cell_vertexes  | chk_cell_vertexes  | chk_cell_vertexes             | chk_cell_vertexes             |
       |             |                   | chk_cell_dimension | chk_cell_dimension | chk_cell_dimension            | chk_cell_dimension            |
       |             |                   |                    | chk_cell_edges     | chk_cell_edges                | chk_cell_edges                |
       |             |                   |                    |                    | check_cell_face_intersections | check_cell_face_intersections |
       |             |                   |                    |                    |                               | check_cell_faces_plainer      |
       |             |                   |                    |                    |                               | check_cell_convex             |
       +-------------+-------------------+--------------------+--------------------+-------------------------------+-------------------------------+
       | vertex      | Perfect           |                    |                    |                               |                               |
       | segment     | Perfect           |                    |                    |                               |                               |
       | triangle    |                   | Perfect            |                    |                               |                               |
       | quad        |                   |                    | Good Enough        |                               | Perfect                       |
       | hexahedron  |                   |                    | Mostly OK          | Good Enough                   | Perfect                       |
       | pyramid     |                   |                    | Good Enough        |                               | Perfect                       |
       +-------------+-------------------+--------------------+--------------------+-------------------------------+-------------------------------+
       \endverbatim


        @tparam chk_point_unique   Only allow unique points on the master point list
        @tparam chk_cell_unique    Only allow unique cells on the master cell list.  Won't prevent cells that are sub-cells of existing cells.
        @tparam chk_cell_vertexes  Do cell vertex checks (See: cell_stat_t)
        @tparam chk_cell_dimension Do cell dimension checks (See: cell_stat_t)
        @tparam chk_cell_edges     Do cell edge checks (See: cell_stat_t)
        @tparam flt_type           Type for floating point values (
        @tparam eps                Epsilon used to detect zero */
  template <bool     chk_point_unique,
            bool     chk_cell_unique,
            bool     chk_cell_vertexes,
            bool     chk_cell_dimension,
            bool     chk_cell_edges,
            typename flt_type,
            flt_type eps
            // double eps
            >
  requires ((std::is_floating_point<flt_type>::value) &&
            (eps > std::numeric_limits<flt_type>::epsilon()))
  class MR_cell_cplx {

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Global Types & Constants. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Universal floating point type -- used for point components, scalar data values, vector components, all computation, etc.... */
      typedef flt_type uft_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** epsilon */
      constexpr static uft_t epsilon      = eps;
      constexpr static uft_t uft_epsilon  = std::numeric_limits<uft_t>::epsilon();
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Node Data.

          Cell complex objects have a central repository of nodes.  Each node has a simple data set consisting of a vector of uft_t values.  This data may be
          used for various purposes; however, three of the values in the vector are always required to specify the node's location in R^3 (i.e. the node's
          "point"). Other uses include normal vectors and RGB colors.  In addition subsets of these values may be named, and extracted by name.  Lists of
          indexes are used to define cells. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a data for a single node. */
      typedef std::vector<uft_t> node_data_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type used to hold an index into a node_data_t or a constant MR_cell_cplx::uft_t. */
      typedef std::variant<int, uft_t> node_data_idx_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a list of data indexes -- one element for scalars, 3 for vectors, etc.... */
      typedef std::vector<node_data_idx_t> node_data_idx_lst_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a list values extracted from node data via a MR_cell_cplx::node_data_idx_lst_t (also see MR_cell_cplx::fvec3_t). */
      typedef node_data_t vdat_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** How to construct a point from a node dataset */
      node_data_idx_lst_t node_data_to_pnt;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Integer type used to indentify/index nodes. */
      typedef int node_idx_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** list of node indexes. */
      typedef std::vector<node_idx_t> node_idx_list_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Hold a tuple of real values defining a point/vector in R^3 (also see MR_cell_cplx::vdat_t).  */
      typedef std::array<uft_t, 3> fvec3_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Hold a list of MR_cell_cplx::fvec3_t objects. */
      typedef std::vector<fvec3_t> fvec3_list_t;
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Named Node Datasets.

          node data may be grouped into named datasets, and written to geometry files. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a dataset name */
      typedef std::string named_data_name_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a list of dataset names */
      typedef std::vector<named_data_name_t> named_data_name_list_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to map names to named data sets (index lists). */
      typedef std::map<named_data_name_t, node_data_idx_lst_t> data_name_to_node_data_idx_lst_t;
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cells. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Cell Status.  */
      enum class cell_stat_t { GOOD,           //!< Looks like a good cell                                       N/A
                               TOO_FEW_PNT,    //!< List of points was empty                                     check_cell_vertexes
                               TOO_MANY_PNT,   //!< List of points was too long (>5)                             check_cell_vertexes
                               NEG_PNT_IDX,    //!< Negative point index (i.e. not on the central node data)     check_cell_vertexes
                               BIG_PNT_IDX,    //!< Point index was too big (i.e. not on the central node data)  check_cell_vertexes
                               DUP_PNT,        //!< At least two points have the same index                      check_cell_vertexes
                               DIM_LOW,        //!< Dimension low (degenerate)                                   check_cell_dimension
                               BAD_EDGEI,      //!< Bad edge-edge intersection                                   check_cell_edge_intersections
                               BAD_FACEI,      //!< Bad face-edge intersection                                   check_cell_face_intersections
                               FACE_BENT,      //!< A face (QUAD, HEXAHEDRON, or PYRAMID) was not plainer        check_cell_faces_plainer
                               CONCAVE,        //!< (QUAD, HEXAHEDRON, or PYRAMID) was concave                   TBD
                             };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Cell Types. */
      enum class cell_kind_t { POINT,
                               SEGMENT,
                               TRIANGLE,
                               QUAD,
                               // TETRAHEDRON
                               HEXAHEDRON,
                               PYRAMID,
                             };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold the vertexes of a poly cell */
      typedef node_idx_list_t cell_verts_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a poly cell -- a list of point indexes */
      // struct cell_t { 
      //     cell_kind_t  type;
      //     cell_verts_t verts;
      // };
      //@}

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Node Data. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** The index of the last point added via the add_node() method. */
      node_idx_t last_point_idx = -1;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** True if the last point given to the add_node() method was new -- i.e. not on the central node data.
          Only updated if chk_point_unique is true.  See: last_point_added_was_new() */
      bool last_point_new = true;
      /** @cond */  // Leaving this out of doxygen because I don't like how it shows up...
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Less operator for pnt_lst_uniq_t.
          @return If a & b are close in space, then return false.  Otherwise uses lexicographic ordering. */
      struct pnt_less {
          bool operator()(const fvec3_t& a, const fvec3_t& b) const { return (((std::abs(a[0]-b[0]) > eps) ||
                                                                               (std::abs(a[1]-b[1]) > eps) ||
                                                                               (std::abs(a[2]-b[2]) > eps)) &&
                                                                              ((a[0] < b[0]) ||
                                                                               ((a[0] == b[0]) &&
                                                                                ((a[1] < b[1]) ||
                                                                                 ((a[1] == b[1]) && (a[2] < b[2])))))); }
      };
      /** @endcond */
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Used to uniquify points and assign point index values */
      typedef std::map<fvec3_t, node_idx_t, pnt_less> pnt_to_node_idx_map_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Maps points to point index -- used to detect points with nearly identical geometric points in R^3 */
      pnt_to_node_idx_map_t pnt_to_node_idx_map;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type to hold a all node data. */
      typedef std::vector<node_data_t> node_idx_to_node_data_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Hold all node data sets */
      node_idx_to_node_data_t node_idx_to_node_data;
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Named Node Datasets. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Hold named descriptions of data (scalars & vectors) */
      data_name_to_node_data_idx_lst_t data_name_to_data_idx_lst;
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cells. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type for a list of poly cells */
      typedef std::vector<cell_verts_t> cell_lst_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** House all poly cells */
      cell_lst_t cell_lst;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Used to uniquify cells */
      typedef std::set<cell_verts_t> uniq_cell_lst_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Unique cell list. */
      uniq_cell_lst_t uniq_cell_lst;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** True if the last cell given to the add_cell() method was new -- i.e. not on the master cell list.
          Only updated if chk_cell_unique is true.  Value is invalid if last_cell_stat is NOT cell_stat_t::GOOD.
          See: last_cell_added_was_new() */
      bool last_cell_new = true;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Status of the last cell added via add_cell() method.
          Only updated if (chk_cell_vertexes || chk_cell_dimension | chk_cell_edges) is true. See: status_of_last_cell_added() */
      cell_stat_t last_cell_stat = cell_stat_t::GOOD;
      //@}

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Node Data. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the number of points in the central node data list.
          Note the return type is node_idx_t (a signed integer type) and not a size_t. */
      inline node_idx_t node_count() const {
        return static_cast<node_idx_t>(node_idx_to_node_data.size());  // Yes.  We mean node_idx_to_node_data.
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return point given a node index */
      inline fvec3_t get_pnt(node_idx_t pnt_idx) const { 
        return fvec3_from_node_data(node_data_to_pnt, node_idx_to_node_data[pnt_idx]); 
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Retruns the index of the last point given to the add_node() method. */
      inline node_idx_t idx_of_last_point_added() const {
        return last_point_idx;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Retruns true if the last point given to the add_node() method was a new point. */
      inline node_idx_t last_point_added_was_new() const {
        return last_point_new;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add a point to the central node data store.

          The geometric point coordinates are extracted, and tested before the node data is added leading to three cases:

          Case 1: If the geometric point contains a NaN, then the point is NOT added to the central node data store, last_point_idx=-1, and
          last_point_new=false.  

          Case 2: If the geometric point is already on the list and chk_point_unique is true, then the point is NOT added to the central node data store,
          last_point_idx is stet to the existing point's index, and last_point_new=false.  Only the geometric point components are used to determine if a
          point exists in the store.

          Case 3: If the geometric point is not on the list, then the point is added to the central node data store, last_point_idx is set to the new point's
          index, and last_point_new=true Note that last_point_idx is always the return value.

          @param node_data The point to add! */
      inline node_idx_t add_node(node_data_t node_data) {
        if (node_data_to_pnt.empty())
          node_data_to_pnt = {0, 1, 2};
        fvec3_t new_pnt = fvec3_from_node_data(node_data_to_pnt, node_data);
        if (pnt_has_nan(new_pnt)) {
          last_point_idx = -1;
          last_point_new = false;
        } else {
          if constexpr (chk_point_unique) {
            if (pnt_to_node_idx_map.contains(new_pnt)) {
              /* Point is already in list */
              last_point_idx = pnt_to_node_idx_map[new_pnt];
              last_point_new = false;
            } else {
              /* Point is not already in list */
              last_point_idx = node_count();
              pnt_to_node_idx_map[new_pnt] = last_point_idx;
              //node_idx_to_pnt.push_back(new_pnt);
              node_idx_to_node_data.push_back(node_data);
              last_point_new = true;
            }
          } else {
            /* Add point without regard to uniqueness */
            last_point_idx = node_count();
            //node_idx_to_pnt.push_back(new_pnt);
            node_idx_to_node_data.push_back(node_data);
            last_point_new = true;
          }
        }
        return last_point_idx;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Construct a fvec3_t from a node_data_t and a node_data_idx_lst_t.
          This is a speciazlied form of vector_from_node_data() that produces fvec3_t objects.  Primarly used to extract the location coordinates from node data.
          Note data_idx_lst must contain at least three values, but only the first three are used.  If the node_data_idx_lst_t has fewer than three elements, then
          this function will throw.
          @param data_idx_lst  Index list used to identify the data to pull
          @param node_data      node data vector from which to pull data */
      inline fvec3_t fvec3_from_node_data(node_data_idx_lst_t data_idx_lst, node_data_t node_data) const {
        fvec3_t v;
        for(int i=0; i<3; ++i)
          if (data_idx_lst[i].index() == 0)
            v[i] = node_data[std::get<int>(data_idx_lst.at(i))];
          else
            v[i] = std::get<uft_t>(data_idx_lst.at(i));
        return v;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Construct a vdat_t from a node_data_t and a node_data_idx_lst_t.
          @param data_idx_lst  Index list used to identify the data to pull
          @param node_data      node data vector from which to pull data */
      inline vdat_t vector_from_node_data(node_data_idx_lst_t data_idx_lst, node_data_t node_data) const {
        vdat_t v;
        for(decltype(data_idx_lst.size()) i=0; i<data_idx_lst.size(); ++i)
          if (data_idx_lst[i].index() == 0)
            v.push_back(node_data[std::get<int>(data_idx_lst.at(i))]);
          else
            v.push_back(std::get<uft_t>(data_idx_lst.at(i)));
        return v;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Extract a uft_t from a node_data_t and a node_data_idx_t.
          @param data_idx  Index list used to identify the data to pull
          @param node_data  node data vector from which to pull data */
      inline uft_t scalar_from_node_data(node_data_idx_t data_idx, node_data_t node_data) const {
        if (data_idx.index() == 0)
          return (node_data[std::get<int>(data_idx)]);
        else
          return (std::get<int>(data_idx));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert a fvec3_t to a string representation
          @param pnt_idx The point to convert.*/
      inline std::string node_to_string(node_idx_t pnt_idx) const {
        std::ostringstream convert;
        if (pnt_idx >= 0) {
          convert << "[ ";
          for(auto c: node_idx_to_node_data[pnt_idx]) 
            convert << std::setprecision(5) << static_cast<uft_t>(c) << " ";
          convert << "]";
        } else {
          convert << "[ DNE ]";
        }
        return(convert.str());
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Named Node Datasets. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the number of named scalar datasets */
      inline int named_scalar_datasets_count() const {
        return static_cast<int>(std::count_if(data_name_to_data_idx_lst.cbegin(), data_name_to_data_idx_lst.cend(), [](auto x) { return x.second.size()==1; }));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the number of named vector datasets */
      inline int named_vector_datasets_count() const {
        return static_cast<int>(std::count_if(data_name_to_data_idx_lst.cbegin(), data_name_to_data_idx_lst.cend(), [](auto x) { return x.second.size()!=1; }));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return the number of named (both vector and scalar) datasets */
      inline int named_datasets_count() const {
        return static_cast<int>(data_name_to_data_idx_lst.size());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Create named data sets by setting raw member.
          @param names value to set. */
      inline void set_data_name_to_data_idx_lst(data_name_to_node_data_idx_lst_t names) {
        data_name_to_data_idx_lst = names;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Create named data sets.
          Simple way to create named scalar datasets for the first n elements in each point's data.
          @param scalar_name_strings Scalar data set names to apply, in order, to node data elements. */
      inline void create_named_datasets(named_data_name_list_t scalar_name_strings) {
        data_name_to_data_idx_lst.clear();
        for(int i=0; i<static_cast<int>(scalar_name_strings.size()); ++i)
          data_name_to_data_idx_lst[scalar_name_strings[i]] = {i};
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** @overload
          @param scalar_name_strings Scalar data set names to apply, in order, to node data elements. 
          @param names               data_name_to_data_idx entries to create -- frequently used to add a named vector dataset. */
      inline void create_named_datasets(named_data_name_list_t scalar_name_strings, data_name_to_node_data_idx_lst_t names) {
        data_name_to_data_idx_lst.clear();
        for(int i=0; i<static_cast<int>(scalar_name_strings.size()); ++i)
          data_name_to_data_idx_lst[scalar_name_strings[i]] = {i};
        for(auto kv : names)
          data_name_to_data_idx_lst[kv.first] = kv.second;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Create mapping beteen node data and geometric pionts.
          Geometric points are constructed from node data via a node_data_idx_lst_t -- just like named data sets.  Each time a geometric point is required,
          it is dynamically pulled from the point's data via this mapping.
          @param point_data_index A node_data_idx_lst_t defining how to construct a geometric point. */
      inline void create_dataset_to_point_mapping(node_data_idx_lst_t point_data_index) {
        node_data_to_pnt = point_data_index;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name 3D Vector Computations. 

          These functions operate on fvec3_t objects (arrays of 3 real numbers), but treat them as 3-vectors. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the magnitude */
      inline uft_t vec3_two_norm(const fvec3_t& pnt) const {
        return std::sqrt(vec3_self_dot_product(pnt));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the self dot product -- i.e. magnitude squared */
      inline uft_t vec3_self_dot_product(const fvec3_t& pnt) const {
        uft_t tmp = 0.0;
        for(int i=0; i<3; ++i)
          tmp += pnt[i]*pnt[i];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the cross prodcut */
      inline uft_t vec3_dot_product(const fvec3_t& pnt1, const fvec3_t& pnt2) const {
        uft_t tmp = 0.0;
        for(int i=0; i<3; ++i)
          tmp += pnt1[i]*pnt2[i];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the cross prodcut */
      inline fvec3_t vec3_cross_product(const fvec3_t& pnt1, const fvec3_t& pnt2) const {
        fvec3_t tmp;
        tmp[0] = pnt1[1]*pnt2[2]-pnt1[2]*pnt2[1];
        tmp[1] = pnt1[2]*pnt2[0]-pnt1[0]*pnt2[2];
        tmp[2] = pnt1[0]*pnt2[1]-pnt1[1]*pnt2[0];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the vector diffrence */
      inline fvec3_t vec3_diff(const fvec3_t& pnt1, const fvec3_t& pnt2) const {
        fvec3_t tmp;
        for(int i=0; i<3; ++i)
          tmp[i] = pnt1[i] - pnt2[i];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the vector diffrence */
      inline uft_t vec3_scalar_triple_product(const fvec3_t& pnt1, const fvec3_t& pnt2, const fvec3_t& pnt3) const {
        return vec3_dot_product(pnt1, vec3_cross_product(pnt2, pnt3));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Unitize the given point in place.  Return true if the result is valid, and false otherwise */
      inline bool vec3_unitize(fvec3_t& pnt) const {
        uft_t length = vec3_two_norm(pnt);
        if (std::abs(length) > eps) {
          for(int i=0; i<3; ++i)
            pnt[i] = pnt[i]/length;
          return true;
        } else {
          return false;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the linear combination */
      inline fvec3_t vec3_linear_combination(uft_t scl1, const fvec3_t& pnt1, uft_t scl2, const fvec3_t& pnt2) const {
        fvec3_t tmp;
        for(int i=0; i<3; ++i)
          tmp[i] = scl1 * pnt1[i] + scl2 * pnt2[i];
        return tmp;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Determinant of 3x3 matrix with given vectors as columns/rows. */
      inline uft_t vec3_det3(const fvec3_t& pnt1, const fvec3_t& pnt2, const fvec3_t& pnt3) const {
        //  MJR TODO NOTE <2024-07-22T15:48:31-0500> vec3_det3: UNTESTED! UNTESTED! UNTESTED! UNTESTED!
        static_assert(false, "vec3_det3: Under active development.  Untested!");
        return (pnt1[0] * pnt2[1] * pnt3[2] -
                pnt1[0] * pnt2[2] * pnt3[1] -
                pnt1[1] * pnt2[0] * pnt3[2] +
                pnt1[1] * pnt2[2] * pnt3[0] +
                pnt1[2] * pnt2[0] * pnt3[1] -
                pnt1[2] * pnt2[1] * pnt3[0]);
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name 3D Geometry Computations.

          Methods with geomi_ prefix take point index types while methods with a geomr_ prefix take double tuples.  */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Segment intersection types */
      enum class seg_isect_t { C0_EMPTY,        //!< Cardinality=0.        Intersection is the empty set.
                               C1_VERTEX1,      //!< Cardinality=1.        Intersection is a single, shared vertex.
                               C1_INTERIOR,     //!< Cardinality=1.        Intersection is a single point in interior of at least one segment.
                               CI_VERTEX2,      //!< Cardinality=infinite. Intersection is a segment -- equal to input segments.
                               CI_VERTEX1,      //!< Cardinality=infinite. Intersection is a segment -- exactially one vertex in the intersection.
                               CI_VERTEX0,      //!< Cardinality=infinite. Intersection is a segment -- no vertexes included
                               BAD_SEGMENT,     //!< At least one of the segments was degenerate
                             };
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert an seg_isect_t to a string */
      inline std::string seg_isect_to_string(seg_isect_t seg_isect) const {
        switch(seg_isect) {
          case seg_isect_t::C0_EMPTY:     return std::string("C0_EMPTY");    break;
          case seg_isect_t::C1_VERTEX1:   return std::string("C1_VERTEX1");  break;
          case seg_isect_t::C1_INTERIOR:  return std::string("C1_INTERIOR"); break;
          case seg_isect_t::CI_VERTEX2:   return std::string("CI_VERTEX2");  break;
          case seg_isect_t::CI_VERTEX1:   return std::string("CI_VERTEX1");  break;
          case seg_isect_t::CI_VERTEX0:   return std::string("CI_VERTEX0");  break;
          case seg_isect_t::BAD_SEGMENT:  return std::string("BAD_SEGMENT"); break;
        }
        return std::string(""); // Never get here, but some compilers can't figure that out. ;)
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Determine the nature of the intersection between two line segments */
      seg_isect_t geomi_seg_isect_type(node_idx_t ilin1pnt1, node_idx_t ilin1pnt2, node_idx_t ilin2pnt1, node_idx_t ilin2pnt2) const {
        //  MJR TODO NOTE <2024-08-02T09:41:48-0500> geomi_seg_isect_type: Repeated point look-up slows things down
        //  MJR TODO NOTE <2024-08-02T09:41:48-0500> geomi_seg_isect_type: Add unit tests for each code branch
        // Check for degenerate segments
        if (ilin1pnt1 == ilin1pnt2)
          return seg_isect_t::BAD_SEGMENT;
        if (ilin2pnt1 == ilin2pnt2)
          return seg_isect_t::BAD_SEGMENT;
        // Count unique points & break into cases
        std::set<node_idx_t> points_sorted;
        points_sorted.insert(ilin1pnt1);
        points_sorted.insert(ilin1pnt2);
        points_sorted.insert(ilin2pnt1);
        points_sorted.insert(ilin2pnt2);
        if (points_sorted.size() == 4) { // ...................................................... REMAINING CASES: C0_EMPTY, C1_INTERIOR, CI_VERTEX0
          if (geomi_pts_colinear(ilin1pnt1, ilin1pnt2, ilin2pnt1, ilin2pnt2)) { // ............... REMAINING CASES: C0_EMPTY,              CI_VERTEX0
            if ( (geomi_pnt_line_distance(ilin1pnt1, ilin1pnt2, ilin2pnt1, true) < eps) ||
                 (geomi_pnt_line_distance(ilin1pnt1, ilin1pnt2, ilin2pnt2, true) < eps) ||
                 (geomi_pnt_line_distance(ilin2pnt1, ilin2pnt2, ilin1pnt1, true) < eps) ||
                 (geomi_pnt_line_distance(ilin2pnt1, ilin2pnt2, ilin1pnt2, true) < eps) ) { // .. REMAINING CASES: CI_VERTEX0
              return seg_isect_t::CI_VERTEX0;
            } else { // ......................................................................... REMAINING CASES: C0_EMPTY
              return seg_isect_t::C0_EMPTY;
            }
          } else { // ........................................................................... REMAINING CASES: C0_EMPTY, C1_INTERIOR
            if (geomi_seg_isect1(ilin1pnt1, ilin1pnt2, ilin2pnt1, ilin2pnt2)) { // .............. REMAINING CASES: C1_INTERIOR
              return seg_isect_t::C1_INTERIOR;
            } else { // ......................................................................... REMAINING CASES: C0_EMPTY
              return seg_isect_t::C0_EMPTY;
            }
          }
          return seg_isect_t::C0_EMPTY;
        } else if (points_sorted.size() == 3) { // .............................................. REMAINING CASES: C1_VERTEX1, CI_VERTEX1
          node_idx_t ipnt1, ipnt2, ipntc;
          if (ilin1pnt1 == ilin2pnt1) {
            ipntc = ilin1pnt1; ipnt1 = ilin1pnt2; ipnt2 = ilin2pnt2;
          } else if (ilin1pnt1 == ilin2pnt2) {
            ipntc = ilin1pnt1; ipnt1 = ilin1pnt2; ipnt2 = ilin2pnt1;
          } else if (ilin1pnt2 == ilin2pnt1) {
            ipntc = ilin1pnt2; ipnt1 = ilin1pnt1; ipnt2 = ilin2pnt2;
          } else if (ilin1pnt2 == ilin2pnt2) {
            ipntc = ilin1pnt2; ipnt1 = ilin1pnt1; ipnt2 = ilin2pnt1;
          } else { // Never get here.  Silences compiler warnings.
            ipntc = 0;         ipnt1 = 0;         ipnt2 = 0;
          }
          if (geomi_pts_colinear(ipnt1, ipnt2, ipntc)) { // ..................................... REMAINING CASES: C1_VERTEX1, CI_VERTEX1
            if ( (geomi_pnt_line_distance(ipnt1, ipntc, ipnt2, true) < eps) ||
                 (geomi_pnt_line_distance(ipnt2, ipntc, ipnt1, true) < eps) ) { // .............. REMAINING CASES: CI_VERTEX1
              return seg_isect_t::CI_VERTEX1;
            } else { // ......................................................................... REMAINING CASES: CI_VERTEX1
              return seg_isect_t::C1_VERTEX1;
            }
          } else { // ........................................................................... REMAINING CASES: C1_VERTEX1
            return seg_isect_t::C1_VERTEX1;
          }
        } else if (points_sorted.size() == 2) { // .............................................. REMAINING CASES: CI_VERTEX2
          return seg_isect_t::CI_VERTEX2;
        } else { // (points_sorted.size() == 1) which is an error. // ........................... REMAINING CASES: CI_VERTEX2
          return seg_isect_t::CI_VERTEX2;
        }
        return seg_isect_t::C0_EMPTY;  // Should never get here...
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if two line segments intersect in a single point */
      inline bool geomi_seg_isect1(node_idx_t ilin1pnt1, node_idx_t ilin1pnt2, node_idx_t ilin2pnt1, node_idx_t ilin2pnt2) const {
        return geomr_seg_isect1(get_pnt(ilin1pnt1), get_pnt(ilin1pnt2), get_pnt(ilin2pnt1), get_pnt(ilin2pnt2));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if two line segments intersect in a single point */
      inline bool geomr_seg_isect1(const fvec3_t& lin1pnt1, const fvec3_t& lin1pnt2, const fvec3_t& lin2pnt1, const fvec3_t& lin2pnt2) const {
        uft_t denom =
          lin1pnt1[0] * lin2pnt1[1] - lin1pnt1[0] * lin2pnt2[1] - lin1pnt1[1] * lin2pnt1[0] + lin1pnt1[1] * lin2pnt2[0] -
          lin1pnt2[0] * lin2pnt1[1] + lin1pnt2[0] * lin2pnt2[1] + lin1pnt2[1] * lin2pnt1[0] - lin1pnt2[1] * lin2pnt2[0];
        if (std::abs(denom) < eps) // Lines are parallel
          return false;
        uft_t numera =
          lin1pnt1[0]*lin2pnt1[1] - lin1pnt1[0]*lin2pnt2[1] -
          lin1pnt1[1]*lin2pnt1[0] + lin1pnt1[1]*lin2pnt2[0] +
          lin2pnt1[0]*lin2pnt2[1] - lin2pnt1[1]*lin2pnt2[0];
        uft_t numerb =
          -(lin1pnt1[0]*lin1pnt2[1] - lin1pnt1[0]*lin2pnt1[1] -
            lin1pnt1[1]*lin1pnt2[0] + lin1pnt1[1]*lin2pnt1[0] +
            lin1pnt2[0]*lin2pnt1[1] - lin1pnt2[1]*lin2pnt1[0]);
        uft_t ua = numera/denom;
        uft_t ub = numerb/denom;
        if ( (ua < 0) || (ub < 0) || (ua > 1) || (ub > 1) )
          return false;
        uft_t eq3 = numera/denom * (lin1pnt2[2] - lin1pnt1[2]) - numerb/denom * (lin2pnt2[2] - lin2pnt1[2]) + lin2pnt1[2] + lin1pnt1[2];
        if (std::abs(eq3) < eps) // Equation in third coordinate is satisfied
          return true;
        else
          return false;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Distance between a point and a line.
          See: geomr_pnt_line_distance(). */
      inline uft_t geomi_pnt_line_distance(node_idx_t ilinpnt1, node_idx_t ilinpnt2, node_idx_t ipnt, bool seg_distance) const {
        return geomr_pnt_line_distance(get_pnt(ilinpnt1), get_pnt(ilinpnt2), get_pnt(ipnt), seg_distance);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Distance between a point and a line.
          The result depends on the value of seg_distance:
          - seg_distance==false: Distance between the point pnt and the line containing linpnt1 & linpnt2.
          - seg_distance==true:  Distance between the point pnt and the line segment with endpoints linpnt1 & linpnt2. */
      uft_t geomr_pnt_line_distance(const fvec3_t& linpnt1, const fvec3_t& linpnt2, const fvec3_t& pnt, bool seg_distance) const {
        uft_t segd = geomr_pnt_pnt_distance(linpnt1, linpnt2);
        fvec3_t d, p;
        uft_t t = 0;
        for(int i=0; i<3; i++) {
          d[i] = (linpnt2[i] - linpnt1[i]) / segd;
          t += (pnt[i] - linpnt1[i])*d[i];
        }
        for(int i=0; i<3; i++) {
          p[i] = linpnt1[i] + t * d[i];
        }
        // p is the point on the line nearest pnt
        if (seg_distance) {
          uft_t dp1 = geomr_pnt_pnt_distance(linpnt1, p);
          uft_t dp2 = geomr_pnt_pnt_distance(linpnt2, p);
          //  MJR TODO NOTE <2024-08-02T09:42:19-0500> geomr_pnt_line_distance: check logic -- dp1>segd || dp2>segd?
          if (std::abs((dp1+dp2)-segd) > eps)
            return std::min(geomr_pnt_pnt_distance(linpnt1, pnt), geomr_pnt_pnt_distance(linpnt2, pnt));
        }
        return geomr_pnt_pnt_distance(p, pnt);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the euculedian (2 norm) distance between two points */
      inline uft_t geomr_pnt_pnt_distance(const fvec3_t& pnt1, const fvec3_t& pnt2) const {
        return vec3_two_norm(vec3_diff(pnt1, pnt2));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the normal of a triangle (or a plane defined by 3 points) */
      inline fvec3_t geomr_tri_normal(const fvec3_t& tripnt1, const fvec3_t& tripnt2, const fvec3_t& tripnt3, bool unit) const {
        fvec3_t basisv1 = vec3_diff(tripnt1, tripnt2);  // basis vectors for pln containing triagnel
        fvec3_t basisv2 = vec3_diff(tripnt3, tripnt2);  // basis vectors for pln containing triagnel
        fvec3_t normal = vec3_cross_product(basisv1, basisv2); // normal vector for tri. n=pld1xpld2
        if (unit)
          vec3_unitize(normal);
        return normal;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the distance between a point and a plane */
      inline uft_t geomr_pnt_pln_distance(const fvec3_t& plnpnt1, const fvec3_t& plnpnt2, const fvec3_t& plnpnt3, const fvec3_t& pnt) const {
        fvec3_t n = geomr_tri_normal(plnpnt1, plnpnt2, plnpnt3, true);
        return std::abs(vec3_dot_product(n, pnt) - vec3_dot_product(n, plnpnt2));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute the distance between a point and a triangle */
      inline uft_t geomr_pnt_tri_distance(const fvec3_t& tripnt1, const fvec3_t& tripnt2, const fvec3_t& tripnt3, const fvec3_t& pnt) const {
        //  MJR TODO NOTE <2024-07-22T15:48:31-0500> geomr_pnt_tri_distance: UNTESTED! UNTESTED! UNTESTED! UNTESTED! UNTESTED!
        fvec3_t basisv1 = vec3_diff(tripnt1, tripnt2);  // basis vectors for pln containing triagnel
        fvec3_t basisv2 = vec3_diff(tripnt3, tripnt2);  // basis vectors for pln containing triagnel
        fvec3_t normal = vec3_cross_product(basisv1, basisv2); // normal vector for tri. ax+by+cz+d=0, a=normal[0], b=normal[1], c=normal[2]
        vec3_unitize(normal);
        uft_t d = -vec3_dot_product(normal, tripnt2);            // ax+by+cz+d=0
        uft_t lambda = vec3_dot_product(normal, pnt) + d;
        fvec3_t q = vec3_diff(vec3_linear_combination(1.0, pnt, lambda, normal), tripnt2); // q is the point in the plane closest to pnt
        uft_t denom =  basisv1[0] * basisv2[1] - basisv2[0] * basisv1[1]; // If zero, then triangle is broken!
        uft_t u     = (q[0]       * basisv2[1] - basisv2[0] *       q[1]) / denom;
        uft_t v     = (basisv1[0] *       q[1] -       q[0] * basisv1[1]) / denom;
        uft_t dd    = std::abs(u*basisv1[2] + v*basisv2[2] - q[2]);
        if ( (u>=0) && (v>=0) && ((u+v)<=1) && (dd<eps) ) { // q is in the triangle =>  Use the plane distance
          return std::abs(lambda);
        } else {                                            // q is not in the triangle =>  Distance will be minimum distance to an edge.
          uft_t d1 = geomr_pnt_line_distance(tripnt1, tripnt2, pnt, true);
          uft_t d2 = geomr_pnt_line_distance(tripnt2, tripnt3, pnt, true);
          uft_t d3 = geomr_pnt_line_distance(tripnt3, tripnt1, pnt, true);
          return std::min(std::min(d1, d2), d3);
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points (just one point for this function) are zero.
          Checks if the infinity norm is less than epsilon.*/
      inline bool geomr_pnt_zero(const fvec3_t& p1) const {
        return ((std::abs(p1[0]) < eps) &&
                (std::abs(p1[1]) < eps) &&
                (std::abs(p1[2]) < eps));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are colinear */
      inline bool geomi_pts_colinear(node_idx_t pi1, node_idx_t pi2, node_idx_t pi3, node_idx_t pi4) const {
        return ( geomr_pts_colinear(get_pnt(pi1), get_pnt(pi2), get_pnt(pi3)) &&
                 geomr_pts_colinear(get_pnt(pi1), get_pnt(pi2), get_pnt(pi4)) );
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are colinear. */
      inline bool geomi_pts_colinear(node_idx_t pi1, node_idx_t pi2, node_idx_t pi3) const {
        return geomr_pts_colinear(get_pnt(pi1), get_pnt(pi2), get_pnt(pi3));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are colinear */
      inline bool geomr_pts_colinear(fvec3_t p1, fvec3_t p2, fvec3_t p3) const {
        return geomr_pnt_zero(vec3_cross_product(vec3_diff(p1, p2), vec3_diff(p1, p3)));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are coplanar */
      inline bool geomi_pts_coplanar(const node_idx_list_t& pnt_list) const {
        if (pnt_list.size() > 3) {
          if ( !(geomi_pts_coplanar(pnt_list[0], pnt_list[1], pnt_list[2], pnt_list[3])))
            return false;
          for(decltype(pnt_list.size()) i=4; i<pnt_list.size(); i++)
            if ( !(geomi_pts_coplanar(pnt_list[0], pnt_list[1], pnt_list[2], pnt_list[i])))
              return false;
        }
        return true;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are coplanar */
      inline bool geomi_pts_coplanar(node_idx_t pi1, node_idx_t pi2, node_idx_t pi3, node_idx_t pi4) const {
        return geomr_pts_coplanar(get_pnt(pi1), get_pnt(pi2), get_pnt(pi3), get_pnt(pi4));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if points are coplanar */
      inline bool geomr_pts_coplanar(fvec3_t p1, fvec3_t p2, fvec3_t p3, fvec3_t p4) const {
        return (std::abs(vec3_scalar_triple_product(vec3_diff(p3, p1), vec3_diff(p2, p1), vec3_diff(p4, p1))) < eps);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Compute intersection of a segment and a triangle */
      inline bool geomr_seg_tri_intersection(fvec3_t tp1, fvec3_t tp2, fvec3_t tp3, fvec3_t sp1, fvec3_t sp2) const {
        //  MJR TODO NOTE <2024-07-11T16:08:16-0500> geomr_seg_tri_intersection: implement
        //  MJR TODO NOTE <2024-07-11T16:08:27-0500> geomr_seg_tri_intersection: Should this be a bool or an enum?
        static_assert(false, "geomr_seg_tri_intersection: Not yet implemented!");
        return (tp1[0]+tp2[0]+tp3[0]+sp1[0]+sp2[0]>1.0);
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Class utilities. */
      //@{
      void clear() {
        last_point_idx = -1;
        last_point_new = true;
        pnt_to_node_idx_map.clear();
        node_idx_to_node_data.clear();
        data_name_to_data_idx_lst.clear();
        node_data_to_pnt.clear();
        cell_lst.clear();
        uniq_cell_lst.clear();
        last_cell_new = true;
        last_cell_stat = cell_stat_t::GOOD;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cells. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Cell segment/face/etc structure.  */
      typedef cell_lst_t cell_structure_t;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Number of cells. */
      inline int num_cells() const {
        return static_cast<int>(cell_lst.size());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_kind_t value to a list of vertexes. */
      inline cell_structure_t& cell_kind_to_structure(cell_kind_t cell_kind, int dimension) const {
        //  MJR TODO NOTE <2024-07-11T16:06:42-0500> cell_kind_to_structure: make sure polygons are all oriented correctly
        int logical_dim = cell_kind_to_dimension(cell_kind);
        if ( (dimension < 0) || (dimension > logical_dim) )
          dimension = logical_dim;
        if (logical_dim > 3) {
          std::cerr << "MR_cell_cplx API USAGE ERROR: Maximum supported dimension is 3" << std::endl;
          exit(0);
        }
        int idx = 0;
        switch(cell_kind) {
          case cell_kind_t::POINT:       idx = 0; break;
          case cell_kind_t::SEGMENT:     idx = 1; break;
          case cell_kind_t::TRIANGLE:    idx = 2; break;
          case cell_kind_t::QUAD:        idx = 3; break;
          case cell_kind_t::PYRAMID:     idx = 4; break;
          case cell_kind_t::HEXAHEDRON:  idx = 5; break;
        }
        static std::vector<std::vector<cell_structure_t>> cst = {{{{0}},                             // vertex
                                                                  {{0},{1}},                         // segment
                                                                  {{0},{1},{2}},                     // triangle
                                                                  {{0},{1},{2},{3}},                 // Quad
                                                                  {{0},{1},{2},{3},                  // Pyramid: Base vertexes
                                                                   {4}},                             // Pyramid: Tip vertex
                                                                  {{0},{1},{2},{3},                  // Hexahedron: Back vertexes
                                                                   {4},{5},{6},{7}}},                // Hexahedron: Front vertexes
                                                                 {{},                                // vertex
                                                                  {{0,1}},                           // segment
                                                                  {{0,1},{1,2},{2,0}},               // triangle
                                                                  {{0,1},{1,2},{2,3},{3,0}},         // Quad
                                                                  {{0,1},{1,2},{2,3},{3,0},          // Pyramid: Base segments
                                                                   {0,4},{1,4},{2,4},{3,4}},         // Pyramid: Side sets
                                                                  {{0,1},{1,2},{2,3},{3,0},          // Hexahedron: Back segments
                                                                   {4,5},{5,6},{6,7},{7,4},          // Hexahedron: Front segments
                                                                   {0,4},{1,5},{2,6},{3,7}}},        // Hexahedron: Back to front segments
                                                                 {{},                                // vertex
                                                                  {},                                // segment
                                                                  {{0,1,2}},                         // triangle
                                                                  {{0,1,2,3}},                       // Quad
                                                                  {{0,1,2,3},                        // Pyramid: Base face
                                                                   {0,1,4},                          // Pyramid: Back face
                                                                   {1,2,4},                          // Pyramid: Left face
                                                                   {2,3,4},                          // Pyramid: Front face
                                                                   {3,0,4}},                         // Pyramid: Right faces
                                                                  {{0,1,2,3},                        // Hexahedron: Back face
                                                                   {4,5,6,7},                        // Hexahedron: Front face
                                                                   {0,3,7,4},                        // Hexahedron: Left face
                                                                   {2,3,7,6},                        // Hexahedron: Top face
                                                                   {1,2,6,5},                        // Hexahedron: Right face
                                                                   {0,1,4,5}}},                      // Hexahedron: Bottom face
                                                                 {{},                                // vertex
                                                                  {},                                // segment
                                                                  {},                                // triangle
                                                                  {},                                // Quad
                                                                  {{0,1,2,3,4}},                     // Pyramid
                                                                  {{0,1,2,3,4,5,6,7}}}};             // Hexahedron
        return (cst[dimension][idx]);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_kind_t value to the logical dimension of the cell.

          Note that the vertexes of a cell_kind_t::QUAD cell might not be coplainer; however, the cell is logically a 2D cell.  Similarly, the vertexes of a
          cell_kind_t::PYRAMID might be coplainer; however, the cell is logically a 3D cell.  */
      inline int cell_kind_to_dimension(cell_kind_t cell_kind) const {
        switch(cell_kind) {
          case cell_kind_t::POINT:       return (0); break;
          case cell_kind_t::SEGMENT:     return (1); break;
          case cell_kind_t::TRIANGLE:    return (2); break;
          case cell_kind_t::QUAD:        return (2); break;
          case cell_kind_t::PYRAMID:     return (3); break;
          case cell_kind_t::HEXAHEDRON:  return (3); break;
        }
        return -1;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_kind_t value to the number of points required for cell type. */
      inline int cell_kind_to_req_pt_cnt(cell_kind_t cell_kind) const {
        switch(cell_kind) {
          case cell_kind_t::POINT:       return (1); break;
          case cell_kind_t::SEGMENT:     return (2); break;
          case cell_kind_t::TRIANGLE:    return (3); break;
          case cell_kind_t::QUAD:        return (4); break;
          case cell_kind_t::PYRAMID:     return (5); break;
          case cell_kind_t::HEXAHEDRON:  return (8); break;
        }
        return -1;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_kind_t value to the VTK cell type integer. */
      inline int cell_kind_to_vtk_type(cell_kind_t cell_kind) const {
        switch(cell_kind) {
          case cell_kind_t::POINT:       return ( 1); break;
          case cell_kind_t::SEGMENT:     return ( 3); break;
          case cell_kind_t::TRIANGLE:    return ( 5); break;
          case cell_kind_t::QUAD:        return ( 9); break;
          case cell_kind_t::HEXAHEDRON:  return (12); break;
          case cell_kind_t::PYRAMID:     return (14); break;
        }
        return -1;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_kind_t value to the VTK cell type integer. */
      inline std::string cell_kind_to_string(cell_kind_t cell_kind) const {
        switch(cell_kind) {
          case cell_kind_t::POINT:       return ("POINT"     ); break;
          case cell_kind_t::SEGMENT:     return ("SEGMENT"   ); break;
          case cell_kind_t::TRIANGLE:    return ("TRIANGLE"  ); break;
          case cell_kind_t::QUAD:        return ("QUAD"      ); break;
          case cell_kind_t::HEXAHEDRON:  return ("HEXAHEDRON"); break;
          case cell_kind_t::PYRAMID:     return ("PYRAMID"   ); break;
        }
        return ""; // Never get here.
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert number of points in a cell to cell_kind_t. */
      inline cell_kind_t req_pt_cnt_to_cell_kind(std::vector<int>::size_type node_count) const {
        switch(node_count) {
          case 1: return (cell_kind_t::POINT);      break;
          case 2: return (cell_kind_t::SEGMENT);    break;
          case 3: return (cell_kind_t::TRIANGLE);   break;
          case 4: return (cell_kind_t::QUAD);       break;
          case 5: return (cell_kind_t::PYRAMID);    break;
          case 8: return (cell_kind_t::HEXAHEDRON); break;
        }
        return (cell_kind_t::POINT); // Never get here
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_stat_t to a bool (true if GOOD). */
      inline bool cell_stat_is_good(cell_stat_t cell_stat) const {
        return (cell_stat == cell_stat_t::GOOD);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_stat_t to a bool (true if BAD). */
      inline bool cell_stat_is_bad(cell_stat_t cell_stat) const {
        return (cell_stat != cell_stat_t::GOOD);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Convert cell_stat_t enum value to a string. */
      inline std::string cell_stat_to_string(cell_stat_t cell_stat) const {
        switch(cell_stat) {
          case cell_stat_t::GOOD:            return (std::string("GOOD"));            break;
          case cell_stat_t::TOO_FEW_PNT:     return (std::string("TOO_FEW_PNT"));     break;
          case cell_stat_t::TOO_MANY_PNT:    return (std::string("TOO_MANY_PNT"));    break;
          case cell_stat_t::NEG_PNT_IDX:     return (std::string("NEG_PNT_IDX"));     break;
          case cell_stat_t::BIG_PNT_IDX:     return (std::string("BIG_PNT_IDX"));     break;
          case cell_stat_t::DUP_PNT:         return (std::string("DUP_PNT"));         break;
          case cell_stat_t::DIM_LOW:         return (std::string("DIM_LOW"));         break;
          case cell_stat_t::BAD_EDGEI:       return (std::string("BAD_EDGEI"));       break;
          case cell_stat_t::BAD_FACEI:       return (std::string("BAD_FACEI"));       break;
          case cell_stat_t::FACE_BENT:       return (std::string("FACE_BENT"));       break;
          case cell_stat_t::CONCAVE:         return (std::string("CONCAVE"));         break;
        }
        return std::string("ERROR");
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Perform cell vertex checks

          @param cell_kind  The type of the potential cell.
          @param cell_verts The vertexes of the potential cell.  */
      inline cell_stat_t check_cell_vertexes(cell_kind_t cell_kind, cell_verts_t cell_verts) const {
        // Check number of points
        std::vector<int>::size_type num_verts   = cell_verts.size();
        std::vector<int>::size_type req_num_pts = cell_kind_to_req_pt_cnt(cell_kind);
        if (num_verts < req_num_pts)
          return cell_stat_t::TOO_FEW_PNT;
        if (num_verts > req_num_pts)
          return cell_stat_t::TOO_MANY_PNT;
        // Check for negative point index
        if (std::any_of(cell_verts.cbegin(), cell_verts.cend(), [](node_idx_t i) { return (i < 0); }))
          return cell_stat_t::NEG_PNT_IDX;
        // Check for too big point index
        if (std::any_of(cell_verts.cbegin(), cell_verts.cend(), [this](node_idx_t i) { return (i >= node_count()); }))
          return cell_stat_t::BIG_PNT_IDX;
        // Check for duplicate points
        if (num_verts > 1) {
          std::set<node_idx_t> cell_verts_sorted;
          for(node_idx_t pnt_idx: cell_verts) {
            if (cell_verts_sorted.contains(pnt_idx))
              return cell_stat_t::DUP_PNT;
            cell_verts_sorted.insert(pnt_idx);
          }
        }
        // Return GOOD
        return cell_stat_t::GOOD;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Perform cell dimension (2D cells must not be colinear & 3D cells must not be coplainer)

          @warning This function assumes that check_cell_vertexes() would have returned GOOD for the cell being checked!

          @warning Will not detect degenerate cell_kind_t::SEGMENT as that is handled by check_cell_vertexes.

          @param cell_kind  The type of the potential cell.
          @param cell_verts The vertexes of the potential cell.  */
      inline cell_stat_t check_cell_dimension(cell_kind_t cell_kind, cell_verts_t cell_verts) const {
        if (cell_kind == cell_kind_t::TRIANGLE) {
          if (geomi_pts_colinear(cell_verts[0], cell_verts[1], cell_verts[2]))
            return cell_stat_t::DIM_LOW;
        } else if (cell_kind == cell_kind_t::QUAD) {
          if (geomi_pts_colinear(cell_verts[0], cell_verts[1], cell_verts[2], cell_verts[3]))
            return cell_stat_t::DIM_LOW;
        } else if (cell_kind == cell_kind_t::HEXAHEDRON) {
          if ( geomi_pts_coplanar(cell_verts))
            return cell_stat_t::DIM_LOW;
        } else if (cell_kind == cell_kind_t::PYRAMID) {
          if ( geomi_pts_coplanar(cell_verts))
            return cell_stat_t::DIM_LOW;
        }
        return cell_stat_t::GOOD;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Checks that cell edges have expected intersections.

          Checks that every pair of cell edges has the correct intersection type (i.e. intersect in a single vertex or an empty intersection).  If a bad
          intersection is detected, then a cell_stat_t of cell_stat_t::BAD_EDGEI will be returned.  Otherwise cell_stat_t::GOOD will be returned.

          @warning This function assumes that check_cell_vertexes() would have returned GOOD for the cell being checked!

          Note that this function will detect some conditions caught by other checks.  For example, this function will return cell_stat_t::BAD_EDGEI for
          cell_kind_t::SEG cells for which check_cell_dimension() returns cell_stat_t::DUP_PNT.  Note the cell_stat_t values are different depending upon which
          check function is called.

          @param cell_kind  The type of the potential cell.
          @param cell_verts The vertexes of the potential cell.  */
      inline cell_stat_t check_cell_edge_intersections(cell_kind_t cell_kind, cell_verts_t cell_verts) const {
        cell_structure_t& segs = cell_kind_to_structure(cell_kind, 1);
        if ( !(segs.empty())) {
          for(decltype(segs.size()) i=0; i<segs.size()-1; i++) {
            for(decltype(segs.size()) j=i+1; j<segs.size(); j++) {
              std::set<node_idx_t> points_sorted;
              points_sorted.insert(segs[i][0]);
              points_sorted.insert(segs[i][1]);
              points_sorted.insert(segs[j][0]);
              points_sorted.insert(segs[j][1]);
              auto it = geomi_seg_isect_type(cell_verts[segs[i][0]], cell_verts[segs[i][1]], cell_verts[segs[j][0]], cell_verts[segs[j][1]]);
              if(points_sorted.size() == 4) {
                if (it != seg_isect_t::C0_EMPTY)
                  return cell_stat_t::BAD_EDGEI;
              } else if(points_sorted.size() == 3) {
                if (it != seg_isect_t::C1_VERTEX1)
                  return cell_stat_t::BAD_EDGEI;
              } else {
                return cell_stat_t::BAD_EDGEI;
              }
            }
          }
        }
        return cell_stat_t::GOOD;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Checks that cell faces have expected intersections.

          @param cell_kind  The type of the potential cell.
          @param cell_verts The vertexes of the potential cell.  */
      inline cell_stat_t check_cell_face_intersections(cell_kind_t cell_kind, cell_verts_t cell_verts) const {
        //  MJR TODO NOTE <2024-08-02T09:42:38-0500> check_cell_face_intersections: Implement
        if (cell_kind == cell_kind_t::HEXAHEDRON) {
          if ( geomi_pts_coplanar(cell_verts))
            return cell_stat_t::DIM_LOW;
        } else if (cell_kind == cell_kind_t::PYRAMID) {
          if ( geomi_pts_coplanar(cell_verts))
            return cell_stat_t::DIM_LOW;
        }
        return cell_stat_t::GOOD;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Check if the vertexes of a cell are coplainer.

          Returns cell_stat_t::FACE_BENT if the cell is *not* plainer, and cell_stat_t::GOOD if it is.

          Note that most visualization packages are tolerant of non-plainer faces, and can render them just fine.  For example, most applications will
          automatically split a non-plainer cell_kind_t::QUAD into two triangles.  That said such faces can be an issue when using a tessellation for
          computation.  Such issues are most severe for things like FEM, but they can also show up for common visualization computations like the extraction
          of level curves/surfaces.

          @param cell_kind  The type of the potential cell.
          @param cell_verts The vertexes of the potential cell.  */
      inline cell_stat_t check_cell_faces_plainer(cell_kind_t cell_kind, cell_verts_t cell_verts) const {
        const cell_structure_t& face_structures = cell_kind_to_structure(cell_kind, 2);
        for(auto face_structure: face_structures) {
          cell_verts_t face;
          for(auto idx: face_structure)
            face.push_back(cell_verts[idx]);
          if ( !(geomi_pts_coplanar(face)))
            return cell_stat_t::FACE_BENT;
        }
        return cell_stat_t::GOOD;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add parts of a cell of the specified dimension.

          Example 1: a valid cell_kind_t::PYRAMID cell added wth dimension == 2 results in three cell_kind_t::TRIANGLE cells and one cell_kind_t::QUAD cell.

          Example 2: a valid cell_kind_t::PYRAMID cell added wth dimension == 3 results in one cell_kind_t::PYRAMID cell added

          Example 3: a valid cell_kind_t::PYRAMID cell added wth dimension == 1 results in 8 cell_kind_t::SEG cells added

          @param new_cell_kind   The type of the new cell
          @param new_cell_verts  The vertexes of the new cell
          @param dimension The dimension of the parts to add.
          @return Number of cells added */
      inline int add_cell(cell_kind_t new_cell_kind, cell_verts_t new_cell_verts, int dimension) {
        int num_added = 0;
        if ( (dimension < 0) || (dimension >= cell_kind_to_dimension(new_cell_kind)) ) {
          if (add_cell(new_cell_kind, new_cell_verts))
            num_added++;
        } else { // We need to break the cell up into lower dimensional bits, and add the bits.
          const cell_structure_t& cell_parts = cell_kind_to_structure(new_cell_kind, dimension);
          for(auto cell_part: cell_parts) {
            cell_verts_t newer_cell_verts;
            for(auto idx: cell_part)
              newer_cell_verts.push_back(new_cell_verts[idx]);
            if (add_cell(req_pt_cnt_to_cell_kind(newer_cell_verts.size()), newer_cell_verts))
              num_added++;
          }
        }
        return num_added;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add a cell
          @param new_cell_kind   The type of the new cell
          @param new_cell_verts  The vertexes of the new cell
          @return A boolean indicateing success
          @retval true  The cell was added or had been added previously
          @retval false The cell could not be added (because of a failed geometric check) */
      inline bool add_cell(cell_kind_t new_cell_kind, cell_verts_t new_cell_verts) {
        // Check vertexes if required
        if constexpr (chk_cell_vertexes) {
          last_cell_stat = check_cell_vertexes(new_cell_kind, new_cell_verts);
          if (cell_stat_is_bad(last_cell_stat))
            return false;
        }
        // Check dimension if required
        if constexpr (chk_cell_dimension) {
          last_cell_stat = check_cell_dimension(new_cell_kind, new_cell_verts);
          if (cell_stat_is_bad(last_cell_stat))
            return false;
        }
        // Check edges
        if constexpr (chk_cell_edges) {
          last_cell_stat = check_cell_edge_intersections(new_cell_kind, new_cell_verts);
          if (cell_stat_is_bad(last_cell_stat))
            return false;
        }
        // Geom was good or we didn't need to check.
        if constexpr (chk_cell_unique) {
          cell_verts_t new_cell_sorted_verts = new_cell_verts;
          std::sort(new_cell_sorted_verts.begin(), new_cell_sorted_verts.end());
          if (uniq_cell_lst.contains(new_cell_sorted_verts)) {
            last_cell_new = false;
          } else {
            last_cell_new = true;
            cell_lst.push_back(new_cell_verts);
            uniq_cell_lst.insert(new_cell_sorted_verts);
          }
        } else {
          cell_lst.push_back(new_cell_verts);
        }
        return true;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Retruns the status of the last cell given to the add_cell() method.
          If (chk_cell_vertexes || chk_cell_dimension | chk_cell_edges) is true, this value is updated each time add_cell() is called.
          Otherwise its value is always cell_stat_t::GOOD. */
      inline cell_stat_t status_of_last_cell_added() const {
        return last_cell_stat;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Retruns true if the last cell given to the add_cell() was new and the return value of status_of_last_cell_added() is cell_stat_t::GOOD.
          If chk_cell_unique is true, this value is updated each time add_cell() is called.  Otherwise its value is always true. */
      inline bool last_cell_added_was_new() const {
        return last_cell_new;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name I/O. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Type returned by I/O functions */
      typedef int io_result;
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump to an XML VTK unstructured grid file

          Note: A point vector data set named "NORMALS" will be used for normal vectors.

          @param file_name   The name of the output file
          @param description This is included as a file comment at the start of the file.
          @return 0 if everything worked, and non-zero otherwise */
      io_result write_xml_vtk(std::string file_name, std::string description) {
        /* Check that we have data */
        if (node_count() <= 0) {
          std::cout << "ERROR(write_xml_vtk): No points!" << std::endl;
          return 1;
        }
        if (num_cells() <= 0) {
          std::cout << "ERROR(write_xml_vtk): No cells!" << std::endl;
          return 2;
        }
        /* Looks like we have data.  Let's open our file */
        std::ofstream out_stream;
        out_stream.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);
        if (out_stream.is_open()) {
          out_stream.imbue(std::locale::classic());
        } else {
          std::cout << "ERROR(write_xml_vtk): Could not open file!" << std::endl;
          return 3;
        }
        out_stream << "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>" << std::endl;
        out_stream << "<!-- " << description << " -->" << std::endl;
        out_stream << "  <UnstructuredGrid>" << std::endl;
        out_stream << "    <Piece NumberOfPoints='" << node_count() << "' NumberOfCells='" << num_cells() << "'>" << std::endl;
        if ( !(data_name_to_data_idx_lst.empty())) {
          std::string scalars_attr_value, vectors_attr_value, normals_attr_value;
          for (auto& kv : data_name_to_data_idx_lst)
            if (kv.second.size()==1)
              scalars_attr_value += (scalars_attr_value.empty() ? "" : " ") + kv.first;
            else
              if (kv.first == "NORMALS")
                normals_attr_value = "NORMALS";
              else
                vectors_attr_value += (vectors_attr_value.empty() ? "" : " ") + kv.first;
          out_stream << "      <PointData";
          if ( !(scalars_attr_value.empty()))
            out_stream << " Scalars='" << scalars_attr_value << "'";
          if ( !(normals_attr_value.empty()))
            out_stream << " Normals='" << normals_attr_value << "'";
          if ( !(vectors_attr_value.empty()))
            out_stream << " Vectors='" << vectors_attr_value << "'";
          out_stream << ">" << std::endl;
          for (auto& kv : data_name_to_data_idx_lst) {
            out_stream << "        <DataArray Name='" << kv.first << "' type='Float64' format='ascii' NumberOfComponents='" << kv.second.size() << "'>" << std::endl;
            out_stream << "          ";
            for (const auto& dv : node_idx_to_node_data) {
              for (auto& idx : kv.second)
                out_stream << std::setprecision(10) << scalar_from_node_data(idx, dv) << " ";
            }
            out_stream << std::endl << "        </DataArray>" << std::endl;
          }
          out_stream << "      </PointData>" << std::endl;
        }
        out_stream << "      <Points>" << std::endl;
        out_stream << "        <DataArray Name='Points' type='Float64' format='ascii' NumberOfComponents='3'>" << std::endl;
        for(node_idx_t pnt_idx=0; pnt_idx<static_cast<node_idx_t>(node_idx_to_node_data.size()); pnt_idx++) {
          fvec3_t pnt = get_pnt(pnt_idx);
          out_stream << "          " << std::setprecision(10) << pnt[0] << " " << pnt[1] << " " << pnt[2] << std::endl;
        }
        out_stream << "        </DataArray>" << std::endl;
        out_stream << "      </Points>" << std::endl;
        out_stream << "      <Cells>" << std::endl;
        out_stream << "        <DataArray type='Int32' Name='connectivity' format='ascii'>" << std::endl;
        for(auto& poly: cell_lst) {
          out_stream << "          " ;
          for(auto& vert: poly)
            out_stream << vert << " ";
          out_stream << std::endl;
        }
        out_stream << "        </DataArray>" << std::endl;
        out_stream << "        <DataArray type='Int32' Name='offsets' format='ascii'>" << std::endl;
        out_stream << "          ";
        std::vector<int>::size_type j = 0;
        for(auto& poly: cell_lst) {
          j += poly.size();
          out_stream << j << " ";
        }
        out_stream << std::endl;
        out_stream << "        </DataArray>" << std::endl;
        out_stream << "        <DataArray type='Int8' Name='types' format='ascii'>" << std::endl;
        out_stream << "          ";
        for(auto& poly: cell_lst)
          out_stream << cell_kind_to_vtk_type(req_pt_cnt_to_cell_kind(poly.size())) << " ";
        out_stream << std::endl;
        out_stream << "        </DataArray>" << std::endl;
        out_stream << "      </Cells>" << std::endl;
        out_stream << "    </Piece>" << std::endl;
        out_stream << "  </UnstructuredGrid>" << std::endl;
        out_stream << "</VTKFile>" << std::endl;

        /* Final newline */
        out_stream << std::endl;
        out_stream.close();
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump to a legacy VTK UNSTRUCTUREDGRID file.

          Note: A point vector data set named "NORMALS" will be used for normal vectors.
          Note: A point vector data set named "COLORS" will be used as COLOR_SCALARS.  Components should be in [0.0, 1.0].

          @param file_name   The name of the output file
          @param description This is the file description.
          @return 0 if everything worked, and non-zero otherwise */
      io_result write_legacy_vtk(std::string file_name, std::string description) {
        /* Check that we have data */
        if (node_count() <= 0) {
          std::cout << "ERROR(write_legacy_vtk): No points!" << std::endl;
          return 1;
        }
        if (num_cells() <= 0) {
          std::cout << "ERROR(write_legacy_vtk): No cells!" << std::endl;
          return 2;
        }
        /* Looks like we have data.  Let's open our file */
        std::ofstream out_stream;
        out_stream.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);
        if (out_stream.is_open()) {
          out_stream.imbue(std::locale::classic());
        } else {
          std::cout << "ERROR(write_legacy_vtk): Could not open file!" << std::endl;
          return 3;
        }
        /* Dump the header */
        out_stream << "# vtk DataFile Version 3.0" << std::endl;
        out_stream << description << std::endl;
        out_stream << "ASCII" << std::endl;
        out_stream << "DATASET UNSTRUCTURED_GRID" << std::endl;
        /* Dump the points */
        out_stream << "POINTS " << node_count() << " double" << std::endl;
        //for (const auto& pnt : node_idx_to_pnt)
        for(node_idx_t pnt_idx=0; pnt_idx<static_cast<node_idx_t>(node_idx_to_node_data.size()); pnt_idx++) {
          fvec3_t pnt = get_pnt(pnt_idx);
          out_stream << std::setprecision(10) << pnt[0] << " " << pnt[1] << " " << pnt[2] << std::endl;
        }
        /* Dump the cell data */
        std::vector<int>::size_type total_cells_ints = 0;
        for(auto& cell: cell_lst)
          total_cells_ints += (1+cell.size());
        out_stream << "CELLS " << num_cells() << " " << total_cells_ints << std::endl;
        for(auto& poly: cell_lst) {
          out_stream << poly.size() << " ";
          for(auto& vert: poly) {
            out_stream << vert << " ";
          }
          out_stream << std::endl;
        }
        out_stream << "CELL_KINDS " << num_cells() << std::endl;
        for(auto& poly: cell_lst)
          out_stream << cell_kind_to_vtk_type(req_pt_cnt_to_cell_kind(poly.size())) << std::endl;
        /* Dump point scalar data */
        if (data_name_to_data_idx_lst.size() > 0) {
          out_stream << "POINT_DATA " << node_count() << std::endl;
          if (named_scalar_datasets_count() > 0) {
            for (auto& kv : data_name_to_data_idx_lst) {
              if (kv.second.size() == 1) {
                out_stream << "SCALARS " << kv.first << " double 1" << std::endl;
                out_stream << "LOOKUP_TABLE default" << std::endl;
                for (const auto& dv : node_idx_to_node_data) {
                  uft_t v = scalar_from_node_data(kv.second[0], dv);
                  out_stream << std::setprecision(10) << v << std::endl;
                }
              }
            }
          }
          if (named_vector_datasets_count() > 0) {
            for (auto& kv : data_name_to_data_idx_lst) {
              if (kv.second.size() == 3) { // Legacy format only supports three element vectors
                if ("NORMALS" == kv.first) {
                  out_stream << "NORMALS " << kv.first << " double" << std::endl;
                } else if ("COLORS" == kv.first) {
                  out_stream << "COLOR_SCALARS " << kv.first << " 3" << std::endl;
                } else {
                  out_stream << "VECTORS " << kv.first << " double" << std::endl; 
                }
                for (const auto& dv : node_idx_to_node_data) {
                  vdat_t v = vector_from_node_data(kv.second, dv);
                  out_stream << std::setprecision(10) << v[0] << " " << v[1] << " " << v[2] << std::endl;
                }
              }
            }
          }
        }
        /* Final newline */
        out_stream << std::endl;
        out_stream.close();
        return 0;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump object to STDOUT.

          @param max_num_print Maximum number of points/cells to print.  Use 0 to print all points/cells. */
      void dump_cplx(int max_num_print) const {
        std::cout << "Meta Data" << std::endl;
        std::cout << "  Points .................. " << node_count() << std::endl;
        std::cout << "  Data Scalars Per Point .. " << node_idx_to_node_data.size() << std::endl;
        std::cout << "  Named Data Sets ......... " << named_datasets_count() << std::endl;
        std::cout << "    Scalar Data Sets ...... " << named_scalar_datasets_count() << std::endl;
        std::cout << "    Vector Data Sets ...... " << named_vector_datasets_count() << std::endl;
        std::cout << "  Cells ................... " << num_cells() << std::endl;
        if (node_count() > 0) {
          int num_nodes_printed = 0;
          std::cout << "NODES BEGIN (" << node_count() << ")" << std::endl;
          for(node_idx_t pnt_idx = 0; pnt_idx<node_count(); ++pnt_idx) {
            std::cout << "  " << pnt_idx << ": " << node_to_string(pnt_idx) << std::endl;
            num_nodes_printed++;
            if ((max_num_print > 0) && (num_nodes_printed >= max_num_print)) {
              std::cout << "  Maximum number of nodes reached.  Halting tree dump." << std::endl;
              break;
            }
          }
          std::cout << "NODES END" << std::endl;
        }
        if (num_cells() > 0) {
          int num_cells_printed = 0;
          std::cout << "CELLS BEGIN (" << num_cells() << ")" << std::endl;
          for(int i=0; i<num_cells(); i++) {
            std::cout << "  ";
            for(auto& vert: cell_lst[i]) 
              std::cout << vert << " ";
            std::cout << "   " << cell_kind_to_string(req_pt_cnt_to_cell_kind(cell_lst[i].size())) << std::endl;
            num_cells_printed++;
            if ((max_num_print > 0) && (num_cells_printed >= max_num_print)) {
              std::cout << "  Maximum number of cells reached.  Halting tree dump." << std::endl;
              break;
            }
          }
          std::cout << "CELLS END" << std::endl;
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Dump to a PLY file

          Note: If one of the point vector data sets is named "NORMALS", then it will be used for normal vectors.
          Note: If one of the point vector data sets is named "COLORS", then it will be used for color data.  Components should be in [0.0, 1.0].

          @param file_name   The name of the output file
          @param description For legacy files, this is the file description.  For XML files this is included as a file comment at the start of the file.
          @return 0 if everything worked, and non-zero otherwise */
      io_result write_ply(std::string file_name, std::string description) {
        /* Check that we have data */
        if (node_count() <= 0) {
          std::cout << "ERROR(write_ply): No points!" << std::endl;
          return 1;
        }
        if (num_cells() <= 0) {
          std::cout << "ERROR(write_ply): No cells!" << std::endl;
          return 2;
        }
        for(auto& poly: cell_lst) {
          cell_kind_t cell_kind = req_pt_cnt_to_cell_kind(poly.size());
          if ( !((cell_kind == cell_kind_t::TRIANGLE) || (cell_kind == cell_kind_t::QUAD))) {
            std::cout << "ERROR(write_ply): Cells must all be 2D (triangles or quads)!" << std::endl;
            return 2;
          }
        }
        /* Looks like we have data.  Let's open our file */
        std::ofstream out_stream;
        out_stream.open(file_name, std::ios::out | std::ios::binary | std::ios::trunc);
        if (out_stream.is_open()) {
          out_stream.imbue(std::locale::classic());
        } else {
          std::cout << "ERROR(write_ply): Could not open file!" << std::endl;
          return 3;
        }
        bool have_colors_data = data_name_to_data_idx_lst.contains("COLORS");
        bool have_normal_data = data_name_to_data_idx_lst.contains("NORMALS");
        out_stream << "ply" << std::endl;
        out_stream << "format ascii 1.0" << std::endl;
        out_stream << "comment software: Mitch Richling's MR_rect_tree package" << std::endl;
        out_stream << "comment note: " << description << std::endl;
        out_stream << "element vertex " << node_count() << std::endl;
        out_stream << "property float x" << std::endl;
        out_stream << "property float y" << std::endl;
        out_stream << "property float z" << std::endl;
        if (have_colors_data) {
          out_stream << "property uchar red" << std::endl;
          out_stream << "property uchar green" << std::endl;
          out_stream << "property uchar blue" << std::endl;
        }
        if (have_normal_data) {
          out_stream << "property float nx" << std::endl;
          out_stream << "property float ny" << std::endl;
          out_stream << "property float nz" << std::endl;
        }
        out_stream << "element face " << num_cells() << std::endl; // May need to be adjusted if cells are not triangles..
        out_stream << "property list uchar int vertex_index" << std::endl;
        out_stream << "end_header" << std::endl;
        // Dump Vertex Data
        for (int i=0; i<node_count(); i++) {
          fvec3_t pnt = get_pnt(i);
          out_stream << std::setprecision(10) << pnt[0] << " " << pnt[1] << " " << pnt[2];
          if (have_colors_data) {
            fvec3_t clr = fvec3_from_node_data(data_name_to_data_idx_lst["COLORS"], node_idx_to_node_data[i]);
            out_stream << " " << static_cast<int>(255*clr[0]) << " " << static_cast<int>(255*clr[1]) << " " << static_cast<int>(255*clr[2]);
          }
          if (have_normal_data) {
            fvec3_t nml = fvec3_from_node_data(data_name_to_data_idx_lst["NORMALS"], node_idx_to_node_data[i]);
            vec3_unitize(nml);
            out_stream << " " << std::setprecision(10) << nml[0] << " " << nml[1] << " " << nml[2];
          }
          out_stream << std::endl;
        }
        // Dump Cells
        for(auto& poly: cell_lst) {
          out_stream << poly.size() << " ";
          for(auto& vert: poly) {
            out_stream << vert << " ";
          }
          out_stream << std::endl;
        }
        /* Final newline */
        out_stream << std::endl;
        out_stream.close();
        return 0;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Point Tests & Predicates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return true if any component of the point has a NaN.
          @param test_pnt The point to test */
      inline bool pnt_has_nan(const fvec3_t& test_pnt) {
        return (std::isnan(test_pnt[0]) || std::isnan(test_pnt[1]) || std::isnan(test_pnt[2]));
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return 0 if point is close to the level, -1 if it is below the level, and 1 if it is above the level.
          @param test_pnt      The point to test
          @param axis_index    Which axis to compare to the level
          @param level         Level to test aginst
          @param close_epsilon Epsilon used to check for "closeness". */
      inline int pnt_vs_level(const fvec3_t& test_pnt, int axis_index, uft_t level, uft_t close_epsilon=epsilon) {
        if (std::abs(test_pnt[axis_index]-level) < close_epsilon)
          return 0;
        else if(test_pnt[axis_index] < level)
          return -1;
        else
          return 1;
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Complex Computation. */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Function that takes and returns a node_data_t (Node Data Transform) */
      typedef std::function<node_data_t(const node_data_t&)> p2data_func_t;                // pd2pd
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Function that takes node data and returns a float (Node Data SDF) */
      typedef std::function<uft_t(const node_data_t&)>       p2real_func_t;                // pd2r
      /** Function that takes node data and returns a bool (Node Data Predicate) */
      typedef std::function<bool(const node_data_t&)>        p2bool_func_t;                // pd2b
      /** Function that takes a cell and returns a bool (Cell Predicate) */
      typedef std::function<bool(const cell_verts_t&)>             c2bool_func_t;                // c2b
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Delete cells matching a predicate function.

          @param func       Predicate function -- if true we get rid of the cell */
      int cull_cells(c2bool_func_t func) {
        int idx_last_good = -1;
        int start_size = num_cells();
        for(int i=0; i<start_size; i++)
          if ( !(func(cell_lst[i]))) {
            idx_last_good++;
            if (idx_last_good != i)
              cell_lst[idx_last_good] = cell_lst[i];
          }
        idx_last_good++;
        cell_lst.resize(idx_last_good);
        return (start_size-num_cells());
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Add new cells with points from existing cells with the given coordinates negated.

          When zero_epsilon is positive, some existing points may be adjusted such that components on the flipped axis near zero will become zero.  When
          chk_point_unique is true, both the original point coordinates and the new point coordinates will be mapped to the same point index -- that is to say if
          you add a new point to the complex with the old coordinates or the new coordinates, you will get the same index.  Note that a point collision may
          occur if the adjusted point is equal-ish to another, existing point -- a message is printed when this occurs.

          @param flip_list            A vector of booleans (0 or 1).  1 indicates the coordinate in a node_data_t vector should be negated.
          @param zero_epsilon         If non-negative, will collapse flipped coordinates near zero to be precisely zero.
          @param reverse_vertex_order Reverse the order of vertexes in each cell -- useful if vertex order impacts rendering. */
      void mirror(std::vector<int> flip_list, uft_t zero_epsilon=epsilon*1000, bool reverse_vertex_order=true) {
        int num_start_cells = num_cells();
        for(int cell_idx=0; cell_idx<num_start_cells; ++cell_idx) {
          cell_verts_t new_cell;
          for(auto pidx: cell_lst[cell_idx]) {
            node_data_t od = node_idx_to_node_data[pidx];
            for(int flip_list_idx=0; flip_list_idx<static_cast<int>(flip_list.size()); ++flip_list_idx)
              if ((flip_list[flip_list_idx]) && (std::abs(od[flip_list_idx]) < zero_epsilon))
                od[flip_list_idx] = 0;
            node_idx_to_node_data[pidx] = od;
            fvec3_t new_old_pnt = fvec3_from_node_data(node_data_to_pnt, od);
            if constexpr (chk_point_unique) {
              if (pnt_to_node_idx_map.contains(new_old_pnt))
                if (pnt_to_node_idx_map[new_old_pnt] != pidx)
                  std::cout << "ERROR(mirror): Collapse caused collision!" << std::endl;
              pnt_to_node_idx_map[new_old_pnt] = pidx;
            }
            node_data_t nd = od;
            for(int flip_list_idx=0; flip_list_idx<static_cast<int>(flip_list.size()); ++flip_list_idx)
              if (flip_list[flip_list_idx])
                nd[flip_list_idx] = -nd[flip_list_idx];
            node_idx_t p = add_node(nd);
            new_cell.push_back(p);
          }
          if (reverse_vertex_order)
            std::reverse(new_cell.begin(), new_cell.end());
          add_cell(req_pt_cnt_to_cell_kind(new_cell.size()), new_cell);
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Clear the solution cache stored in edge_solver_sdf(). */
      void clear_cache_edge_solver_sdf() {
        edge_solver_sdf(nullptr, 0, 0, nullptr, 0.0);
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Solve an SDF function for zero on a line between two node data sets, and add the solution to the central node data.

          The SDF must have different signs on the given points.  If this is not the case, then the endpoint with SDF closest to zero will be returned.

          Uses bisection.  Halting conditions:
          - The SDF value at the most recent guess is within solve_epsilon of zero
          - The absolute difference in SDF values at the end points is less than solve_epsilon

          This function stores a cache of previous solution results.  This cache may be cleared by clear_cache_edge_solver_sdf().

          @param pnt_idx1       First edge vertex
          @param pnt_idx2       Second edge vertex
          @param sdf_func       Data function (if nullptr, then clear solver cache and return immediately)
          @param solve_epsilon  Used to detect SDF value near zero
          @param dat_func       Produce the node data vector for the newly solved point. */
      node_idx_t edge_solver_sdf(p2data_func_t dat_func, node_idx_t pnt_idx1, node_idx_t pnt_idx2, p2real_func_t sdf_func, uft_t solve_epsilon=epsilon/10) {
        // Solver cache.
        static std::unordered_map<node_idx_t, std::unordered_map<node_idx_t, node_idx_t>> edge_solver_cache;
        if (sdf_func == nullptr) {
          edge_solver_cache.clear();
          return pnt_idx1;
        }
        node_idx_t cache_key1 = std::min(pnt_idx1, pnt_idx2);
        node_idx_t cache_key2 = std::max(pnt_idx1, pnt_idx2);
        // Check to see if we solved this one before
        if (edge_solver_cache.contains(cache_key1))
          if (edge_solver_cache[cache_key1].contains(cache_key2))
            return edge_solver_cache[cache_key1][cache_key2];
        // Gotta do the work and solve...
        // Init neg&pos with the points that have eng&pos sdf values
        node_data_t neg_node_data;
        uft_t      neg_pnt_sdfv;
        node_data_t pos_node_data = node_idx_to_node_data[pnt_idx1];
        uft_t      pos_pnt_sdfv = sdf_func(pos_node_data);
        if (pos_pnt_sdfv > 0) {
          neg_node_data = node_idx_to_node_data[pnt_idx2];
          neg_pnt_sdfv = sdf_func(neg_node_data);
        } else {
          neg_node_data = pos_node_data;
          neg_pnt_sdfv = pos_pnt_sdfv;
          pos_node_data = node_idx_to_node_data[pnt_idx2];
          pos_pnt_sdfv = sdf_func(pos_node_data);
        }
        // Init sol_node_data to end point with sdf value nearest zero
        node_data_t sol_node_data;
        uft_t      sol_pnt_sdfv;
        if (std::abs(pos_pnt_sdfv) < std::abs(neg_pnt_sdfv)) {
          sol_node_data = pos_node_data;
          sol_pnt_sdfv = pos_pnt_sdfv;
        } else {
          sol_node_data = neg_node_data;
          sol_pnt_sdfv = neg_pnt_sdfv;
        }
        if (neg_pnt_sdfv < 0) {  // Just to make sure pos&neg are pos&neg...
          while ((std::abs(sol_pnt_sdfv) > solve_epsilon) && (pos_pnt_sdfv-neg_pnt_sdfv) > solve_epsilon) {
            for(decltype(pos_node_data.size()) i=0; i<pos_node_data.size(); i++)
              sol_node_data[i] = (pos_node_data[i] + neg_node_data[i])/static_cast<uft_t>(2.0);
            sol_pnt_sdfv = sdf_func(sol_node_data);
            if (sol_pnt_sdfv > 0) {
              pos_node_data = sol_node_data;
              pos_pnt_sdfv = sol_pnt_sdfv;
            } else {
              neg_node_data = sol_node_data;
              neg_pnt_sdfv = sol_pnt_sdfv;
            }
          }
        }
        // Add out point, update solver cache, and return index
        node_idx_t sol_pnt_idx = add_node(dat_func(sol_node_data));
        edge_solver_cache[cache_key1][cache_key2] = sol_pnt_idx;
        return sol_pnt_idx;
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Fold triangles that cross over an SDF function.
          @param data_func      Data function
          @param sdf_func       SDF function
          @param solve_epsilon  Used to detect SDF value near zero */
      void triangle_folder(p2data_func_t data_func, p2real_func_t sdf_func, uft_t solve_epsilon=epsilon/10) {
        clear_cache_edge_solver_sdf();
        int num_start_cells = num_cells();
        for(int cell_idx=0; cell_idx<num_start_cells; ++cell_idx) {
          if (static_cast<int>(cell_lst[cell_idx].size()) == cell_kind_to_req_pt_cnt(cell_kind_t::TRIANGLE)) {
            auto& cur_cell = cell_lst[cell_idx];
            // Find and count zeros, positive, and negative vertexes
            int zero_cnt= 0, plus_cnt= 0, negv_cnt= 0;
            int zero_loc=-1, plus_loc=-1, negv_loc=-1;
            for(int i=0; i<3; i++) {
              uft_t sdf_val = sdf_func(node_idx_to_node_data[cur_cell[i]]);
              if (std::abs(sdf_val) <= solve_epsilon) {
                zero_cnt++;
                zero_loc = i;
              } else {
                if (sdf_val < static_cast<uft_t>(0.0)) {
                  plus_cnt++;
                  plus_loc = i;
                } else {
                  negv_cnt++;
                  negv_loc = i;
                }
              }
            }
            const std::vector<std::array<int, 3>> pmat { {0, 1, 2}, {1, 2, 0}, {2, 0, 1}};
            std::array<int, 3> p;
            if ((zero_cnt == 0) && (plus_cnt > 0) && (negv_cnt >0)) { // three triangles
              if (plus_cnt == 1)
                p = pmat[plus_loc];
              else
                p = pmat[negv_loc];
              auto orgv0 = cur_cell[p[0]];
              auto orgv1 = cur_cell[p[1]];
              auto orgv2 = cur_cell[p[2]];
              auto newv1 = edge_solver_sdf(data_func, orgv0, orgv1, sdf_func, solve_epsilon);
              auto newv2 = edge_solver_sdf(data_func, orgv0, orgv2, sdf_func, solve_epsilon);
              if ((newv1 != orgv0) && (newv1 != orgv1) && (newv2 != orgv0) && (newv2 != orgv2)) {
                cur_cell[p[1]] = newv1; // Modify current triangle in place
                cur_cell[p[2]] = newv2; // Modify current triangle in place
                add_cell(cell_kind_t::TRIANGLE, {newv1, orgv1, newv2});  // Add new triangle
                add_cell(cell_kind_t::TRIANGLE, {orgv1, orgv2, newv2});  // Add new triangle
              }
            } else if ((zero_cnt == 1) && (plus_cnt == 1) && (negv_cnt == 1)) { // two triangles
              p = pmat[zero_loc];
              auto orgv0 = cur_cell[p[0]];
              auto orgv1 = cur_cell[p[1]];
              auto orgv2 = cur_cell[p[2]];
              auto newv0 = edge_solver_sdf(data_func, orgv1, orgv2, sdf_func, solve_epsilon);
              if ((newv0 != orgv1) && (newv0 != orgv2)) {
                cur_cell[2] = newv0; // Modify current triangle in place
                add_cell(cell_kind_t::TRIANGLE, {orgv0, newv0, orgv2});  // Add new triangle
              }
            }
          }
        }
      }
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Fold segments that cross over an SDF function.
          @param data_func      Data function
          @param sdf_func       SDF function
          @param solve_epsilon  Used to detect SDF value near zero */
      void segment_folder(p2data_func_t data_func, p2real_func_t sdf_func, uft_t solve_epsilon=epsilon/10) {
        clear_cache_edge_solver_sdf();
        int num_start_cells = num_cells();
        for(int cell_idx=0; cell_idx<num_start_cells; ++cell_idx) {
          if (static_cast<int>(cell_lst[cell_idx].size()) == cell_kind_to_req_pt_cnt(cell_kind_t::SEGMENT)) {
            auto& cur_cell = cell_lst[cell_idx];
            int plus_cnt=0,  negv_cnt=0;
            for(int i=0; i<2; i++) {
              uft_t sdf_val = sdf_func(node_idx_to_node_data[cur_cell[i]]);
              if (sdf_val < static_cast<uft_t>(0.0))
                plus_cnt++;
              else if (sdf_val > static_cast<uft_t>(0.0))
                negv_cnt++;
            }
            if ((plus_cnt == 1) && (negv_cnt == 1)) {
              auto orgv0 = cur_cell[0];
              auto orgv1 = cur_cell[1];
              auto newv  = edge_solver_sdf(data_func, orgv0, orgv1, sdf_func, solve_epsilon);
              if ((newv != orgv0) && (newv != orgv1)) {
                cur_cell[1] = newv; // Modify current segment in place
                add_cell(cell_kind_t::SEGMENT, {newv, orgv1});  // Add new segment
              }
            }
          }
        }
      }
      //@}

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      /** @name Cell Predicates */
      //@{
      //--------------------------------------------------------------------------------------------------------------------------------------------------------
      /** Return true if ALL vertexes are above the level by more than epsilon */
      bool cell_above_level(const cell_verts_t cell_verts, int level_index, uft_t level, uft_t level_epsilon=epsilon) {
        return std::all_of(cell_verts.cbegin(), cell_verts.cend(), [this, level_index, level, level_epsilon](int v) { return (node_idx_to_node_data[v][level_index] > level+level_epsilon); });
      }
      /** Return true if ALL vertexes are below the level by more than epsilon */
      bool cell_below_level(const cell_verts_t cell_verts, int level_index, uft_t level, uft_t level_epsilon=epsilon) {
        return std::all_of(cell_verts.cbegin(), cell_verts.cend(), [this, level_index, level, level_epsilon](int v) { return (node_idx_to_node_data[v][level_index] < level+level_epsilon); });
      }
      /** Return true if cell is near SDF boundry */
      bool cell_near_sdf_boundry(const cell_verts_t cell_verts, p2real_func_t sdf_function, uft_t sdf_epsilon=epsilon) {
        int pos_cnt=0, neg_cnt=0;
        for(auto v: cell_verts) {
          uft_t sv = sdf_function(node_idx_to_node_data[v]);
          if (mjr::math::fc::near_zero(sv, sdf_epsilon))
            return true;
          if (sv < 0) {
            neg_cnt++;
          } else {
            pos_cnt++;
          }
          if ((neg_cnt > 0) && (pos_cnt > 0))
            return true;
        }
        return false;
      }
      //@}

  };

  typedef MR_cell_cplx<false, false, false, false, false, double, 1.0e-3> MRccF3;
  typedef MR_cell_cplx<true,  true,  true,  true,  true,  double, 1.0e-3> MRccT3;

  typedef MR_cell_cplx<false, false, false, false, false, double, 1.0e-5> MRccF5;
  typedef MR_cell_cplx<true,  true,  true,  true,  true,  double, 1.0e-5> MRccT5;

  typedef MR_cell_cplx<false, false, false, false, false, double, 1.0e-9> MRccF9;
  typedef MR_cell_cplx<true,  true,  true,  true,  true,  double, 1.0e-9> MRccT9;
}

#define MJR_INCLUDE_MR_cell_cplx
#endif
