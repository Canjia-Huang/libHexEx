/*
 * Copyright 2016 Computer Graphics Group, RWTH Aachen University
 * Author: Max Lyon <lyon@cs.rwth-aachen.de>
 *
 * This file is part of HexEx.
 *
 * HexEx is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * HexEx is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HexEx.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <array>
#include <string>
#include <vector>

#include "Typedefs.hh"
#include "HexEx/Config/Export.hh"

namespace HexEx
{

/**
 * @brief The trajectories of the extracted hex mesh elements through the input
 *        tetrahedral mesh.
 *
 * HexEx extracts a hex mesh from an integer-grid map by tracing grid edges
 * (darts) through the tetrahedra. This structure records where the extracted
 * elements came from: it contains a *refined mesh* that is a subdivision of the
 * input tetrahedral mesh, obtained by cutting every tet with the integer grid
 * planes, and for every element of the refined mesh it stores which element of
 * the extracted hex mesh (if any) it belongs to.
 *
 * The refined mesh is tetrahedralized: every piece `tet ∩ unit cube` is split
 * into tetrahedra by coning from the piece's centroid to the triangles of its
 * triangulated boundary. This makes all elements of the refined mesh storeable
 * in standard formats (a geogram mesh stores vertices, edges, triangular facets
 * and tetrahedra). The boundary triangles of a piece are triangulated
 * canonically (by the vertex set), so the refined facets of neighbouring pieces
 * match and the refined mesh is conforming.
 *
 * The refined mesh's *facet* and *edge* stores contain exactly the piece
 * boundaries (the interfaces between pieces and the domain boundary), which is
 * where the hex face/edge trajectories live; the facets/edges that are interior
 * to a piece are implied by its tetrahedra.
 *
 * The correspondence is established without a global integer grid frame: a
 * global frame does not exist in general, because the frame of an integer-grid
 * map rotates around the singular edges of the parametrization. The layout
 * therefore groups the pieces of every unit cube (the *groups*) by propagating
 * the local integer coordinates through the extractor's transition functions,
 * which are exact per face, and identifies the hex elements through the
 * extractor's own intermediate hex mesh.
 *
 * Semantics of the correspondence:
 *
 *  - every refined vertex that is a corner of its piece's unit cube is matched
 *    to a vertex of the hex mesh (`vertex_hex_vertex`), otherwise -1.
 *  - every refined edge that lies on two integer grid planes (i.e. on a grid
 *    line) belongs to a hex mesh edge (`edge_hex_edge`) and knows its position
 *    on that edge's polyline (`edge_hex_edge_step`), otherwise -1.
 *  - every refined facet that lies on an integer grid plane belongs to a hex
 *    mesh face (`facet_hex_face`), otherwise -1.
 *  - every refined cell belongs to a piece (`cell_piece`) and to the hex mesh
 *    cell of that piece (`cell_hex_cell`), -1 if the extractor produced no hex
 *    cell there (e.g. inside the holes it reports).
 *
 * @see HexExtractor::getHexMeshWithTrajectory
 */
struct HEXEX_EXPORT HexTrajectoryLayout
{
    // ------------------------------------------------------------------
    // The refined mesh (subdivision of the input tetrahedral mesh)
    // ------------------------------------------------------------------

    /// 3d position of every refined vertex
    std::vector<Vec3d> vertex_position;
    /**
     * Provenance of every refined vertex w.r.t. the input tet mesh:
     * 0 = inside a tet, 1 = in the interior of a tet face, 2 = in the
     * interior of a tet edge, 3 = a tet vertex. (This is the number of input
     * mesh faces the vertex lies on, capped at 3.)
     */
    std::vector<int> vertex_kind;
    std::vector<int> vertex_tet_cell;   ///< an incident input cell, -1
    std::vector<int> vertex_tet_face;   ///< the input face, -1
    std::vector<int> vertex_tet_edge;   ///< the input edge, -1
    std::vector<int> vertex_tet_vertex; ///< the input vertex, -1
    std::vector<int> vertex_hex_vertex; ///< the hex mesh vertex, -1
    /// 1 if the vertex is a corner of its piece's unit cube
    std::vector<unsigned char> vertex_is_grid;

    /// edge endpoints (2 vertex indices per edge)
    std::vector<std::array<int, 2> > edge_vertices;
    std::vector<int> edge_hex_edge;      ///< the hex mesh edge, -1
    std::vector<int> edge_hex_edge_step; ///< position of the edge on the edge's polyline, -1
    std::vector<int> edge_tet_edge;      ///< the input tet edge the edge lies on, -1
    std::vector<int> edge_piece;         ///< a piece incident to the edge, -1
    std::vector<unsigned char> edge_on_grid; ///< 1 if the edge lies on a grid line

    /// triangular facets (3 vertex indices per facet)
    std::vector<std::array<int, 3> > facet_vertices;
    std::vector<int> facet_hex_face;   ///< the hex mesh face, -1
    std::vector<int> facet_tet_face;   ///< the input face the facet lies on, -1
    /**
     * 0 = on the boundary of the input tet the piece lives in,
     * 1 = on a hex face (integer grid plane), 2 = interior of a piece.
     */
    std::vector<int> facet_kind;
    std::vector<int> facet_piece_a;    ///< the piece on one side
    std::vector<int> facet_piece_b;    ///< the piece on the other side, -1 on the domain boundary

    /// tetrahedra (4 vertex indices per cell)
    std::vector<std::array<int, 4> > cell_vertices;
    std::vector<int> cell_piece;    ///< the piece (index into pieces)
    std::vector<int> cell_tet_cell; ///< the input tet cell
    std::vector<int> cell_hex_cell; ///< the hex mesh cell, -1

    /**
     * @brief One piece: the intersection of one input tet with one unit cube
     *        of the integer grid.
     */
    struct Piece
    {
        int tet_cell = -1;               ///< the input tet
        std::array<int, 3> cube;         ///< the cube's minimal corner, in the tet's parameter space
        double volume = 0;               ///< the volume in 3d space
        int hex_cell = -1;               ///< the corresponding hex mesh cell, -1
        int group = -1;                  ///< the group of pieces that form one unit cube
        int first_cell = -1;             ///< index of the first refined cell (tetrahedron)
        int n_cells = 0;                 ///< number of refined cells
        int first_vertex = -1;           ///< index of the first refined vertex
    };
    std::vector<Piece> pieces;

    /// a group of pieces: all pieces of one unit cube of the integer grid
    struct Group
    {
        int hex_cell = -1;               ///< the corresponding hex mesh cell, -1
        std::vector<int> pieces;         ///< its pieces
        /// hex mesh vertex of the eight cube corners (corner i has x = i&1, y = (i>>1)&1, z = (i>>2)&1), -1
        std::array<int, 8> corners;
    };
    std::vector<Group> groups;

    // ------------------------------------------------------------------
    // The correspondence, per element of the extracted hex mesh
    // ------------------------------------------------------------------

    /// number of refined vertices/edges/facets/cells that belong to each hex element
    std::vector<int> hex_vertex_n_refined;
    std::vector<int> hex_edge_n_refined;
    std::vector<int> hex_face_n_refined;
    std::vector<int> hex_cell_n_refined;
    std::vector<int> hex_cell_n_tets;

    /// the refined elements that belong to each hex element
    std::vector<std::vector<int> > hex_vertex_vertices;
    std::vector<std::vector<int> > hex_edge_edges;  ///< ordered by `edge_hex_edge_step`
    std::vector<std::vector<int> > hex_edge_polyline;
    std::vector<std::vector<int> > hex_face_facets;
    std::vector<std::vector<int> > hex_cell_pieces;
    std::vector<std::vector<int> > hex_cell_cells;

    // ------------------------------------------------------------------
    // Diagnostics
    // ------------------------------------------------------------------

    size_t n_input_tets = 0;
    /// input tets whose parametrization is degenerate (zero volume); they are skipped
    size_t n_degenerate_tets = 0;
    /// pieces (tet ∩ unit cube) with positive volume
    size_t n_pieces = 0;
    /// pieces whose volume is numerically zero and which were dropped
    size_t n_degenerate_pieces = 0;
    /// pieces whose two distinct boundary vertices were welded together
    size_t n_vertex_weld_conflicts = 0;
    /// pieces whose tetrahedralization does not reproduce the piece's volume
    size_t n_pieces_with_bad_tetrahedralization = 0;
    size_t n_refined_vertices = 0;
    size_t n_refined_edges = 0;
    size_t n_refined_facets = 0;
    size_t n_refined_cells = 0;
    /// groups of pieces that form one unit cube
    size_t n_groups = 0;
    /// groups whose eight corner hex vertices could not all be determined
    size_t n_groups_without_all_corners = 0;
    /// groups that do not correspond to a hex mesh cell (holes)
    size_t n_groups_without_hex_cell = 0;
    /// groups that map to a hex mesh cell another group already claimed (folds)
    size_t n_duplicate_hex_cells = 0;
    /// groups whose corners were assigned different hex vertices by different pieces
    size_t n_corner_conflicts = 0;
    /// refined grid vertices that could not be matched to an extracted vertex
    size_t n_unmatched_grid_vertices = 0;
    /// refined edges on a grid line that could not be identified with a hex edge
    size_t n_unresolved_grid_edges = 0;
    /// refined facets on a grid plane that could not be identified with a hex face
    size_t n_unresolved_grid_facets = 0;
    /// hex mesh elements that have no trajectory element at all
    size_t n_hex_vertices_without_trajectory = 0;
    size_t n_hex_edges_without_trajectory = 0;
    size_t n_hex_faces_without_trajectory = 0;
    size_t n_hex_cells_without_trajectory = 0;
    /// hex edges whose polyline branches (a folded parametrization visits the edge twice)
    size_t n_hex_edges_with_branched_polyline = 0;
    /// refined edges with only one incident facet (the domain boundary, or a crack)
    size_t n_open_edges = 0;
    /// open refined edges that are not on the domain boundary
    size_t n_crack_edges = 0;
    /// sum of the volumes of the input tets (non-degenerate ones) and of the pieces
    double volume_input_tets = 0;
    double volume_pieces = 0;

    /// remove all data
    void clear();
    /// human readable summary, one line per item
    std::string describe() const;
};

/**
 * @brief Plain connectivity arrays of an extracted hex mesh.
 *
 * Used to connect the trajectory layout with the hex mesh without making the
 * layout depend on a particular mesh type.
 */
struct HEXEX_EXPORT TrajectoryMeshView
{
    int n_vertices = 0;
    std::vector<std::vector<int> > cells; ///< per cell: its vertices
    std::vector<std::vector<int> > faces; ///< per face: its vertices
    std::vector<std::array<int, 2> > edges; ///< per edge: its two vertices
};

/**
 * @brief Extracts the vertex positions and the connectivity of a hex mesh.
 *
 * This is the input the geogram bridge needs; it keeps the bridge independent
 * of the concrete mesh type.
 */
template <typename HexMeshT>
void makeTrajectoryMesh(const HexMeshT& hexMesh, std::vector<Vec3d>& vertices,
                        TrajectoryMeshView& connectivity)
{
    vertices.resize((size_t)hexMesh.n_vertices());
    for (auto vh : hexMesh.vertices())
    {
        const auto p = hexMesh.vertex(vh);
        vertices[(size_t)vh.idx()] = Vec3d(p[0], p[1], p[2]);
    }

    connectivity.n_vertices = (int)hexMesh.n_vertices();

    connectivity.cells.resize((size_t)hexMesh.n_cells());
    for (auto ch : hexMesh.cells())
    {
        std::vector<int>& cell = connectivity.cells[(size_t)ch.idx()];
        const auto range = hexMesh.cell_vertices(ch);
        for (auto it = range.first; it != range.second; ++it)
            cell.push_back((int)it->idx());
    }

    connectivity.faces.resize((size_t)hexMesh.n_faces());
    for (auto fh : hexMesh.faces())
    {
        const std::vector<VertexHandle> faceVertices =
                hexMesh.get_halfface_vertices(hexMesh.halfface_handle(fh, 0));
        std::vector<int>& face = connectivity.faces[(size_t)fh.idx()];
        for (size_t i = 0; i < faceVertices.size(); ++i)
            face.push_back((int)faceVertices[i].idx());
    }

    connectivity.edges.resize((size_t)hexMesh.n_edges());
    for (auto eh : hexMesh.edges())
    {
        const HalfEdgeHandle heh = hexMesh.halfedge_handle(eh, 0);
        connectivity.edges[(size_t)eh.idx()][0] = (int)hexMesh.from_vertex_handle(heh).idx();
        connectivity.edges[(size_t)eh.idx()][1] = (int)hexMesh.to_vertex_handle(heh).idx();
    }
}

} // namespace HexEx
