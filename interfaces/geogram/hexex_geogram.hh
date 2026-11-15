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

/// @file
#pragma once

#include <HexEx/TrajectoryLayout.hh>

#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

#include <string>
#include <vector>

#ifndef HEXEX_GEOGRAM_EXPORT
#  if defined(_WIN32) && defined(HEXEX_GEOGRAM_EXPORT_SYMBOLS)
#    define HEXEX_GEOGRAM_EXPORT __declspec(dllexport)
#  elif defined(_WIN32)
#    define HEXEX_GEOGRAM_EXPORT __declspec(dllimport)
#  else
#    define HEXEX_GEOGRAM_EXPORT
#  endif
#endif

namespace HexEx
{

/**
 * @brief Conversion of the trajectories into geogram meshes.
 *
 * All information is stored in geogram *attributes*, so that the result can be
 * written with geogram's mesh IO (.geogram, .mesh, ...) without losing
 * anything. The attribute names are short on purpose: the element an attribute
 * belongs to (vertex, edge, facet, cell) is part of its identity. An attribute
 * that can be "not there" is an `int` holding -1; attributes that always have
 * a value are `unsigned int`.
 *
 * @par toRefinedMesh() - the refined mesh, i.e. the subdivision of the input
 *      tetrahedral mesh along the integer grid planes. It contains vertices,
 *      edges, triangular facets *and* tetrahedra:
 *
 *      | element | attribute   | type         | meaning                                              |
 *      |---------|-------------|--------------|------------------------------------------------------|
 *      | vertex  | `kind`      | unsigned int | 0 inside a tet, 1 on a tet face, 2 on a tet edge,    |
 *      |         |             |              | 3 on a tet vertex                                    |
 *      | vertex  | `tet_vertex`| int          | the input vertex, -1                                 |
 *      | vertex  | `tet_edge`  | int          | the input edge, -1                                   |
 *      | vertex  | `tet_face`  | int          | the input face, -1                                   |
 *      | vertex  | `tet_cell`  | int          | an incident input cell, -1                           |
 *      | vertex  | `hex_vertex`| int          | the hex mesh vertex, -1                              |
 *      | vertex  | `grid`      | unsigned int | 1 if the vertex is a corner of its piece's unit cube |
 *      | edge    | `hex_edge`  | int          | the hex mesh edge, -1                                |
 *      | edge    | `hex_edge_step` | int      | position on that hex edge's polyline, -1             |
 *      | edge    | `tet_edge`  | int          | the input edge the edge lies on, -1                  |
 *      | edge    | `piece`     | int          | an incident piece (see the cells' `piece`), -1       |
 *      | edge    | `grid`      | unsigned int | 1 if the edge lies on a grid line                    |
 *      | facet   | `hex_face`  | int          | the hex mesh face, -1                                |
 *      | facet   | `tet_face`  | int          | the input face the facet lies on, -1                 |
 *      | facet   | `piece`     | int          | the piece on one side (see the cells' `piece`), -1   |
 *      | facet   | `kind`      | unsigned int | 0 tet boundary, 1 hex face, 2 piece interior         |
 *      | cell    | `piece`     | unsigned int | the piece (tet x unit cube)                          |
 *      | cell    | `tet_cell`  | unsigned int | the input tet                                        |
 *      | cell    | `hex_cell`  | int          | the hex mesh cell, -1                                |
 *
 *      The facet store contains exactly the piece boundaries (the interfaces
 *      between pieces and the domain boundary); the facets that are interior to
 *      a piece are implied by its tetrahedra.
 *
 * @par toHexMesh() - the extracted hex mesh:
 *
 *      | element | attribute   | type         | meaning                                              |
 *      |---------|-------------|--------------|------------------------------------------------------|
 *      | vertex  | `hex_vertex`| unsigned int | the vertex's own index                               |
 *      | vertex  | `n_refined` | unsigned int | number of refined vertices with this hex vertex      |
 *      | edge    | `hex_edge`  | unsigned int | the edge's own index                                 |
 *      | edge    | `n_refined` | unsigned int | number of refined edges on it                        |
 *      | edge    | `border`    | unsigned int | 1 if it has no trajectory (hole / boundary)          |
 *      | face    | `hex_face`  | unsigned int | the face's own index                                 |
 *      | face    | `n_refined` | unsigned int | number of refined facets on it                       |
 *      | face    | `border`    | unsigned int | 1 if it has no trajectory (hole / boundary)          |
 *      | cell    | `hex_cell`  | unsigned int | the cell's own index                                 |
 *      | cell    | `n_refined` | unsigned int | number of refined tetrahedra inside the cell         |
 *      | cell    | `n_pieces`  | unsigned int | number of pieces inside the cell                     |
 *
 * @see HexTrajectoryLayout, HexExtractor::getHexMeshWithTrajectory
 */
class HEXEX_GEOGRAM_EXPORT GeogramBridge
{
public:
    /**
     * @brief Initializes geogram.
     *
     * Called automatically by the conversion functions; call it explicitly if
     * you want to use geogram before/after the conversion.
     */
    static void initialize();

    /**
     * @brief Fills @p out with the refined mesh, including all attributes
     *        listed above.
     *
     * @param withCells if false, the refined tetrahedra are omitted and only
     *        the vertices, edges and facets (the trajectories of the hex
     *        vertices, edges and faces) are written
     *
     * The cells are connected at the end (`out.cells.connect()`), because the
     * .geogram writer stores the cell adjacency and its loader does not rebuild
     * it. The facet adjacency is not computed: the facets are the piece
     * boundaries, which are non-manifold where pieces meet.
     */
    static void toRefinedMesh(const HexTrajectoryLayout& layout, GEO::Mesh& out,
                              bool withCells = true);

    /**
     * @brief Fills @p out with the extracted hex mesh and its attributes.
     *
     * The cells are reordered into geogram's hexahedron convention
     * (`v0 = (0,0,0), v1 = (1,0,0), v2 = (0,1,0), v3 = (1,1,0), v4 = (0,0,1),
     * v5 = (1,0,1), v6 = (0,1,1), v7 = (1,1,1)`): the extractor's own order
     * comes from the local UV cubes of the parametrization and would make
     * geogram treat face diagonals as edges.
     *
     * The cells are connected at the end (`out.cells.connect()`), because the
     * .geogram writer stores the cell adjacency and its loader does not rebuild
     * it. The facet adjacency is not computed: the facets are the piece
     * boundaries, which are non-manifold where pieces meet.
     *
     * @param layout the layout the hex mesh was extracted with
     * @param vertices the hex mesh's vertex positions
     * @param connectivity the hex mesh's cells, faces and edges
     */
    static void toHexMesh(const HexTrajectoryLayout& layout,
                          const std::vector<Vec3d>& vertices,
                          const TrajectoryMeshView& connectivity,
                          GEO::Mesh& out);

    /**
     * @brief Writes a mesh with geogram's mesh IO.
     *
     * The format is chosen from the file name's extension, e.g. `.geogram`,
     * `.mesh`. Attributes survive in the formats that can store them
     * (`.geogram` and `.mesh` can, `.obj` and `.stl` cannot).
     *
     * @return true on success
     */
    static bool save(const GEO::Mesh& mesh, const std::string& filename);

    /// Convenience: refine + save (see toRefinedMesh() and save()).
    static bool saveRefinedMesh(const HexTrajectoryLayout& layout, const std::string& filename,
                                bool withCells = true);

    /// Convenience: hex mesh + save (see toHexMesh() and save()).
    static bool saveHexMesh(const HexTrajectoryLayout& layout,
                            const std::vector<Vec3d>& vertices,
                            const TrajectoryMeshView& connectivity,
                            const std::string& filename);

    /**
     * @brief Human readable list of the attributes of a mesh, one per line.
     *
     * Useful to check what actually made it into a file.
     */
    static std::string describe(const GEO::Mesh& mesh);
};

} // namespace HexEx
