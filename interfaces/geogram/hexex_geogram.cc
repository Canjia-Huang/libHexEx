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

#include "hexex_geogram.hh"

#include <geogram/basic/attributes.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <geogram/mesh/mesh_io.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <unordered_set>

namespace HexEx
{

namespace
{

const char* const ATTR_KIND = "kind";
const char* const ATTR_TET_VERTEX = "tet_vertex";
const char* const ATTR_TET_EDGE = "tet_edge";
const char* const ATTR_TET_FACE = "tet_face";
const char* const ATTR_TET_CELL = "tet_cell";
const char* const ATTR_HEX_VERTEX = "hex_vertex";
const char* const ATTR_HEX_EDGE = "hex_edge";
const char* const ATTR_HEX_EDGE_STEP = "hex_edge_step";
const char* const ATTR_HEX_FACE = "hex_face";
const char* const ATTR_HEX_CELL = "hex_cell";
const char* const ATTR_GRID = "grid";
const char* const ATTR_GRID_AXIS = "grid_axis";
const char* const ATTR_PIECE = "piece";
const char* const ATTR_BORDER = "border";
const char* const ATTR_N_REFINED = "n_refined";
const char* const ATTR_N_PIECES = "n_pieces";

inline unsigned int asUnsigned(int v)
{
    geo_debug_assert(v >= 0);
    return v < 0 ? 0u : (unsigned int)v;
}

inline uint64_t edgeKey(int a, int b)
{
    if (a > b) std::swap(a, b);
    return ((uint64_t)(uint32_t)a << 32) | (uint64_t)(uint32_t)b;
}

/**
 * @brief Reorders the eight vertices of a hexahedron into geogram's convention.
 *
 * geogram's hexahedron (see GEO::hex_descriptor) assumes the vertex order
 *
 *     v0 = (0,0,0)   v1 = (1,0,0)   v2 = (0,1,0)   v3 = (1,1,0)
 *     v4 = (0,0,1)   v5 = (1,0,1)   v6 = (0,1,1)   v7 = (1,1,1)
 *
 * i.e. v1, v2 and v4 are the three edge neighbours of v0, and v3, v5 and v6
 * are the neighbours of (v1,v2), (v1,v4) and (v2,v4). The extractor's own
 * order comes from the local UV cubes of the parametrization and does not
 * follow this convention, so the cells have to be reordered before they are
 * handed to geogram - otherwise geogram interprets face diagonals as edges and
 * the cells are twisted.
 *
 * The three neighbours of v0 are assigned to v1, v2 and v4 such that the frame
 * is right-handed; the cell's shape does not depend on that choice.
 */
bool reorderHexCell(const std::array<int, 8>& cell,
                    const std::unordered_set<uint64_t>& edges,
                    const std::vector<Vec3d>& positions,
                    std::array<int, 8>& order)
{
    const auto adjacent = [&](int a, int b) {
        return edges.find(edgeKey(a, b)) != edges.end();
    };

    const int v0 = *std::min_element(cell.begin(), cell.end());

    std::vector<int> neighbors;
    for (size_t i = 0; i < cell.size(); ++i)
        if (cell[i] != v0 && adjacent(v0, cell[i])) neighbors.push_back(cell[i]);
    if (neighbors.size() != 3) return false;
    std::sort(neighbors.begin(), neighbors.end());

    // a right-handed assignment of the three neighbours to (v1, v2, v4)
    static const int permutations[6][3] = {
        {0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}
    };
    int assignment[3] = { 0, 1, 2 };
    for (int p = 0; p < 6; ++p)
    {
        const Vec3d d1 = positions[neighbors[permutations[p][0]]] - positions[v0];
        const Vec3d d2 = positions[neighbors[permutations[p][1]]] - positions[v0];
        const Vec3d d4 = positions[neighbors[permutations[p][2]]] - positions[v0];
        if (((d1 % d2) | d4) > 0)
        {
            for (int i = 0; i < 3; ++i) assignment[i] = permutations[p][i];
            break;
        }
    }

    const int v1 = neighbors[assignment[0]];
    const int v2 = neighbors[assignment[1]];
    const int v4 = neighbors[assignment[2]];

    const auto commonNeighbor = [&](int a, int b) -> int {
        for (size_t i = 0; i < cell.size(); ++i)
        {
            const int v = cell[i];
            if (v == v0 || v == a || v == b) continue;
            if (adjacent(a, v) && adjacent(b, v)) return v;
        }
        return -1;
    };

    const int v3 = commonNeighbor(v1, v2); // neighbour of v1 and v2
    const int v5 = commonNeighbor(v1, v4); // neighbour of v1 and v4
    const int v6 = commonNeighbor(v2, v4); // neighbour of v2 and v4
    if (v3 < 0 || v5 < 0 || v6 < 0) return false;

    int v7 = -1;
    for (size_t i = 0; i < cell.size(); ++i)
    {
        const int v = cell[i];
        if (v == v0 || v == v1 || v == v2 || v == v3 || v == v4 || v == v5 || v == v6) continue;
        v7 = v;
        break;
    }
    if (v7 < 0) return false;
    if (!adjacent(v3, v7) || !adjacent(v5, v7) || !adjacent(v6, v7)) return false;

    order = std::array<int, 8>{{ v0, v1, v2, v3, v4, v5, v6, v7 }};
    return true;
}

} // namespace

void GeogramBridge::initialize()
{
    static bool initialized = false;
    if (!initialized)
    {
        GEO::initialize();
        /*
         * geogram's ".geogram" writer asks for "sys:compression_level" (and
         * other system settings) while saving. Those variables only exist once
         * the corresponding argument group has been imported, and geogram
         * asserts if one of them is missing - so import the standard group.
         */
        GEO::CmdLine::import_arg_group("standard");
        initialized = true;
    }
}

void GeogramBridge::toRefinedMesh(const HexTrajectoryLayout& layout, GEO::Mesh& out, bool withCells)
{
    initialize();
    out.clear(false);

    // ---- vertices ------------------------------------------------------
    const GEO::index_t nVertices = (GEO::index_t)layout.vertex_position.size();
    for (GEO::index_t v = 0; v < nVertices; ++v)
    {
        const Vec3d& p = layout.vertex_position[v];
        out.vertices.create_vertex(p.data());
    }

    GEO::Attribute<unsigned int> vertexKind(out.vertices.attributes(), ATTR_KIND);
    GEO::Attribute<int> vertexTetVertex(out.vertices.attributes(), ATTR_TET_VERTEX);
    GEO::Attribute<int> vertexTetEdge(out.vertices.attributes(), ATTR_TET_EDGE);
    GEO::Attribute<int> vertexTetFace(out.vertices.attributes(), ATTR_TET_FACE);
    GEO::Attribute<int> vertexTetCell(out.vertices.attributes(), ATTR_TET_CELL);
    GEO::Attribute<int> vertexHexVertex(out.vertices.attributes(), ATTR_HEX_VERTEX);
    GEO::Attribute<unsigned int> vertexGrid(out.vertices.attributes(), ATTR_GRID);

    for (GEO::index_t v = 0; v < nVertices; ++v)
    {
        vertexKind[v] = asUnsigned(layout.vertex_kind[v]);
        vertexTetVertex[v] = layout.vertex_tet_vertex[v];
        vertexTetEdge[v] = layout.vertex_tet_edge[v];
        vertexTetFace[v] = layout.vertex_tet_face[v];
        vertexTetCell[v] = layout.vertex_tet_cell[v];
        vertexHexVertex[v] = layout.vertex_hex_vertex[v];
        vertexGrid[v] = layout.vertex_is_grid[v] ? 1u : 0u;
    }

    // ---- edges ---------------------------------------------------------
    GEO::Attribute<int> edgeHexEdge(out.edges.attributes(), ATTR_HEX_EDGE);
    GEO::Attribute<int> edgeHexEdgeStep(out.edges.attributes(), ATTR_HEX_EDGE_STEP);
    GEO::Attribute<int> edgeTetEdge(out.edges.attributes(), ATTR_TET_EDGE);
    GEO::Attribute<int> edgePiece(out.edges.attributes(), ATTR_PIECE);
    GEO::Attribute<unsigned int> edgeGrid(out.edges.attributes(), ATTR_GRID);
    GEO::Attribute<int> edgeGridAxis(out.edges.attributes(), ATTR_GRID_AXIS);

    for (size_t e = 0; e < layout.edge_vertices.size(); ++e)
    {
        const GEO::index_t index = out.edges.create_edge(
                (GEO::index_t)layout.edge_vertices[e][0],
                (GEO::index_t)layout.edge_vertices[e][1]);
        edgeHexEdge[index] = layout.edge_hex_edge[e];
        edgeHexEdgeStep[index] = layout.edge_hex_edge_step[e];
        edgeTetEdge[index] = layout.edge_tet_edge[e];
        edgePiece[index] = layout.edge_piece[e];
        edgeGrid[index] = layout.edge_on_grid[e] ? 1u : 0u;
        edgeGridAxis[index] = layout.edge_grid_axis[e];
    }

    // ---- facets --------------------------------------------------------
    GEO::Attribute<int> facetHexFace(out.facets.attributes(), ATTR_HEX_FACE);
    GEO::Attribute<int> facetTetFace(out.facets.attributes(), ATTR_TET_FACE);
    GEO::Attribute<int> facetPiece(out.facets.attributes(), ATTR_PIECE);
    GEO::Attribute<unsigned int> facetKind(out.facets.attributes(), ATTR_KIND);

    for (size_t f = 0; f < layout.facet_vertices.size(); ++f)
    {
        const GEO::index_t index = out.facets.create_triangle(
                (GEO::index_t)layout.facet_vertices[f][0],
                (GEO::index_t)layout.facet_vertices[f][1],
                (GEO::index_t)layout.facet_vertices[f][2]);
        facetHexFace[index] = layout.facet_hex_face[f];
        facetTetFace[index] = layout.facet_tet_face[f];
        facetPiece[index] = layout.facet_piece_a[f];
        facetKind[index] = asUnsigned(layout.facet_kind[f]);
    }

    // ---- cells ---------------------------------------------------------
    const GEO::index_t nCells = withCells ? (GEO::index_t)layout.cell_vertices.size() : 0;
    if (nCells > 0)
    {
        const GEO::index_t first = out.cells.create_tets(nCells);
        for (GEO::index_t c = 0; c < nCells; ++c)
            for (GEO::index_t i = 0; i < 4; ++i)
                out.cells.set_vertex(first + c, i, (GEO::index_t)layout.cell_vertices[c][i]);

        /*
         * Build the cell adjacency: geogram does not derive it automatically,
         * and the .geogram writer stores it (the loader skips connect() for
         * that format). Sliver removal is disabled: no cell may be removed,
         * because the attributes are indexed by cell.
         */
        out.cells.connect(false);

        GEO::Attribute<unsigned int> cellPiece(out.cells.attributes(), ATTR_PIECE);
        GEO::Attribute<unsigned int> cellTetCell(out.cells.attributes(), ATTR_TET_CELL);
        GEO::Attribute<int> cellHexCell(out.cells.attributes(), ATTR_HEX_CELL);
        for (GEO::index_t c = 0; c < nCells; ++c)
        {
            cellPiece[first + c] = asUnsigned(layout.cell_piece[c]);
            cellTetCell[first + c] = asUnsigned(layout.cell_tet_cell[c]);
            cellHexCell[first + c] = layout.cell_hex_cell[c];
        }
    }
}

void GeogramBridge::toHexMesh(const HexTrajectoryLayout& layout,
                              const std::vector<Vec3d>& vertices,
                              const TrajectoryMeshView& connectivity,
                              GEO::Mesh& out)
{
    initialize();
    out.clear(false);

    for (size_t v = 0; v < vertices.size(); ++v)
        out.vertices.create_vertex(vertices[v].data());

    GEO::Attribute<unsigned int> vertexHexVertex(out.vertices.attributes(), ATTR_HEX_VERTEX);
    GEO::Attribute<unsigned int> vertexNRefined(out.vertices.attributes(), ATTR_N_REFINED);
    for (GEO::index_t v = 0; v < out.vertices.nb(); ++v)
    {
        vertexHexVertex[v] = (unsigned int)v;
        vertexNRefined[v] = (v < layout.hex_vertex_n_refined.size())
                ? asUnsigned(layout.hex_vertex_n_refined[v]) : 0u;
    }

    GEO::Attribute<unsigned int> edgeHexEdge(out.edges.attributes(), ATTR_HEX_EDGE);
    GEO::Attribute<unsigned int> edgeNRefined(out.edges.attributes(), ATTR_N_REFINED);
    GEO::Attribute<unsigned int> edgeBorder(out.edges.attributes(), ATTR_BORDER);
    for (size_t e = 0; e < connectivity.edges.size(); ++e)
    {
        const GEO::index_t index = out.edges.create_edge(
                (GEO::index_t)connectivity.edges[e][0],
                (GEO::index_t)connectivity.edges[e][1]);
        edgeHexEdge[index] = (unsigned int)e;
        const bool hasTrajectory = (e < layout.hex_edge_n_refined.size() && layout.hex_edge_n_refined[e] > 0);
        edgeNRefined[index] = asUnsigned(e < layout.hex_edge_n_refined.size() ? layout.hex_edge_n_refined[e] : 0);
        edgeBorder[index] = hasTrajectory ? 0u : 1u;
    }

    GEO::Attribute<unsigned int> facetHexFace(out.facets.attributes(), ATTR_HEX_FACE);
    GEO::Attribute<unsigned int> facetNRefined(out.facets.attributes(), ATTR_N_REFINED);
    GEO::Attribute<unsigned int> facetBorder(out.facets.attributes(), ATTR_BORDER);
    for (size_t f = 0; f < connectivity.faces.size(); ++f)
    {
        const std::vector<int>& face = connectivity.faces[f];
        if (face.size() < 3) continue;
        GEO::index_t index;
        if (face.size() == 3)
            index = out.facets.create_triangle((GEO::index_t)face[0], (GEO::index_t)face[1], (GEO::index_t)face[2]);
        else if (face.size() == 4)
            index = out.facets.create_quad((GEO::index_t)face[0], (GEO::index_t)face[1],
                                           (GEO::index_t)face[2], (GEO::index_t)face[3]);
        else
        {
            GEO::vector<GEO::index_t> polygon(face.size());
            for (size_t i = 0; i < face.size(); ++i) polygon[(GEO::index_t)i] = (GEO::index_t)face[i];
            index = out.facets.create_polygon(polygon);
        }
        facetHexFace[index] = (unsigned int)f;
        const bool hasTrajectory = (f < layout.hex_face_n_refined.size() && layout.hex_face_n_refined[f] > 0);
        facetNRefined[index] = asUnsigned(f < layout.hex_face_n_refined.size() ? layout.hex_face_n_refined[f] : 0);
        facetBorder[index] = hasTrajectory ? 0u : 1u;
    }

    const GEO::index_t nCells = (GEO::index_t)connectivity.cells.size();
    if (nCells > 0)
    {
        // geogram's hexahedra use a different vertex order than the extractor
        const std::unordered_set<uint64_t> edges = [&]() {
            std::unordered_set<uint64_t> result;
            result.reserve(connectivity.edges.size() * 2);
            for (size_t e = 0; e < connectivity.edges.size(); ++e)
                result.insert(edgeKey(connectivity.edges[e][0], connectivity.edges[e][1]));
            return result;
        }();

        const GEO::index_t first = out.cells.create_hexes(nCells);
        size_t nFallbacks = 0;
        for (GEO::index_t c = 0; c < nCells; ++c)
        {
            const std::vector<int>& cell = connectivity.cells[c];
            std::array<int, 8> order;
            if (cell.size() != 8 ||
                !reorderHexCell(std::array<int, 8>{{ cell[0], cell[1], cell[2], cell[3],
                                                     cell[4], cell[5], cell[6], cell[7] }},
                                edges, vertices, order))
            {
                ++nFallbacks;
                for (GEO::index_t i = 0; i < 8 && i < cell.size(); ++i)
                    out.cells.set_vertex(first + c, i, (GEO::index_t)cell[i]);
                continue;
            }
            for (GEO::index_t i = 0; i < 8; ++i)
                out.cells.set_vertex(first + c, i, (GEO::index_t)order[i]);
        }

        /*
         * Build the cell adjacency: geogram does not derive it automatically,
         * and the .geogram writer stores it (the loader skips connect() for
         * that format). Sliver removal is disabled: no cell may be removed,
         * because the attributes are indexed by cell.
         */
        out.cells.connect(false);
        if (nFallbacks > 0)
            std::cerr << "GeogramBridge: warning: " << nFallbacks
                      << " hexahedra could not be reordered into geogram's vertex convention" << std::endl;
        GEO::Attribute<unsigned int> cellHexCell(out.cells.attributes(), ATTR_HEX_CELL);
        GEO::Attribute<unsigned int> cellNRefined(out.cells.attributes(), ATTR_N_REFINED);
        GEO::Attribute<unsigned int> cellNPieces(out.cells.attributes(), ATTR_N_PIECES);
        for (GEO::index_t c = 0; c < nCells; ++c)
        {
            cellHexCell[first + c] = (unsigned int)c;
            cellNRefined[first + c] = asUnsigned(c < layout.hex_cell_n_refined.size()
                    ? layout.hex_cell_n_refined[c] : 0);
            cellNPieces[first + c] = asUnsigned(c < layout.hex_cell_pieces.size()
                    ? (int)layout.hex_cell_pieces[c].size() : 0);
        }
    }
}

bool GeogramBridge::save(const GEO::Mesh& mesh, const std::string& filename)
{
    initialize();
    return GEO::mesh_save(mesh, filename);
}

bool GeogramBridge::saveRefinedMesh(const HexTrajectoryLayout& layout, const std::string& filename,
                                    bool withCells)
{
    GEO::Mesh mesh;
    toRefinedMesh(layout, mesh, withCells);
    return save(mesh, filename);
}

bool GeogramBridge::saveHexMesh(const HexTrajectoryLayout& layout,
                                const std::vector<Vec3d>& vertices,
                                const TrajectoryMeshView& connectivity,
                                const std::string& filename)
{
    GEO::Mesh mesh;
    toHexMesh(layout, vertices, connectivity, mesh);
    return save(mesh, filename);
}

std::string GeogramBridge::describe(const GEO::Mesh& mesh)
{
    std::ostringstream os;
    os << "vertices: " << mesh.vertices.nb() << ", edges: " << mesh.edges.nb()
       << ", facets: " << mesh.facets.nb() << ", facet corners: " << mesh.facet_corners.nb()
       << ", cells: " << mesh.cells.nb() << "\n";

    struct StoreListing
    {
        static void append(std::ostream& os, const char* what, const GEO::AttributesManager& manager)
        {
            GEO::vector<std::string> names;
            manager.list_attribute_names(names);
            std::sort(names.begin(), names.end());
            os << "  " << what << ":";
            for (GEO::index_t i = 0; i < names.size(); ++i)
            {
                const std::string name = names[i];
                const GEO::AttributeStore* store = manager.find_attribute_store(name);
                os << " " << name;
                if (store) os << "[" << store->dimension() << "]";
            }
            os << "\n";
        }
    };

    StoreListing::append(os, "vertex attributes", mesh.vertices.attributes());
    StoreListing::append(os, "edge attributes", mesh.edges.attributes());
    StoreListing::append(os, "facet attributes", mesh.facets.attributes());
    StoreListing::append(os, "cell attributes", mesh.cells.attributes());

    return os.str();
}

} // namespace HexEx
