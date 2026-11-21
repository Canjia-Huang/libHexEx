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

/**
 * @file
 *
 * Extracts a hex mesh from an integer-grid map and writes the trajectories of
 * its vertices, edges, facets and cells through the input tet mesh as two
 * geogram meshes:
 *
 *   <out>_refined.geogram  the refined mesh (subdivision of the input tet mesh)
 *   <out>_hex.geogram      the extracted hex mesh
 *
 * Both meshes carry the correspondence as geogram attributes, see
 * HexEx::GeogramBridge. After writing, both files are read back and the
 * attributes are verified.
 */

#include <HexEx/HexEx.hh>
#include <HexEx/FileAccessor.hh>
#include <HexEx/HexExtractor.hh>
#include <HexEx/TrajectoryLayout.hh>

#include <hexex_geogram.hh>

#include <geogram/basic/attributes.h>
#include <geogram/mesh/mesh_io.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace HexEx;

namespace
{

int g_failures = 0;
int g_warnings = 0;

void check(bool condition, const std::string& what)
{
    if (condition)
        std::cout << "  [ ok ] " << what << std::endl;
    else
    {
        std::cout << "  [FAIL] " << what << std::endl;
        ++g_failures;
    }
}

void warn(bool condition, const std::string& what)
{
    if (condition)
        std::cout << "  [ ok ] " << what << std::endl;
    else
    {
        std::cout << "  [warn] " << what << std::endl;
        ++g_warnings;
    }
}

void printUsage(const char* progname)
{
    std::cout << "usage: " << progname << " <in.hexex> <out-prefix>\n"
              << "\n"
              << "Extracts a hex mesh from the integer-grid map in <in.hexex> and\n"
              << "writes the trajectories of its elements as geogram meshes:\n"
              << "  <out-prefix>_refined.geogram   refined mesh (subdivision of the tet mesh)\n"
              << "  <out-prefix>_hex.geogram       extracted hex mesh\n"
              << "\n"
              << "Set HEXEX_TRAJECTORY_DIAGNOSTICS=1 for extra diagnostics.\n";
}

/// verifies that reading the refined mesh back yields the expected attributes
void verifyRefinedRoundTrip(const std::string& filename, const HexTrajectoryLayout& layout)
{
    std::cout << "verifying " << filename << std::endl;

    GEO::Mesh mesh;
    if (!GEO::mesh_load(filename, mesh))
    {
        std::cout << "  [FAIL] could not read the file back" << std::endl;
        ++g_failures;
        return;
    }

    check(mesh.vertices.nb() == layout.vertex_position.size(), "vertex count");
    check(mesh.edges.nb() == layout.edge_vertices.size(), "edge count");
    check(mesh.facets.nb() == layout.facet_vertices.size(), "facet count");
    check(mesh.cells.nb() == layout.cell_vertices.size(), "cell count");
    check(mesh.cells.are_simplices(), "cells are tetrahedra");

    GEO::Attribute<int> vertexHexVertex(mesh.vertices.attributes(), "hex_vertex");
    GEO::Attribute<unsigned int> vertexGrid(mesh.vertices.attributes(), "grid");
    GEO::Attribute<int> edgeHexEdge(mesh.edges.attributes(), "hex_edge");
    GEO::Attribute<int> edgeStep(mesh.edges.attributes(), "hex_edge_step");
    GEO::Attribute<int> edgePiece(mesh.edges.attributes(), "piece");
    GEO::Attribute<int> edgeTetEdge(mesh.edges.attributes(), "tet_edge");
    GEO::Attribute<unsigned int> edgeGrid(mesh.edges.attributes(), "grid");
    GEO::Attribute<int> edgeGridAxis(mesh.edges.attributes(), "grid_axis");
    GEO::Attribute<int> facetHexFace(mesh.facets.attributes(), "hex_face");
    GEO::Attribute<int> facetPiece(mesh.facets.attributes(), "piece");
    GEO::Attribute<int> facetTetFace(mesh.facets.attributes(), "tet_face");
    GEO::Attribute<unsigned int> facetKind(mesh.facets.attributes(), "kind");
    GEO::Attribute<int> cellHexCell(mesh.cells.attributes(), "hex_cell");
    GEO::Attribute<unsigned int> cellPiece(mesh.cells.attributes(), "piece");

    check(vertexHexVertex.is_bound() && vertexGrid.is_bound(), "vertex attributes survived");
    check(edgeHexEdge.is_bound() && edgeStep.is_bound() && edgePiece.is_bound() && edgeTetEdge.is_bound() &&
          edgeGrid.is_bound() && edgeGridAxis.is_bound(),
          "edge attributes survived");
    check(facetHexFace.is_bound() && facetPiece.is_bound() && facetTetFace.is_bound() && facetKind.is_bound(),
          "facet attributes survived");
    check(cellHexCell.is_bound() && cellPiece.is_bound(), "cell attributes survived");
    if (g_failures) return;

    size_t mismatches = 0;
    for (GEO::index_t v = 0; v < mesh.vertices.nb(); ++v)
    {
        if (vertexHexVertex[v] != layout.vertex_hex_vertex[v]) ++mismatches;
        if ((vertexGrid[v] != 0) != (layout.vertex_is_grid[v] != 0)) ++mismatches;
    }
    for (GEO::index_t e = 0; e < mesh.edges.nb(); ++e)
    {
        if (edgeHexEdge[e] != layout.edge_hex_edge[e]) ++mismatches;
        if (edgeStep[e] != layout.edge_hex_edge_step[e]) ++mismatches;
        if (edgePiece[e] != layout.edge_piece[e]) ++mismatches;
        if (edgeTetEdge[e] != layout.edge_tet_edge[e]) ++mismatches;
        if ((edgeGrid[e] != 0) != (layout.edge_on_grid[e] != 0)) ++mismatches;
        if (edgeGridAxis[e] != layout.edge_grid_axis[e]) ++mismatches;
    }
    for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f)
    {
        if (facetHexFace[f] != layout.facet_hex_face[f]) ++mismatches;
        if (facetPiece[f] != layout.facet_piece_a[f]) ++mismatches;
        if (facetTetFace[f] != layout.facet_tet_face[f]) ++mismatches;
        if (facetKind[f] != (unsigned int)layout.facet_kind[f]) ++mismatches;
    }
    for (GEO::index_t c = 0; c < mesh.cells.nb(); ++c)
    {
        if (cellHexCell[c] != layout.cell_hex_cell[c]) ++mismatches;
        if (cellPiece[c] != (unsigned int)layout.cell_piece[c]) ++mismatches;
    }
    check(mismatches == 0, "attribute values match the layout");
    std::cout << GeogramBridge::describe(mesh);
}

/// verifies that reading the hex mesh back yields the expected attributes
void verifyHexRoundTrip(const std::string& filename,
                        const HexTrajectoryLayout& layout,
                        const HexahedralMesh& hexMesh)
{
    std::cout << "verifying " << filename << std::endl;

    GEO::Mesh mesh;
    if (!GEO::mesh_load(filename, mesh))
    {
        std::cout << "  [FAIL] could not read the file back" << std::endl;
        ++g_failures;
        return;
    }

    check(mesh.vertices.nb() == hexMesh.n_vertices(), "vertex count");
    check(mesh.edges.nb() == hexMesh.n_edges(), "edge count");
    check(mesh.facets.nb() == hexMesh.n_faces(), "facet count");
    check(mesh.cells.nb() == hexMesh.n_cells(), "cell count");

    GEO::Attribute<unsigned int> vertexNRefined(mesh.vertices.attributes(), "n_refined");
    GEO::Attribute<unsigned int> edgeNRefined(mesh.edges.attributes(), "n_refined");
    GEO::Attribute<unsigned int> facetNRefined(mesh.facets.attributes(), "n_refined");
    GEO::Attribute<unsigned int> cellNRefined(mesh.cells.attributes(), "n_refined");
    GEO::Attribute<unsigned int> cellNPieces(mesh.cells.attributes(), "n_pieces");

    check(vertexNRefined.is_bound(), "vertex n_refined survived");
    check(edgeNRefined.is_bound(), "edge n_refined survived");
    check(facetNRefined.is_bound(), "facet n_refined survived");
    check(cellNRefined.is_bound() && cellNPieces.is_bound(), "cell attributes survived");
    if (g_failures) return;

    size_t mismatches = 0;
    for (GEO::index_t v = 0; v < mesh.vertices.nb(); ++v)
        if (vertexNRefined[v] != (unsigned int)layout.hex_vertex_n_refined[v]) ++mismatches;
    for (GEO::index_t e = 0; e < mesh.edges.nb(); ++e)
        if (edgeNRefined[e] != (unsigned int)layout.hex_edge_n_refined[e]) ++mismatches;
    for (GEO::index_t f = 0; f < mesh.facets.nb(); ++f)
        if (facetNRefined[f] != (unsigned int)layout.hex_face_n_refined[f]) ++mismatches;
    for (GEO::index_t c = 0; c < mesh.cells.nb(); ++c)
    {
        if (cellNRefined[c] != (unsigned int)layout.hex_cell_n_refined[c]) ++mismatches;
        if (cellNPieces[c] != (unsigned int)layout.hex_cell_pieces[c].size()) ++mismatches;
    }
    check(mismatches == 0, "attribute values match the layout");
    std::cout << GeogramBridge::describe(mesh);
}

/// structural checks of the layout
void verifyLayout(const HexTrajectoryLayout& layout, const HexahedralMesh& hexMesh)
{
    std::cout << "checking the trajectory layout" << std::endl;

    // 1. the pieces partition the input tet mesh
    const double relativeVolumeError = std::fabs(layout.volume_pieces - layout.volume_input_tets) /
            std::max(1e-30, layout.volume_input_tets);
    check(relativeVolumeError < 1e-9, "the pieces partition the input tet mesh (relative volume error " +
            std::to_string(relativeVolumeError) + ")");

    // 2. every refined cell has a piece and belongs to a tet
    bool allCellsHavePieces = true;
    for (size_t c = 0; c < layout.cell_vertices.size(); ++c)
        if (layout.cell_piece[c] < 0 || layout.cell_tet_cell[c] < 0) allCellsHavePieces = false;
    check(allCellsHavePieces, "every refined tetrahedron has a piece and an input tet");

    // 3. refined element indices are consistent
    bool indicesOk = true;
    for (size_t i = 0; i < layout.vertex_hex_vertex.size(); ++i)
        if (layout.vertex_hex_vertex[i] >= (int)hexMesh.n_vertices()) indicesOk = false;
    for (size_t i = 0; i < layout.edge_hex_edge.size(); ++i)
        if (layout.edge_hex_edge[i] >= (int)hexMesh.n_edges()) indicesOk = false;
    for (size_t i = 0; i < layout.facet_hex_face.size(); ++i)
        if (layout.facet_hex_face[i] >= (int)hexMesh.n_faces()) indicesOk = false;
    for (size_t i = 0; i < layout.cell_hex_cell.size(); ++i)
        if (layout.cell_hex_cell[i] >= (int)hexMesh.n_cells()) indicesOk = false;
    check(indicesOk, "all hex element indices are in range");

    // 4. the refined mesh has no cracks
    check(layout.n_crack_edges == 0, "the refined mesh has no cracks (" +
            std::to_string(layout.n_crack_edges) + " crack edges)");

    // 5. the hex elements are covered by trajectories (holes are reported, not fatal)
    warn(layout.n_hex_vertices_without_trajectory == 0,
         "every hex vertex has trajectory vertices (" +
         std::to_string(layout.n_hex_vertices_without_trajectory) + " without)");
    warn(layout.n_hex_edges_without_trajectory == 0,
         "every hex edge has a polyline (" +
         std::to_string(layout.n_hex_edges_without_trajectory) + " without)");
    warn(layout.n_hex_faces_without_trajectory == 0,
         "every hex face has trajectory facets (" +
         std::to_string(layout.n_hex_faces_without_trajectory) + " without)");
    warn(layout.n_hex_cells_without_trajectory == 0,
         "every hex cell has trajectory tetrahedra (" +
         std::to_string(layout.n_hex_cells_without_trajectory) + " without)");

    // 6. every hex edge's polyline is ordered consistently
    size_t badPolylines = 0;
    for (size_t h = 0; h < layout.hex_edge_edges.size(); ++h)
    {
        const std::vector<int>& edges = layout.hex_edge_edges[h];
        for (size_t i = 0; i < edges.size(); ++i)
            if (layout.edge_hex_edge_step[edges[i]] != (int)i) ++badPolylines;
    }
    check(badPolylines == 0, "hex edge polylines are ordered (" +
            std::to_string(badPolylines) + " inconsistencies)");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        printUsage(argv[0]);
        return argc == 1 ? 0 : 2;
    }

    const std::string inFilename = argv[1];
    const std::string outPrefix = argv[2];

    TetrahedralMesh tetMesh;
    PerCellVertexProperty<Parameter> parameters(
            tetMesh.request_cell_property<HexEx::VertexMapProp<HexEx::Parameter>>());
    if (!HexEx::readFromFile(inFilename, tetMesh, parameters))
    {
        std::cerr << "could not read " << inFilename << std::endl;
        return 1;
    }

    std::cout << "extracting hex mesh from " << inFilename << std::endl;
    HexExtractor extractor(tetMesh, parameters);
    extractor.extract();

    HexahedralMesh hexMesh;
    HexTrajectoryLayout layout;
    extractor.getHexMeshWithTrajectory(hexMesh, layout);

    std::cout << "\nhex mesh: " << hexMesh.n_vertices() << " vertices, "
              << hexMesh.n_edges() << " edges, " << hexMesh.n_faces() << " faces, "
              << hexMesh.n_cells() << " cells" << std::endl;
    std::cout << "\n" << layout.describe() << std::endl;

    verifyLayout(layout, hexMesh);

    // the two geogram meshes
    std::vector<Vec3d> hexVertices;
    TrajectoryMeshView hexConnectivity;
    makeTrajectoryMesh(hexMesh, hexVertices, hexConnectivity);

    const std::string refinedFilename = outPrefix + "_refined.geogram";
    const std::string hexFilename = outPrefix + "_hex.geogram";

    std::cout << "writing " << refinedFilename << std::endl;
    if (!GeogramBridge::saveRefinedMesh(layout, refinedFilename))
    {
        std::cerr << "could not write " << refinedFilename << std::endl;
        return 1;
    }
    std::cout << "writing " << hexFilename << std::endl;
    if (!GeogramBridge::saveHexMesh(layout, hexVertices, hexConnectivity, hexFilename))
    {
        std::cerr << "could not write " << hexFilename << std::endl;
        return 1;
    }

    std::cout << "\n" << std::endl;
    verifyRefinedRoundTrip(refinedFilename, layout);
    std::cout << "\n" << std::endl;
    verifyHexRoundTrip(hexFilename, layout, hexMesh);

    std::cout << "\n" << (g_failures == 0 ? "All checks passed" : "CHECKS FAILED")
              << " (" << g_failures << " failures, " << g_warnings << " warnings)" << std::endl;

    return g_failures == 0 ? 0 : 1;
}
