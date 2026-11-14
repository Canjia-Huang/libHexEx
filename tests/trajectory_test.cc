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

#include <gtest/gtest.h>

#include "common.hh"

#include <HexExtractor.hh>
#include <FileAccessor.hh>
#include <TrajectoryLayout.hh>

#include <cmath>
#include <string>

using namespace HexEx;

namespace
{

/**
 * @brief Runs the extraction and the trajectory construction and checks the
 *        invariants that must always hold.
 *
 * @param requireFullCoverage whether every hex element must have trajectory
 *        elements (true for clean integer-grid maps)
 * @return the number of hex elements without trajectory elements
 */
size_t checkTrajectories(TetrahedralMesh& mesh,
                         PerCellVertexProperty<Parameter>& parametrization,
                         bool requireFullCoverage)
{
    EXPECT_TRUE(mesh.n_cells() > 0);

    HexExtractor extractor(mesh, parametrization);
    extractor.extract();

    HexahedralMesh hexMesh;
    HexTrajectoryLayout layout;
    extractor.getHexMeshWithTrajectory(hexMesh, layout);

    EXPECT_TRUE(hexMesh.n_cells() > 0);
    EXPECT_TRUE(layout.n_refined_vertices > 0);
    EXPECT_TRUE(layout.n_refined_edges > 0);
    EXPECT_TRUE(layout.n_refined_facets > 0);
    EXPECT_TRUE(layout.n_refined_cells > 0);
    EXPECT_TRUE(layout.n_pieces > 0);

    // the pieces partition the input tet mesh
    const double relativeVolumeError = std::fabs(layout.volume_pieces - layout.volume_input_tets) /
            std::max(1e-30, layout.volume_input_tets);
    EXPECT_LT(relativeVolumeError, 1e-9);

    // every refined tetrahedron belongs to a piece and to an input tet
    for (size_t c = 0; c < layout.cell_vertices.size(); ++c)
    {
        EXPECT_GE(layout.cell_piece[c], 0);
        EXPECT_GE(layout.cell_tet_cell[c], 0);
    }

    // the refined mesh has no cracks
    EXPECT_EQ(layout.n_crack_edges, 0u);
    EXPECT_EQ(layout.n_vertex_weld_conflicts, 0u);

    // all hex element indices are in range
    for (size_t v = 0; v < layout.vertex_hex_vertex.size(); ++v)
        EXPECT_LT(layout.vertex_hex_vertex[v], (int)hexMesh.n_vertices());
    for (size_t e = 0; e < layout.edge_hex_edge.size(); ++e)
        EXPECT_LT(layout.edge_hex_edge[e], (int)hexMesh.n_edges());
    for (size_t f = 0; f < layout.facet_hex_face.size(); ++f)
        EXPECT_LT(layout.facet_hex_face[f], (int)hexMesh.n_faces());
    for (size_t c = 0; c < layout.cell_hex_cell.size(); ++c)
        EXPECT_LT(layout.cell_hex_cell[c], (int)hexMesh.n_cells());

    // the polylines of the hex edges are ordered consistently
    for (size_t h = 0; h < layout.hex_edge_edges.size(); ++h)
    {
        const std::vector<int>& edges = layout.hex_edge_edges[h];
        for (size_t i = 0; i < edges.size(); ++i)
            EXPECT_EQ(layout.edge_hex_edge_step[edges[i]], (int)i);
    }

    const size_t withoutTrajectory =
            layout.n_hex_vertices_without_trajectory +
            layout.n_hex_edges_without_trajectory +
            layout.n_hex_faces_without_trajectory +
            layout.n_hex_cells_without_trajectory;
    if (requireFullCoverage)
        EXPECT_EQ(withoutTrajectory, 0u);

    return withoutTrajectory;
}

} // namespace

TEST(Trajectory, Cube)
{
    // a cube parametrized with an integer grid of side length 2 yields a
    // 2x2x2 hex mesh and a perfect trajectory for every element
    for (auto transition : {false, true})
    {
        TetrahedralMesh mesh;
        if (transition)
            createCubeWithTransition(mesh, 2.0);
        else
            createCube(mesh, 2.0);

        PerCellVertexProperty<Parameter> parametrization =
                mesh.template request_cell_property<HexEx::VertexMapProp<HexEx::Parameter>>("Parametrization");
        mesh.set_persistent(parametrization);

        const size_t without = checkTrajectories(mesh, parametrization, true);
        EXPECT_EQ(without, 0u);
    }
}

TEST(Trajectory, SphereData)
{
    // "sphere.hexex" is an extremely coarse, heavily folded integer-grid map;
    // its trajectories are complete except around the holes of the extractor
    TetrahedralMesh mesh;
    PerCellVertexProperty<Parameter> parametrization =
            mesh.template request_cell_property<HexEx::VertexMapProp<HexEx::Parameter>>();
    ASSERT_TRUE(HexEx::readFromFile("testdata/sphere.hexex", mesh, parametrization));

    const size_t withoutTrajectory = checkTrajectories(mesh, parametrization, false);
    EXPECT_LT(withoutTrajectory, 100u);
}
