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
 * hexex - command line application for HexEx.
 *
 * Reads a tetrahedral mesh with a per-cell integer grid map (.hexex), extracts
 * the hexahedral mesh and writes it in OpenVolumeMesh format. Optionally it
 * also writes the trajectories of the extracted elements - the refined mesh
 * and the hex mesh in geogram format, where every piece of information is kept
 * in mesh attributes.
 *
 * Command line:
 *
 *   hexex --in <input.hexex> [--out <output.ovm>]
 *         [--out-refined <refined.geogram>] [--out-hex <hex.geogram>]
 *         [--out-trajectory <directory-or-prefix>]
 *
 * @see HexEx::HexTrajectoryLayout, HexEx::GeogramBridge
 */

#include <HexEx/HexEx.hh>
#include <HexEx/FileAccessor.hh>
#include <HexEx/HexExtractor.hh>
#include <HexEx/TrajectoryLayout.hh>
#include <HexEx/Verbose.hh>

#include <HexEx/Config/Version.hh>

#ifdef HEXEX_HAVE_GEOGRAM
#include <hexex_geogram.hh>
#endif

#include <OpenVolumeMesh/FileManager/FileManager.hh>

#include <CLI/CLI.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace HexEx;

namespace
{

/// "name.hexex" -> "name"
std::string modelName(const std::string& path)
{
    const std::filesystem::path p(path);
    const std::string filename = p.filename().string();
    const std::size_t pos = filename.find_last_of('.');
    return pos == std::string::npos ? filename : filename.substr(0, pos);
}

/**
 * @brief Turns an output argument into a file name.
 *
 * If @p path is an existing directory the file is named
 * `<path>/<model><suffix>`, otherwise @p path is used verbatim.
 */
std::string resolveOutput(const std::string& path, const std::string& model, const std::string& suffix)
{
    if (path.empty()) return path;
    if (std::filesystem::is_directory(path))
        return (std::filesystem::path(path) / (model + suffix)).string();
    return path;
}

/**
 * @brief The structural checks of the trajectory layout.
 *
 * @return the number of failures
 */
int checkLayout(const HexTrajectoryLayout& layout, const HexahedralMesh& hexMesh)
{
    int failures = 0;
    int warnings = 0;

    const auto report = [&](bool ok, const std::string& what, bool warningOnly) {
        std::cout << "  [" << (ok ? "ok  " : (warningOnly ? "warn" : "FAIL")) << "] " << what << std::endl;
        if (!ok) { if (warningOnly) ++warnings; else ++failures; }
    };

    const double relativeVolumeError = std::fabs(layout.volume_pieces - layout.volume_input_tets) /
            std::max(1e-30, layout.volume_input_tets);
    report(relativeVolumeError < 1e-9,
           "the pieces partition the input tet mesh (relative volume error " +
           std::to_string(relativeVolumeError) + ")", false);

    bool cellsOk = true;
    for (std::size_t c = 0; c < layout.cell_vertices.size(); ++c)
        if (layout.cell_piece[c] < 0 || layout.cell_tet_cell[c] < 0) cellsOk = false;
    report(cellsOk, "every refined tetrahedron has a piece and an input tet", false);

    bool indicesOk = true;
    for (std::size_t v = 0; v < layout.vertex_hex_vertex.size(); ++v)
        if (layout.vertex_hex_vertex[v] >= (int)hexMesh.n_vertices()) indicesOk = false;
    for (std::size_t e = 0; e < layout.edge_hex_edge.size(); ++e)
        if (layout.edge_hex_edge[e] >= (int)hexMesh.n_edges()) indicesOk = false;
    for (std::size_t f = 0; f < layout.facet_hex_face.size(); ++f)
        if (layout.facet_hex_face[f] >= (int)hexMesh.n_faces()) indicesOk = false;
    for (std::size_t c = 0; c < layout.cell_hex_cell.size(); ++c)
        if (layout.cell_hex_cell[c] >= (int)hexMesh.n_cells()) indicesOk = false;
    report(indicesOk, "all hex element indices are in range", false);

    report(layout.n_crack_edges == 0,
           "the refined mesh has no cracks (" + std::to_string(layout.n_crack_edges) + " crack edges)", false);

    report(layout.n_hex_vertices_without_trajectory == 0,
           "every hex vertex has trajectory vertices (" +
           std::to_string(layout.n_hex_vertices_without_trajectory) + " without)", true);
    report(layout.n_hex_edges_without_trajectory == 0,
           "every hex edge has a polyline (" +
           std::to_string(layout.n_hex_edges_without_trajectory) + " without)", true);
    report(layout.n_hex_faces_without_trajectory == 0,
           "every hex face has trajectory facets (" +
           std::to_string(layout.n_hex_faces_without_trajectory) + " without)", true);
    report(layout.n_hex_cells_without_trajectory == 0,
           "every hex cell has trajectory tetrahedra (" +
           std::to_string(layout.n_hex_cells_without_trajectory) + " without)", true);

    std::size_t badPolylines = 0;
    for (std::size_t h = 0; h < layout.hex_edge_edges.size(); ++h)
    {
        const std::vector<int>& edges = layout.hex_edge_edges[h];
        for (std::size_t i = 0; i < edges.size(); ++i)
            if (layout.edge_hex_edge_step[edges[i]] != (int)i) ++badPolylines;
    }
    report(badPolylines == 0,
           "the hex edge polylines are ordered (" + std::to_string(badPolylines) + " inconsistencies)", false);

    std::cout << "  " << failures << " failures, " << warnings << " warnings" << std::endl;
    return failures;
}

} // namespace

int main(int argc, char** argv)
{
    CLI::App app(
            "hexex - extract a hexahedral mesh from a tetrahedral mesh with an integer grid map\n"
            "Based on the original libHexEx command line tool (cmdline_tool): it extracts the hex "
            "mesh exactly as before and additionally records the trajectories of its vertices, "
            "edges, facets and cells through the input tet mesh, which can be written as geogram "
            "meshes with all information stored in mesh attributes (see --out-refined, --out-hex "
            "and --out-trajectory).");
    /*
     * Print the footer verbatim: CLI11's paragraph formatting would re-wrap it and
     * strip the indentation.
     */
    app.footer(
            "SUPPORT:\n"
            " - Developed by huangcanjia\n"
            " - For bug reports or requirements, please contact: Canjia Huang <huangcanjia0214@gmail.com>");
    if (auto* formatter = dynamic_cast<CLI::Formatter*>(app.get_formatter().get()))
        formatter->enable_footer_formatting(false);

    std::string inputFile;
    std::string outputFile;
    std::string refinedFile;
    std::string hexFile;
    std::string trajectoryPath;

    bool verbose = false;
    bool simpleMergePosition = false;
    bool surfaceOnly = false;
    bool check = false;

    app.add_option("--in,-i", inputFile,
            "Tetrahedral mesh with a per-cell integer grid map (.hexex)")->required()->check(CLI::ExistingFile);
    app.add_option("--out,-o", outputFile,
            "Extracted hex mesh in OpenVolumeMesh format (.ovm; a directory derives the file name)");
#ifdef HEXEX_HAVE_GEOGRAM
    app.add_option("--out-refined", refinedFile,
            "Refined mesh (the trajectory of the extracted elements: the input tet mesh cut along "
            "all integer grid planes) in geogram format with all information in mesh attributes");
    app.add_option("--out-hex", hexFile,
            "Extracted hex mesh in geogram format with the trajectory attributes");
    app.add_option("--out-trajectory", trajectoryPath,
            "Write both geogram meshes; a directory derives the file names, anything else is used "
            "as a file name prefix");
#else
    app.add_option("--out-refined", refinedFile,
            "Refined mesh in geogram format (not available in this build)");
    app.add_option("--out-hex", hexFile,
            "Extracted hex mesh in geogram format (not available in this build)");
    app.add_option("--out-trajectory", trajectoryPath,
            "Write both geogram meshes (not available in this build)");
#endif
    app.add_flag("-v,--verbose", verbose, "Print the extractor's progress output");
    app.add_flag("--simple-merge-position", simpleMergePosition,
            "Use the unoptimized merge positions (faster, slightly different vertex positions)");
    app.add_flag("--refined-surface-only", surfaceOnly,
            "Omit the refined tetrahedra, i.e. write only the trajectories of the hex vertices, "
            "edges and faces (smaller files)");
    app.add_flag("--check,--validate", check,
            "Check the trajectory layout (volume conservation, cracks, coverage) and report it");
    app.set_version_flag("-V,--version", std::string(HEXEX_VERSION));

    CLI11_PARSE(app, argc, argv);

#ifndef HEXEX_HAVE_GEOGRAM
    if (!refinedFile.empty() || !hexFile.empty() || !trajectoryPath.empty())
    {
        std::cerr << "This HexEx build has no geogram support: --out-refined, --out-hex and "
                     "--out-trajectory are unavailable. Configure with -DGEOGRAM_ROOT=<path> or "
                     "-DHEXEX_WITH_GEOGRAM=ON."
                  << std::endl;
        return 2;
    }
#endif

    const std::string model = modelName(inputFile);

    // --out-trajectory is a convenience for setting both geogram outputs
    if (!trajectoryPath.empty())
    {
        if (std::filesystem::is_directory(trajectoryPath))
        {
            if (refinedFile.empty())
                refinedFile = (std::filesystem::path(trajectoryPath) / (model + "_refined.geogram")).string();
            if (hexFile.empty())
                hexFile = (std::filesystem::path(trajectoryPath) / (model + "_hex.geogram")).string();
        }
        else
        {
            if (refinedFile.empty()) refinedFile = trajectoryPath + "_refined.geogram";
            if (hexFile.empty()) hexFile = trajectoryPath + "_hex.geogram";
        }
    }

    outputFile = resolveOutput(outputFile, model, "_hex.ovm");
    refinedFile = resolveOutput(refinedFile, model, "_refined.geogram");
    hexFile = resolveOutput(hexFile, model, "_hex.geogram");

    if (outputFile.empty() && refinedFile.empty() && hexFile.empty() && !check)
    {
        std::cerr << "Nothing to do: pass --out, --out-refined, --out-hex, --out-trajectory or --check."
                  << std::endl;
        return 2;
    }

    HexEx::setVerbose(verbose);

    /*
     * ------------------------------------------------------------------
     * Read the input
     * ------------------------------------------------------------------
     */
    TetrahedralMesh tetMesh;
    PerCellVertexProperty<Parameter> parameters(
            tetMesh.request_cell_property<HexEx::VertexMapProp<HexEx::Parameter>>());
    if (!HexEx::readFromFile(inputFile, tetMesh, parameters))
    {
        std::cerr << "Could not read " << inputFile << std::endl;
        return 3;
    }

    /*
     * ------------------------------------------------------------------
     * Extract
     *
     * The trajectory layout (which the geogram outputs and --check are built
     * from) is only computed if it was asked for.
     * ------------------------------------------------------------------
     */
    const bool needLayout = !refinedFile.empty() || !hexFile.empty() || check;

    HexExtractor extractor(tetMesh, parameters);
    extractor.extract();

    HexahedralMesh hexMesh;
    HexTrajectoryLayout layout;
    if (needLayout)
        extractor.getHexMeshWithTrajectory(hexMesh, layout, !simpleMergePosition);
    else
        extractor.getHexMesh(hexMesh, !simpleMergePosition);

    std::cout << "extracted " << hexMesh.n_cells() << " hex cells, " << hexMesh.n_faces()
              << " faces, " << hexMesh.n_edges() << " edges, " << hexMesh.n_vertices()
              << " vertices" << std::endl;

    /*
     * ------------------------------------------------------------------
     * Write the outputs
     * ------------------------------------------------------------------
     */
    if (!outputFile.empty())
    {
        OpenVolumeMesh::IO::FileManager fileManager;
        if (!fileManager.writeFile(outputFile, hexMesh))
        {
            std::cerr << "Could not write " << outputFile << std::endl;
            return 4;
        }
        std::cout << "wrote " << outputFile << " (" << hexMesh.n_cells() << " hex cells)" << std::endl;
    }

    if (needLayout)
        std::cout << "\n" << layout.describe() << std::endl;

    if (check)
    {
        std::cout << "checking the trajectory layout" << std::endl;
        if (checkLayout(layout, hexMesh) != 0)
        {
            std::cerr << "the trajectory layout does not satisfy the structural invariants" << std::endl;
            return 5;
        }
        std::cout << std::endl;
    }

#ifdef HEXEX_HAVE_GEOGRAM
    if (!refinedFile.empty())
    {
        if (!GeogramBridge::saveRefinedMesh(layout, refinedFile, !surfaceOnly))
        {
            std::cerr << "Could not write " << refinedFile << std::endl;
            return 4;
        }
        std::cout << "wrote " << refinedFile << " (" << layout.n_pieces << " pieces, "
                  << (surfaceOnly ? 0 : layout.n_refined_cells) << " refined tetrahedra, "
                  << layout.n_refined_facets << " facets, all information in mesh attributes)"
                  << std::endl;
    }

    if (!hexFile.empty())
    {
        std::vector<Vec3d> vertices;
        TrajectoryMeshView connectivity;
        makeTrajectoryMesh(hexMesh, vertices, connectivity);

        if (!GeogramBridge::saveHexMesh(layout, vertices, connectivity, hexFile))
        {
            std::cerr << "Could not write " << hexFile << std::endl;
            return 4;
        }
        std::cout << "wrote " << hexFile << " (" << hexMesh.n_cells()
                  << " hex cells with the trajectory attributes)" << std::endl;
    }
#endif

    return 0;
}
