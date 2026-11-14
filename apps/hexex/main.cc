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
 * the hexahedral mesh and writes it. The format of the hex mesh is chosen from
 * the output file name: `.ovm` is written with OpenVolumeMesh, every other
 * extension with geogram (.geogram, .mesh, .obj, ...). Optionally the
 * trajectories of the extracted elements - the refined mesh (the input tet mesh
 * cut along all integer grid planes) and the hex mesh - are written as geogram
 * meshes, where every piece of information is kept in mesh attributes.
 *
 * Command line:
 *
 *   hexex --in <input.hexex> --out <hex.ovm|hex.geogram>
 *         [--out-poly <poly.geogram>] [--check]
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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
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

/// the lower case extension of a path, including the dot
std::string lowerExtension(const std::string& path)
{
    std::string extension = std::filesystem::path(path).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(),
                   [](unsigned char c) { return (char)std::tolower(c); });
    return extension;
}

/// whether the file name asks for the OpenVolumeMesh format
bool isOvmFile(const std::string& path)
{
    return lowerExtension(path) == ".ovm";
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
            "edges, facets and cells through the input tet mesh, which are available as a "
            "polygonal mesh with all information stored in mesh attributes (see --out-poly).");
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
    std::string polyFile;

    bool verbose = false;
    bool simpleMergePosition = false;
    bool polySurfaceOnly = false;
    bool check = false;

    app.add_option("--in,-i", inputFile,
            "Tetrahedral mesh with a per-cell integer grid map (.hexex)")->required()->check(CLI::ExistingFile);
    app.add_option("--out,-o", outputFile,
            "Extracted hex mesh. The format is chosen from the file name: .ovm uses "
            "OpenVolumeMesh, every other extension is written with geogram. Only .geogram keeps "
            "the cells and the trajectory attributes (.mesh keeps the cells but not the "
            "attributes, .obj/.ply/.stl only the surface). A directory derives the file name.");
#ifdef HEXEX_HAVE_GEOGRAM
    app.add_option("--out-poly,--o-poly", polyFile,
            "The polygonal mesh - the trajectory of the extracted elements, i.e. the input tet "
            "mesh cut along all integer grid planes - in geogram format. Use .geogram to keep "
            "the tetrahedra/cells and all attributes (a directory derives the file name).");
#else
    app.add_option("--out-poly,--o-poly", polyFile,
            "The polygonal mesh (trajectory of the extracted elements) in geogram format "
            "(not available in this build)");
#endif
    app.add_flag("-v,--verbose", verbose, "Print the extractor's progress output");
    app.add_flag("--simple-merge-position", simpleMergePosition,
            "Use the unoptimized merge positions (faster, slightly different vertex positions)");
    app.add_flag("--poly-surface-only", polySurfaceOnly,
            "Omit the tetrahedra of the polygonal mesh, i.e. write only the trajectories of the "
            "hex vertices, edges and faces (smaller files)");
    app.add_flag("--check,--validate", check,
            "Check the trajectory layout (volume conservation, cracks, coverage) and report it");
    app.set_version_flag("-V,--version", std::string(HEXEX_VERSION));

    CLI11_PARSE(app, argc, argv);

#ifndef HEXEX_HAVE_GEOGRAM
    if (!polyFile.empty() || (!outputFile.empty() && !isOvmFile(outputFile)))
    {
        std::cerr << "This HexEx build has no geogram support: only the OpenVolumeMesh format "
                     "(.ovm) is available, and --out-poly is unavailable. Configure with "
                     "-DGEOGRAM_ROOT=<path> or -DHEXEX_WITH_GEOGRAM=ON."
                  << std::endl;
        return 2;
    }
#endif

    const std::string model = modelName(inputFile);
    outputFile = resolveOutput(outputFile, model, "_hex.ovm");
    polyFile = resolveOutput(polyFile, model, "_poly.geogram");

    if (outputFile.empty() && polyFile.empty() && !check)
    {
        std::cerr << "Nothing to do: pass --out, --out-poly or --check." << std::endl;
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
     * The trajectory layout (which the polygonal output, the geogram hex
     * output and --check are built from) is only computed if it was asked for.
     * ------------------------------------------------------------------
     */
    const bool needLayout = !polyFile.empty() || check || (!outputFile.empty() && !isOvmFile(outputFile));

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
        if (isOvmFile(outputFile))
        {
            OpenVolumeMesh::IO::FileManager fileManager;
            if (!fileManager.writeFile(outputFile, hexMesh))
            {
                std::cerr << "Could not write " << outputFile << std::endl;
                return 4;
            }
        }
#ifdef HEXEX_HAVE_GEOGRAM
        else
        {
            std::vector<Vec3d> vertices;
            TrajectoryMeshView connectivity;
            makeTrajectoryMesh(hexMesh, vertices, connectivity);
            if (!GeogramBridge::saveHexMesh(layout, vertices, connectivity, outputFile))
            {
                std::cerr << "Could not write " << outputFile << std::endl;
                return 4;
            }
        }
#endif
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
    if (!polyFile.empty())
    {
        if (!GeogramBridge::saveRefinedMesh(layout, polyFile, !polySurfaceOnly))
        {
            std::cerr << "Could not write " << polyFile << std::endl;
            return 4;
        }
        std::cout << "wrote " << polyFile << " (" << layout.n_pieces << " pieces, "
                  << (polySurfaceOnly ? 0 : layout.n_refined_cells) << " refined tetrahedra, "
                  << layout.n_refined_facets << " facets, all information in mesh attributes)"
                  << std::endl;
    }
#endif

    return 0;
}
