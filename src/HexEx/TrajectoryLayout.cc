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

#include "HexExtractor.hh"
#include "TrajectoryLayout.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>

namespace HexEx
{

namespace
{

// ---------------------------------------------------------------------------
// tunables
// ---------------------------------------------------------------------------

/// Two refined vertices are the same refined vertex if they are this close
/// (in units of the grid, i.e. one grid cell has size 1).
const double WELD_TOL = 1e-7;

/// Pieces with a smaller volume are numerically empty and are dropped.
const double VOLUME_EPS = 1e-14;

/// A parameter space point lies on a plane if its distance is below this.
const double PLANE_EPS = 1e-9;

/// A refined grid vertex is matched to an extracted hex vertex if their
/// positions differ by less than this.
const double MATCH_TOL = 1e-6;

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

typedef std::array<int64_t, 3> Key3;

struct Key3Hash
{
    size_t operator()(const Key3& k) const
    {
        uint64_t h = (uint64_t)k[0] * 73856093ull;
        h ^= (uint64_t)k[1] * 19349663ull;
        h ^= (uint64_t)k[2] * 83492791ull;
        h ^= (uint64_t)k[0] >> 17;
        h ^= (uint64_t)k[1] >> 19;
        h ^= (uint64_t)k[2] >> 23;
        return (size_t)h;
    }
};

struct PairHash
{
    size_t operator()(const std::pair<int, int>& p) const
    {
        return (size_t)(((uint64_t)(uint32_t)p.first << 32) ^ (uint64_t)(uint32_t)p.second);
    }
};

struct TripleHash
{
    size_t operator()(const std::array<int, 3>& t) const
    {
        uint64_t h = (uint64_t)(uint32_t)t[0] * 73856093ull;
        h ^= (uint64_t)(uint32_t)t[1] * 19349663ull;
        h ^= (uint64_t)(uint32_t)t[2] * 83492791ull;
        return (size_t)h;
    }
};

inline Key3 quantize(const Vec3d& p, double tol)
{
    return Key3{{ (int64_t)std::llround(p[0] / tol),
                  (int64_t)std::llround(p[1] / tol),
                  (int64_t)std::llround(p[2] / tol) }};
}

inline int otherAxis(int a, int b) { return 3 - a - b; }

inline bool debugDiagnostics()
{
    static const bool enabled = (std::getenv("HEXEX_TRAJECTORY_DIAGNOSTICS") != nullptr);
    return enabled;
}

// ---------------------------------------------------------------------------
// A convex polyhedron in the parameter space of one input tet.
// ---------------------------------------------------------------------------

/**
 * @brief A plane in the parameter space of one input tet.
 *
 * `normal . x <= offset` is the halfspace the piece lives in, so `normal`
 * points out of the piece. `type == 0` means the plane is one of the tet's
 * four faces, `type == 1` that it is one of the six faces of the unit cube.
 * For cube planes `axis`/`side` identify the face of the piece's unit cube.
 */
struct PlaneInfo
{
    int type = 0;
    int face = -1;   ///< input mesh face index (type 0)
    int axis = -1;   ///< local grid axis of the cube face (type 1)
    int side = 0;    ///< 0 = the cube face with coordinate `cube[axis]`, 1 = `cube[axis] + 1`
    Vec3d normal;
    double offset = 0;
};

struct Polyhedron
{
    std::vector<Vec3d> vertices;          ///< in the tet's parameter space
    std::vector<std::vector<int> > faces; ///< vertex loops, oriented outwards
    std::vector<int> facePlane;           ///< index into the plane array

    void clear()
    {
        vertices.clear();
        faces.clear();
        facePlane.clear();
    }

    /**
     * @brief Volume of the polyhedron (faces oriented outwards).
     */
    double volume() const
    {
        double v = 0.0;
        for (size_t f = 0; f < faces.size(); ++f)
        {
            const std::vector<int>& face = faces[f];
            if (face.size() < 3) continue;
            const Vec3d& a = vertices[face[0]];
            for (size_t i = 1; i + 1 < face.size(); ++i)
                v += a | (vertices[face[i]] % vertices[face[i + 1]]);
        }
        return v / 6.0;
    }
};

/**
 * @brief Clips a convex polyhedron by the halfspace of @p plane.
 *
 * The cut polygon is added as a new face supported by the plane, so the result
 * stays closed and keeps its outward orientation.
 */
void clipPolyhedron(Polyhedron& poly, const PlaneInfo& plane, int planeIndex)
{
    const Vec3d n = plane.normal;
    const double d = plane.offset;
    const double eps = PLANE_EPS * (1.0 + std::fabs(d));

    std::vector<Vec3d> vertices = poly.vertices;
    std::unordered_map<Key3, int, Key3Hash> indexOf;
    indexOf.reserve(vertices.size() * 2);
    for (size_t i = 0; i < vertices.size(); ++i)
        indexOf.insert(std::make_pair(quantize(vertices[i], PLANE_EPS), (int)i));

    std::vector<std::vector<int> > faces;
    std::vector<int> facePlane;
    std::vector<int> cut;

    for (size_t f = 0; f < poly.faces.size(); ++f)
    {
        const std::vector<int>& face = poly.faces[f];
        const int m = (int)face.size();
        std::vector<int> out;
        for (int i = 0; i < m; ++i)
        {
            const int ia = face[i];
            const int ib = face[(i + 1) % m];
            const double da = (n | vertices[ia]) - d;
            const double db = (n | vertices[ib]) - d;
            const bool ina = da <= eps;
            const bool inb = db <= eps;
            if (ina) out.push_back(ia);
            if (ina != inb)
            {
                const double t = da / (da - db);
                const Vec3d p = vertices[ia] + t * (vertices[ib] - vertices[ia]);
                const Key3 key = quantize(p, PLANE_EPS);
                std::unordered_map<Key3, int, Key3Hash>::iterator it = indexOf.find(key);
                int id;
                if (it == indexOf.end())
                {
                    id = (int)vertices.size();
                    vertices.push_back(p);
                    indexOf.insert(std::make_pair(key, id));
                }
                else
                    id = it->second;
                out.push_back(id);
                cut.push_back(id);
            }
        }
        std::vector<int> cleaned;
        for (size_t i = 0; i < out.size(); ++i)
            if (cleaned.empty() || cleaned.back() != out[i])
                cleaned.push_back(out[i]);
        while (cleaned.size() > 1 && cleaned.front() == cleaned.back())
            cleaned.pop_back();
        if (cleaned.size() >= 3)
        {
            faces.push_back(cleaned);
            facePlane.push_back(poly.facePlane[f]);
        }
    }

    std::sort(cut.begin(), cut.end());
    cut.erase(std::unique(cut.begin(), cut.end()), cut.end());
    if (cut.size() >= 3)
    {
        // if the clipping plane coincides with a face of the polyhedron the cut
        // polygon is that face; do not duplicate it
        bool duplicate = false;
        for (size_t f = 0; f < faces.size() && !duplicate; ++f)
        {
            if (faces[f].size() != cut.size()) continue;
            std::vector<int> a = faces[f];
            std::sort(a.begin(), a.end());
            duplicate = (a == cut);
        }
        if (!duplicate)
        {
            // sort around the centroid such that the polygon's normal is the
            // outward normal of the clipping plane
            const Vec3d helper = (std::fabs(n[0]) < 0.9) ? Vec3d(1, 0, 0) : Vec3d(0, 1, 0);
            const Vec3d u = (helper % n).normalized();
            const Vec3d v = (n % u).normalized();
            Vec3d c(0, 0, 0);
            for (size_t i = 0; i < cut.size(); ++i) c = c + vertices[cut[i]];
            c = c * (1.0 / (double)cut.size());
            std::sort(cut.begin(), cut.end(), [&](int a, int b) {
                const Vec3d pa = vertices[a] - c;
                const Vec3d pb = vertices[b] - c;
                return std::atan2(pa | v, pa | u) < std::atan2(pb | v, pb | u);
            });
            faces.push_back(cut);
            facePlane.push_back(planeIndex);
        }
    }

    poly.vertices.swap(vertices);
    poly.faces.swap(faces);
    poly.facePlane.swap(facePlane);
}

/**
 * @brief Initializes @p poly with the tetrahedron whose vertices are `u[0..3]`.
 *
 * @param tetFaces the tet's four faces as local vertex index triples
 */
void initTetrahedron(Polyhedron& poly, const Vec3d u[4], const int tetFaces[4][3])
{
    poly.clear();
    for (int i = 0; i < 4; ++i) poly.vertices.push_back(u[i]);

    Vec3d center(0, 0, 0);
    for (int i = 0; i < 4; ++i) center = center + u[i];
    center = center * 0.25;

    for (int f = 0; f < 4; ++f)
    {
        std::vector<int> loop;
        loop.push_back(tetFaces[f][0]);
        loop.push_back(tetFaces[f][1]);
        loop.push_back(tetFaces[f][2]);
        const Vec3d a = u[loop[0]];
        const Vec3d b = u[loop[1]];
        const Vec3d c = u[loop[2]];
        const Vec3d n = (b - a) % (c - a);
        const Vec3d fc = (a + b + c) * (1.0 / 3.0);
        if ((n | (center - fc)) > 0) std::reverse(loop.begin(), loop.end());
        poly.faces.push_back(loop);
        poly.facePlane.push_back(f);
    }
}

/**
 * @brief Triangulates a convex polygon by a fan from its smallest vertex.
 *
 * The triangulation only depends on the *set* of vertices, so two pieces that
 * share a polygon produce exactly the same triangles (this is what makes the
 * refined mesh conforming).
 */
void canonicalTriangulation(const std::vector<int>& loop, std::vector<std::array<int, 3> >& triangles)
{
    triangles.clear();
    if (loop.size() < 3) return;

    size_t start = 0;
    for (size_t i = 1; i < loop.size(); ++i)
        if (loop[i] < loop[start]) start = i;

    std::vector<int> rotated;
    rotated.reserve(loop.size());
    rotated.push_back(loop[start]);
    for (size_t i = 1; i < loop.size(); ++i)
        rotated.push_back(loop[(start + i) % loop.size()]);

    for (size_t i = 1; i + 1 < rotated.size(); ++i)
    {
        std::array<int, 3> tri = {{ rotated[0], rotated[i], rotated[i + 1] }};
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]) continue;
        triangles.push_back(tri);
    }
}

/// per refined vertex data that is only needed while building the mesh
struct VertexAux
{
    /// the integer grid planes through the vertex, as (local axis, local coordinate)
    std::vector<std::array<int, 2> > gridPlanes;
};

} // namespace

// ---------------------------------------------------------------------------
// HexTrajectoryLayout
// ---------------------------------------------------------------------------

void HexTrajectoryLayout::clear()
{
    *this = HexTrajectoryLayout();
}

std::string HexTrajectoryLayout::describe() const
{
    std::ostringstream os;

    os << "input tetrahedra:                  " << n_input_tets
       << " (degenerate parametrization: " << n_degenerate_tets
       << ", filled without trajectory: " << n_filler_cells << ")\n";
    os << "pieces (tet x unit cube):          " << n_pieces
       << " (degenerate: " << n_degenerate_pieces
       << ", weld conflicts: " << n_vertex_weld_conflicts
       << ", bad tetrahedralization: " << n_pieces_with_bad_tetrahedralization << ")\n";
    os << "groups (unit cubes):               " << n_groups
       << " (without all corners: " << n_groups_without_all_corners
       << ", without hex cell: " << n_groups_without_hex_cell
       << ", duplicate cells: " << n_duplicate_hex_cells
       << ", corner conflicts: " << n_corner_conflicts << ")\n";
    os << "refined vertices:                  " << n_refined_vertices << "\n";
    os << "refined edges:                     " << n_refined_edges << "\n";
    os << "refined facets:                    " << n_refined_facets << "\n";
    os << "refined cells (tetrahedra):        " << n_refined_cells << "\n";
    os << "volume of the input tets:          " << volume_input_tets << "\n";
    os << "volume of the pieces:              " << volume_pieces
       << " (relative difference " << (volume_input_tets != 0.0
                ? (volume_pieces - volume_input_tets) / volume_input_tets : 0.0) << ")\n";
    os << "hex vertices without trajectory:   " << n_hex_vertices_without_trajectory << "\n";
    os << "hex edges without trajectory:      " << n_hex_edges_without_trajectory << "\n";
    os << "hex faces without trajectory:      " << n_hex_faces_without_trajectory << "\n";
    os << "hex cells without trajectory:      " << n_hex_cells_without_trajectory << "\n";
    os << "unmatched grid vertices:           " << n_unmatched_grid_vertices << "\n";
    os << "unresolved grid edges/facets:      " << n_unresolved_grid_edges
       << " / " << n_unresolved_grid_facets << "\n";
    os << "hex edges with branched polyline:  " << n_hex_edges_with_branched_polyline << "\n";
    os << "open refined edges:                " << n_open_edges
       << " (of which cracks: " << n_crack_edges << ")\n";

    return os.str();
}

// ---------------------------------------------------------------------------
// HexExtractor::buildTrajectoryLayout
// ---------------------------------------------------------------------------

void HexExtractor::buildTrajectoryLayout(HexTrajectoryLayout& layout, const TrajectoryMeshView& view)
{
    layout.clear();

    const size_t nCells = inputMesh.n_cells();
    layout.n_input_tets = nCells;

    const size_t nHexVertices = (size_t)view.n_vertices;
    const size_t nHexEdges = view.edges.size();
    const size_t nHexFaces = view.faces.size();
    const size_t nHexCells = view.cells.size();

    // -----------------------------------------------------------------------
    // 1. the extracted hex mesh's elements, indexed by their vertex sets
    // -----------------------------------------------------------------------
    std::map<std::vector<int>, int> finalCellByVertices;
    std::map<std::vector<int>, int> finalFaceByVertices;
    std::unordered_map<std::pair<int, int>, int, PairHash> finalEdgeByVertices;

    for (size_t c = 0; c < nHexCells; ++c)
    {
        if (view.cells[c].empty()) continue;
        std::vector<int> key = view.cells[c];
        std::sort(key.begin(), key.end());
        finalCellByVertices[key] = (int)c;
    }
    for (size_t f = 0; f < nHexFaces; ++f)
    {
        if (view.faces[f].empty()) continue;
        std::vector<int> key = view.faces[f];
        std::sort(key.begin(), key.end());
        finalFaceByVertices[key] = (int)f;
    }
    for (size_t e = 0; e < nHexEdges; ++e)
    {
        int a = view.edges[e][0];
        int b = view.edges[e][1];
        if (a < 0 || b < 0) continue;
        if (a > b) std::swap(a, b);
        finalEdgeByVertices[std::make_pair(a, b)] = (int)e;
    }

    // -----------------------------------------------------------------------
    // 2. a spatial index of the extracted (intermediate) hex vertices
    //
    //    A grid vertex of the refined mesh is matched to an extracted vertex by
    //    position. This is what makes the correspondence independent of a
    //    global integer grid frame.
    // -----------------------------------------------------------------------
    std::unordered_map<Key3, std::vector<int>, Key3Hash> intVertexBuckets;
    std::vector<Vec3d> intVertexPosition(intermediateHexMesh.n_vertices());
    for (auto vh : intermediateHexMesh.vertices())
    {
        const Vec3d p = intermediateHexMesh.vertex(vh);
        intVertexPosition[vh.idx()] = p;
        intVertexBuckets[quantize(p, MATCH_TOL)].push_back(vh.idx());
    }

    std::vector<int> intVertexHexVertex(intermediateHexMesh.n_vertices(), -1);
    for (auto vh : intermediateHexMesh.vertices())
    {
        const int cls = equivalenceClassIds[vh];
        if (cls < 0 || (size_t)cls >= equivalenceClassVertices.size()) continue;
        if (!equivalenceClassVertices[cls].is_valid()) continue;
        const int hexVertex = equivalenceClassVertices[cls].idx();
        if (hexVertex < 0 || (size_t)hexVertex >= nHexVertices) continue;
        intVertexHexVertex[vh.idx()] = hexVertex;
    }

    auto matchHexVertex = [&](const Vec3d& position) -> int
    {
        const Key3 key = quantize(position, MATCH_TOL);
        int best = -1;
        double bestDist2 = MATCH_TOL * MATCH_TOL;
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                {
                    const Key3 nk = Key3{{ key[0] + dx, key[1] + dy, key[2] + dz }};
                    const std::unordered_map<Key3, std::vector<int>, Key3Hash>::const_iterator it =
                            intVertexBuckets.find(nk);
                    if (it == intVertexBuckets.end()) continue;
                    for (size_t i = 0; i < it->second.size(); ++i)
                    {
                        const int v = it->second[i];
                        const Vec3d d = intVertexPosition[v] - position;
                        const double dist2 = d | d;
                        if (dist2 <= bestDist2) { bestDist2 = dist2; best = v; }
                    }
                }
        if (best < 0) return -1;
        return intVertexHexVertex[best];
    };

    // -----------------------------------------------------------------------
    // 3. the refined mesh
    // -----------------------------------------------------------------------
    std::unordered_map<Key3, std::vector<int>, Key3Hash> weldBuckets;
    std::vector<unsigned char> vertexHasThreeGridPlanes;

    auto findOrAddVertex = [&](const Vec3d& position, bool forceNew, bool& isNew) -> int
    {
        if (!forceNew)
        {
            const Key3 key = quantize(position, WELD_TOL);
            int best = -1;
            double bestDist2 = WELD_TOL * WELD_TOL;
            for (int dx = -1; dx <= 1; ++dx)
                for (int dy = -1; dy <= 1; ++dy)
                    for (int dz = -1; dz <= 1; ++dz)
                    {
                        const Key3 nk = Key3{{ key[0] + dx, key[1] + dy, key[2] + dz }};
                        const std::unordered_map<Key3, std::vector<int>, Key3Hash>::const_iterator it =
                                weldBuckets.find(nk);
                        if (it == weldBuckets.end()) continue;
                        for (size_t i = 0; i < it->second.size(); ++i)
                        {
                            const int v = it->second[i];
                            const Vec3d d = layout.vertex_position[v] - position;
                            const double dist2 = d | d;
                            if (dist2 <= bestDist2) { bestDist2 = dist2; best = v; }
                        }
                    }
            if (best >= 0) { isNew = false; return best; }
        }

        const int id = (int)layout.vertex_position.size();
        layout.vertex_position.push_back(position);
        layout.vertex_kind.push_back(0);
        layout.vertex_tet_cell.push_back(-1);
        layout.vertex_tet_face.push_back(-1);
        layout.vertex_tet_edge.push_back(-1);
        layout.vertex_tet_vertex.push_back(-1);
        layout.vertex_hex_vertex.push_back(-1);
        layout.vertex_is_grid.push_back(0);
        vertexHasThreeGridPlanes.push_back(0);
        weldBuckets[quantize(position, WELD_TOL)].push_back(id);
        isNew = true;
        return id;
    };

    // edges and facets of the refined mesh
    std::vector<std::array<int, 2> > edgeVertices;
    std::vector<int> edgePiece;
    std::vector<int> edgeFacetA;
    std::vector<int> edgeFacetB;
    std::unordered_map<std::pair<int, int>, int, PairHash> edgeIndex;
    std::unordered_map<std::array<int, 3>, int, TripleHash> facetIndex;

    // per facet: the cube plane it lies on (local axis / side), -1 if it is a tet face
    std::vector<int> facetCubePlaneAxis;
    std::vector<int> facetCubePlaneSide;

    // A (piece, welded vertex) pair and the piece-local grid planes through it.
    // The plane coordinates live in the piece's parameter frame, so they must
    // not be merged across pieces.
    struct LocalVertexInfo
    {
        int vertex;
        int nPlanes;
        std::array<std::array<int, 2>, 3> planes;
        std::array<int, 3> corner; ///< corner of the piece's cube, if nPlanes == 3
    };
    std::vector<std::vector<LocalVertexInfo> > pieceVertices;

    // where a refined edge lies on a grid line: the piece-local definition
    struct EdgeLocalInfo
    {
        int piece;
        int varying;
        int firstAxis, firstCoord;
        int secondAxis, secondCoord;
    };
    std::unordered_map<std::pair<int, int>, EdgeLocalInfo, PairHash> edgeLocalByPair;

    for (auto ch : inputMesh.cells())
    {
        const int ci = ch.idx();
        const std::vector<VertexHandle>& tetVertices = cellVertices[ch];
        if (tetVertices.size() != 4) continue;

        Vec3d p3[4];
        Vec3d u[4];
        for (int i = 0; i < 4; ++i)
        {
            p3[i] = inputPosition(tetVertices[i]);
            u[i] = parameter(ch, tetVertices[i]);
        }

        const double tetVolume3D = std::fabs((p3[1] - p3[3]) | ((p3[2] - p3[3]) % (p3[0] - p3[3]))) / 6.0;
        layout.volume_input_tets += tetVolume3D;

        // A tet whose parametrization is degenerate (its parameter space image
        // is flat) has no grid cells, but it still occupies volume in 3d. Fill
        // it with a single cell (and its domain boundary facets) so that the
        // complex stays a partition of the input mesh; it has no trajectory.
        const double tvParam = std::fabs((u[1] - u[3]) | ((u[2] - u[3]) % (u[0] - u[3]))) / 6.0;
        if (!(tvParam > 1e-9))
        {
            ++layout.n_degenerate_tets;
            if (tetVolume3D > VOLUME_EPS)
            {
                const int pieceIndex = (int)layout.pieces.size();
                HexTrajectoryLayout::Piece piece;
                piece.tet_cell = ci;
                piece.cube = std::array<int, 3>{{ -1, -1, -1 }};
                piece.volume = tetVolume3D;
                piece.group = -2; // not part of any unit cube
                piece.first_vertex = (int)layout.vertex_position.size();
                piece.first_cell = (int)layout.cell_vertices.size();
                piece.n_cells = 1;
                layout.pieces.push_back(piece);
                pieceVertices.push_back(std::vector<LocalVertexInfo>());
                layout.volume_pieces += tetVolume3D;
                ++layout.n_filler_cells;

                int ids[4];
                for (int i = 0; i < 4; ++i)
                {
                    bool isNew = false;
                    ids[i] = findOrAddVertex(p3[i], false, isNew);
                    if (isNew)
                    {
                        layout.vertex_kind[ids[i]] = 3;
                        layout.vertex_tet_vertex[ids[i]] = tetVertices[i].idx();
                        layout.vertex_tet_cell[ids[i]] = ci;
                    }
                }

                std::array<int, 4> cell = {{ ids[0], ids[1], ids[2], ids[3] }};
                const Vec3d a = layout.vertex_position[cell[0]];
                const Vec3d b = layout.vertex_position[cell[1]];
                const Vec3d c = layout.vertex_position[cell[2]];
                const Vec3d d = layout.vertex_position[cell[3]];
                if (((b - a) | ((c - a) % (d - a))) < 0) std::swap(cell[2], cell[3]);
                layout.cell_vertices.push_back(cell);
                layout.cell_piece.push_back(pieceIndex);
                layout.cell_tet_cell.push_back(ci);
                layout.cell_hex_cell.push_back(-1);

                // the domain boundary faces of the tet are not covered by any
                // piece, so add them to keep the surface closed
                for (auto hfh : inputMesh.cell(ch).halffaces())
                {
                    const FaceHandle fh = inputMesh.face_handle(hfh);
                    if (!inputMesh.is_boundary(fh)) continue;
                    const std::vector<VertexHandle> fvs = inputMesh.get_halfface_vertices(hfh);
                    if (fvs.size() != 3) continue;
                    int tri[3];
                    for (int i = 0; i < 3; ++i)
                    {
                        tri[i] = -1;
                        for (int k = 0; k < 4; ++k)
                            if (tetVertices[k] == fvs[i]) { tri[i] = ids[k]; break; }
                    }
                    if (tri[0] < 0 || tri[1] < 0 || tri[2] < 0) continue;
                    std::array<int, 3> key = {{ tri[0], tri[1], tri[2] }};
                    std::sort(key.begin(), key.end());
                    if (facetIndex.find(key) != facetIndex.end()) continue;
                    const int facet = (int)layout.facet_vertices.size();
                    facetIndex.insert(std::make_pair(key, facet));
                    layout.facet_vertices.push_back(std::array<int, 3>{{ tri[0], tri[1], tri[2] }});
                    layout.facet_hex_face.push_back(-1);
                    layout.facet_tet_face.push_back(fh.idx());
                    layout.facet_kind.push_back(0);
                    layout.facet_piece_a.push_back(pieceIndex);
                    layout.facet_piece_b.push_back(-1);
                    facetCubePlaneAxis.push_back(-1);
                    facetCubePlaneSide.push_back(-1);
                }
            }
            continue;
        }

        Matrix4x4d paramTo3D =
                getInverseParametrizationMatrix(p3[0], p3[1], p3[2], p3[3], u[0], u[1], u[2], u[3]);
        // ratio between the 3d volume and the parameter space volume
        const double volumeScale = std::fabs(
                (paramTo3D(0,0) * (paramTo3D(1,1) * paramTo3D(2,2) - paramTo3D(1,2) * paramTo3D(2,1))) -
                (paramTo3D(0,1) * (paramTo3D(1,0) * paramTo3D(2,2) - paramTo3D(1,2) * paramTo3D(2,0))) +
                (paramTo3D(0,2) * (paramTo3D(1,0) * paramTo3D(2,1) - paramTo3D(1,1) * paramTo3D(2,0))));

        // the tet's four faces (local vertex indices) and planes
        int tetFaces[4][3];
        std::array<PlaneInfo, 10> planes;
        {
            int f = 0;
            for (auto hfh : inputMesh.cell(ch).halffaces())
            {
                const std::vector<VertexHandle> fvs = inputMesh.get_halfface_vertices(hfh);
                if (fvs.size() != 3) { f = -1; break; }

                int li[3];
                int opp = -1;
                for (int i = 0; i < 3; ++i)
                {
                    li[i] = -1;
                    for (int k = 0; k < 4; ++k)
                        if (tetVertices[k] == fvs[i]) { li[i] = k; break; }
                    if (li[i] < 0) { f = -1; break; }
                }
                if (f < 0) break;
                for (int k = 0; k < 4; ++k)
                {
                    const bool inFace = (k == li[0] || k == li[1] || k == li[2]);
                    if (!inFace) { opp = k; break; }
                }
                if (opp < 0) { f = -1; break; }

                const Vec3d a = u[li[0]];
                const Vec3d b = u[li[1]];
                const Vec3d c = u[li[2]];
                Vec3d n = (b - a) % (c - a);
                const double len = n.length();
                if (!(len > 1e-14)) { f = -1; break; }
                n = n * (1.0 / len);
                if ((n | (u[opp] - a)) > 0) n = n * (-1.0);

                tetFaces[f][0] = li[0];
                tetFaces[f][1] = li[1];
                tetFaces[f][2] = li[2];
                planes[f].type = 0;
                planes[f].face = inputMesh.face_handle(hfh).idx();
                planes[f].normal = n;
                planes[f].offset = n | a;
                ++f;
            }
            if (f != 4) continue;
        }

        // the unit cubes of the integer grid that can intersect this tet
        Vec3d lo(1e300, 1e300, 1e300), hi(-1e300, -1e300, -1e300);
        for (int i = 0; i < 4; ++i)
            for (int k = 0; k < 3; ++k)
            {
                lo[k] = std::min(lo[k], u[i][k]);
                hi[k] = std::max(hi[k], u[i][k]);
            }
        const int cubeLo[3] = { (int)std::floor(lo[0] + 1e-9),
                                (int)std::floor(lo[1] + 1e-9),
                                (int)std::floor(lo[2] + 1e-9) };
        const int cubeHi[3] = { (int)std::ceil(hi[0] - 1e-9),
                                (int)std::ceil(hi[1] - 1e-9),
                                (int)std::ceil(hi[2] - 1e-9) };
        if (cubeHi[0] - cubeLo[0] > 256 || cubeHi[1] - cubeLo[1] > 256 || cubeHi[2] - cubeLo[2] > 256)
            continue; // pathological parametrization, do not attempt

        for (int I = cubeLo[0]; I < cubeHi[0]; ++I)
            for (int J = cubeLo[1]; J < cubeHi[1]; ++J)
                for (int K = cubeLo[2]; K < cubeHi[2]; ++K)
                {
                    const int cube[3] = { I, J, K };

                    // the six cube planes
                    for (int ax = 0; ax < 3; ++ax)
                        for (int side = 0; side < 2; ++side)
                        {
                            PlaneInfo pl;
                            pl.type = 1;
                            pl.axis = ax;
                            pl.side = side;
                            pl.face = cube[ax] + side;
                            Vec3d n(0, 0, 0);
                            n[ax] = side ? 1.0 : -1.0;
                            pl.normal = n;
                            pl.offset = n[ax] * (double)(cube[ax] + side);
                            planes[4 + 2 * ax + side] = pl;
                        }

                    Polyhedron poly;
                    initTetrahedron(poly, u, tetFaces);
                    for (int ip = 4; ip < 10; ++ip)
                        clipPolyhedron(poly, planes[ip], ip);

                    if (poly.faces.empty()) { ++layout.n_degenerate_pieces; continue; }
                    const double pieceParamVolume = poly.volume();
                    if (!(std::fabs(pieceParamVolume) > VOLUME_EPS))
                    {
                        ++layout.n_degenerate_pieces;
                        continue;
                    }
                    const double pieceVolume = std::fabs(pieceParamVolume) * volumeScale;

                    // only the vertices used by a face belong to the polyhedron
                    const int nPolyVertices = (int)poly.vertices.size();
                    std::vector<int> used;
                    std::vector<int> remap(nPolyVertices, -1);
                    for (size_t f = 0; f < poly.faces.size(); ++f)
                        for (size_t i = 0; i < poly.faces[f].size(); ++i)
                        {
                            const int local = poly.faces[f][i];
                            if (remap[local] < 0)
                            {
                                remap[local] = (int)used.size();
                                used.push_back(local);
                            }
                        }
                    const int nUsed = (int)used.size();
                    std::vector<int> welded(nPolyVertices, -1);

                    // ---- boundary vertices ----------------------------------
                    std::vector<int> kindOfVertex(nUsed, 0);
                    std::vector<int> tetFaceOfVertex(nUsed, -1);
                    std::vector<int> tetEdgeOfVertex(nUsed, -1);
                    std::vector<int> tetVertexOfVertex(nUsed, -1);
                    std::vector<LocalVertexInfo> localVertices(nUsed);
                    for (int v = 0; v < nUsed; ++v)
                    {
                        localVertices[v].vertex = -1;
                        localVertices[v].nPlanes = 0;
                    }

                    for (int v = 0; v < nUsed; ++v)
                    {
                        const Vec3d& local = poly.vertices[used[v]];
                        const Vec3d position = paramTo3D.transform_point(local);

                        int tetPlaneIndices[4];
                        int nTetPlanes = 0;
                        int gridPlaneIndices[6];
                        int nGridPlanes = 0;
                        for (int ip = 0; ip < 10; ++ip)
                        {
                            const double eps = PLANE_EPS * (1.0 + std::fabs(planes[ip].offset));
                            const double dist = std::fabs((planes[ip].normal | local) - planes[ip].offset);
                            if (dist > eps) continue;
                            if (ip < 4) { if (nTetPlanes < 4) tetPlaneIndices[nTetPlanes++] = ip; }
                            else { if (nGridPlanes < 6) gridPlaneIndices[nGridPlanes++] = ip; }
                        }

                        if (nTetPlanes >= 3)
                        {
                            kindOfVertex[v] = 3;
                            const std::vector<VertexHandle> f0v = inputMesh.get_halfface_vertices(
                                    inputMesh.halfface_handle(FaceHandle(planes[tetPlaneIndices[0]].face), 0));
                            for (size_t i = 0; i < f0v.size() && tetVertexOfVertex[v] < 0; ++i)
                            {
                                bool onAll = true;
                                for (int j = 1; j < 3; ++j)
                                {
                                    const std::vector<VertexHandle> fjv = inputMesh.get_halfface_vertices(
                                            inputMesh.halfface_handle(FaceHandle(planes[tetPlaneIndices[j]].face), 0));
                                    if (std::find(fjv.begin(), fjv.end(), f0v[i]) == fjv.end())
                                    {
                                        onAll = false;
                                        break;
                                    }
                                }
                                if (onAll) tetVertexOfVertex[v] = f0v[i].idx();
                            }
                        }
                        else if (nTetPlanes == 2)
                        {
                            kindOfVertex[v] = 2;
                            const std::vector<VertexHandle> f0v = inputMesh.get_halfface_vertices(
                                    inputMesh.halfface_handle(FaceHandle(planes[tetPlaneIndices[0]].face), 0));
                            const std::vector<VertexHandle> f1v = inputMesh.get_halfface_vertices(
                                    inputMesh.halfface_handle(FaceHandle(planes[tetPlaneIndices[1]].face), 0));
                            std::vector<VertexHandle> common;
                            for (size_t i = 0; i < f0v.size(); ++i)
                                if (std::find(f1v.begin(), f1v.end(), f0v[i]) != f1v.end())
                                    common.push_back(f0v[i]);
                            if (common.size() == 2)
                            {
                                const HalfEdgeHandle heh = inputMesh.find_halfedge(common[0], common[1]);
                                if (heh.is_valid()) tetEdgeOfVertex[v] = inputMesh.edge_handle(heh).idx();
                            }
                        }
                        else if (nTetPlanes == 1)
                        {
                            kindOfVertex[v] = 1;
                            tetFaceOfVertex[v] = planes[tetPlaneIndices[0]].face;
                        }

                        for (int i = 0; i < nGridPlanes; ++i)
                        {
                            const PlaneInfo& pl = planes[gridPlaneIndices[i]];
                            if (localVertices[v].nPlanes >= 3) break;
                            localVertices[v].planes[localVertices[v].nPlanes][0] = pl.axis;
                            localVertices[v].planes[localVertices[v].nPlanes][1] = cube[pl.axis] + pl.side;
                            ++localVertices[v].nPlanes;
                        }

                        bool isNew = false;
                        const int id = findOrAddVertex(position, false, isNew);
                        welded[used[v]] = id;

                        if (isNew)
                        {
                            layout.vertex_kind[id] = kindOfVertex[v];
                            layout.vertex_tet_cell[id] = ci;
                            layout.vertex_tet_face[id] = tetFaceOfVertex[v];
                            layout.vertex_tet_edge[id] = tetEdgeOfVertex[v];
                            layout.vertex_tet_vertex[id] = tetVertexOfVertex[v];
                        }
                        else
                        {
                            if (layout.vertex_tet_face[id] < 0) layout.vertex_tet_face[id] = tetFaceOfVertex[v];
                            if (layout.vertex_tet_edge[id] < 0) layout.vertex_tet_edge[id] = tetEdgeOfVertex[v];
                            if (layout.vertex_tet_vertex[id] < 0) layout.vertex_tet_vertex[id] = tetVertexOfVertex[v];
                            if (layout.vertex_tet_cell[id] < 0) layout.vertex_tet_cell[id] = ci;
                        }

                        localVertices[v].vertex = id;
                        if (nGridPlanes >= 3)
                        {
                            if ((size_t)id >= vertexHasThreeGridPlanes.size())
                                vertexHasThreeGridPlanes.resize((size_t)id + 1, 0);
                            vertexHasThreeGridPlanes[id] = 1;
                            layout.vertex_is_grid[id] = 1;
                        }
                    }

                    // two distinct boundary vertices of the same piece must not
                    // be welded together
                    bool weldConflict = false;
                    for (int a = 0; a < nUsed && !weldConflict; ++a)
                        for (int b = a + 1; b < nUsed; ++b)
                            if (welded[used[a]] == welded[used[b]]) { weldConflict = true; break; }
                    if (weldConflict) { ++layout.n_vertex_weld_conflicts; continue; }

                    const int pieceIndex = (int)layout.pieces.size();
                    HexTrajectoryLayout::Piece piece;
                    piece.tet_cell = ci;
                    piece.cube = std::array<int, 3>{{ I, J, K }};
                    piece.volume = pieceVolume;
                    piece.first_vertex = (int)layout.vertex_position.size();
                    piece.first_cell = (int)layout.cell_vertices.size();
                    layout.pieces.push_back(piece);
                    pieceVertices.push_back(localVertices);
                    layout.volume_pieces += pieceVolume;

                    // the centroid; the piece is tetrahedralized by coning from
                    // it to the triangles of the triangulated boundary
                    Vec3d localCentroid(0, 0, 0);
                    for (int v = 0; v < nUsed; ++v) localCentroid = localCentroid + poly.vertices[used[v]];
                    localCentroid = localCentroid * (1.0 / (double)nUsed);
                    bool centroidIsNew = false;
                    const int centroid = findOrAddVertex(paramTo3D.transform_point(localCentroid), true, centroidIsNew);
                    if (centroidIsNew)
                    {
                        layout.vertex_kind[centroid] = 0;
                        layout.vertex_tet_cell[centroid] = ci;
                    }

                    // ---- facets --------------------------------------------
                    std::unordered_map<int, int> localIndexOfVertex;
                    for (int v = 0; v < nUsed; ++v)
                        localIndexOfVertex[welded[used[v]]] = v;

                    std::vector<std::array<int, 3> > triangles;
                    for (size_t f = 0; f < poly.faces.size(); ++f)
                    {
                        const std::vector<int>& loopLocal = poly.faces[f];
                        std::vector<int> loop;
                        loop.reserve(loopLocal.size());
                        for (size_t i = 0; i < loopLocal.size(); ++i) loop.push_back(welded[loopLocal[i]]);
                        canonicalTriangulation(loop, triangles);

                        // remember on which grid line the polygon's edges lie
                        for (size_t i = 0; i < loop.size(); ++i)
                        {
                            const int va = loop[i];
                            const int vb = loop[(i + 1) % loop.size()];
                            if (va == vb) continue;
                            const std::unordered_map<int, int>::const_iterator ia = localIndexOfVertex.find(va);
                            const std::unordered_map<int, int>::const_iterator ib = localIndexOfVertex.find(vb);
                            if (ia == localIndexOfVertex.end() || ib == localIndexOfVertex.end()) continue;
                            const LocalVertexInfo& ta = localVertices[ia->second];
                            const LocalVertexInfo& tb = localVertices[ib->second];

                            int axes[3] = { -1, -1, -1 };
                            int coords[3] = { 0, 0, 0 };
                            int n = 0;
                            for (int j = 0; j < ta.nPlanes; ++j)
                                for (int k = 0; k < tb.nPlanes; ++k)
                                    if (ta.planes[j] == tb.planes[k])
                                    {
                                        if (n < 3) { axes[n] = ta.planes[j][0]; coords[n] = ta.planes[j][1]; ++n; }
                                    }
                            if (n < 2) continue;
                            const int firstAxis = axes[0];
                            int secondAxis = -1;
                            for (int j = 1; j < n; ++j)
                                if (axes[j] != firstAxis) { secondAxis = axes[j]; break; }
                            if (secondAxis < 0) continue;

                            int firstCoord = 0, secondCoord = 0;
                            for (int j = 0; j < n; ++j)
                            {
                                if (axes[j] == firstAxis) firstCoord = coords[j];
                                else if (axes[j] == secondAxis) secondCoord = coords[j];
                            }
                            const std::pair<int, int> key(std::min(va, vb), std::max(va, vb));
                            if (edgeLocalByPair.find(key) == edgeLocalByPair.end())
                            {
                                EdgeLocalInfo info;
                                info.piece = pieceIndex;
                                info.varying = otherAxis(firstAxis, secondAxis);
                                info.firstAxis = firstAxis;
                                info.firstCoord = firstCoord;
                                info.secondAxis = secondAxis;
                                info.secondCoord = secondCoord;
                                edgeLocalByPair.insert(std::make_pair(key, info));
                            }
                        }

                        const PlaneInfo& pl = planes[poly.facePlane[f]];
                        const int tetFace = (pl.type == 0) ? pl.face : -1;

                        // A facet can lie on a hex face even if its supporting
                        // polygon plane is a tet face: that happens when the
                        // domain boundary coincides with a grid plane. Detect it
                        // by looking for a cube plane shared by all vertices of
                        // the polygon (in this piece's frame).
                        int cubeAxis = (pl.type == 1) ? pl.axis : -1;
                        int cubeSide = (pl.type == 1) ? pl.side : -1;
                        int cubeCoord = (pl.type == 1) ? pl.face : 0;
                        (void)cubeCoord;
                        if (pl.type == 0 && !loop.empty())
                        {
                            const std::unordered_map<int, int>::const_iterator first =
                                    localIndexOfVertex.find(loop[0]);
                            if (first != localIndexOfVertex.end())
                            {
                                const LocalVertexInfo& firstInfo = localVertices[first->second];
                                for (int j = 0; j < firstInfo.nPlanes; ++j)
                                {
                                    const int axis = firstInfo.planes[j][0];
                                    const int coord = firstInfo.planes[j][1];
                                    bool all = true;
                                    for (size_t i = 1; i < loop.size() && all; ++i)
                                    {
                                        const std::unordered_map<int, int>::const_iterator it =
                                                localIndexOfVertex.find(loop[i]);
                                        if (it == localIndexOfVertex.end()) { all = false; break; }
                                        const LocalVertexInfo& info = localVertices[it->second];
                                        bool found = false;
                                        for (int k = 0; k < info.nPlanes; ++k)
                                            if (info.planes[k][0] == axis && info.planes[k][1] == coord)
                                                found = true;
                                        if (!found) all = false;
                                    }
                                    if (!all) continue;
                                    cubeAxis = axis;
                                    cubeCoord = coord;
                                    cubeSide = coord - cube[axis];
                                    break;
                                }
                            }
                        }

                        for (size_t t = 0; t < triangles.size(); ++t)
                        {
                            std::array<int, 3> tri = triangles[t];
                            std::array<int, 3> key = tri;
                            std::sort(key.begin(), key.end());
                            if (key[0] == key[1] || key[1] == key[2]) continue;

                            std::unordered_map<std::array<int, 3>, int, TripleHash>::iterator it =
                                    facetIndex.find(key);
                            int facet;
                            if (it == facetIndex.end())
                            {
                                facet = (int)layout.facet_vertices.size();
                                facetIndex.insert(std::make_pair(key, facet));
                                layout.facet_vertices.push_back(tri);
                                layout.facet_hex_face.push_back(-1);
                                layout.facet_tet_face.push_back(tetFace);
                                layout.facet_kind.push_back(cubeAxis >= 0 ? 1 : 0);
                                layout.facet_piece_a.push_back(pieceIndex);
                                layout.facet_piece_b.push_back(-1);
                                facetCubePlaneAxis.push_back(cubeAxis);
                                facetCubePlaneSide.push_back(cubeSide);
                            }
                            else
                            {
                                facet = it->second;
                                if (layout.facet_piece_a[facet] < 0) layout.facet_piece_a[facet] = pieceIndex;
                                else if (layout.facet_piece_b[facet] < 0) layout.facet_piece_b[facet] = pieceIndex;
                                if (layout.facet_tet_face[facet] < 0) layout.facet_tet_face[facet] = tetFace;
                                if (cubeAxis >= 0 && facetCubePlaneAxis[facet] < 0)
                                {
                                    facetCubePlaneAxis[facet] = cubeAxis;
                                    facetCubePlaneSide[facet] = cubeSide;
                                    layout.facet_kind[facet] = 1;
                                }
                            }
                        }
                    }

                    // ---- cells: cone from the centroid over every triangle --
                    double tetrahedralizedVolume = 0.0;
                    for (size_t f = 0; f < poly.faces.size(); ++f)
                    {
                        const std::vector<int>& loopLocal = poly.faces[f];
                        std::vector<int> loop;
                        loop.reserve(loopLocal.size());
                        for (size_t i = 0; i < loopLocal.size(); ++i) loop.push_back(welded[loopLocal[i]]);
                        canonicalTriangulation(loop, triangles);
                        for (size_t t = 0; t < triangles.size(); ++t)
                        {
                            std::array<int, 3> tri = triangles[t];
                            if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]) continue;
                            std::array<int, 4> tet = {{ centroid, tri[0], tri[1], tri[2] }};
                            const Vec3d a = layout.vertex_position[tet[0]];
                            const Vec3d b = layout.vertex_position[tet[1]];
                            const Vec3d c = layout.vertex_position[tet[2]];
                            const Vec3d d = layout.vertex_position[tet[3]];
                            const double signedVolume = (b - a) | ((c - a) % (d - a));
                            if (signedVolume < 0) std::swap(tet[2], tet[3]);
                            tetrahedralizedVolume += std::fabs(signedVolume) / 6.0;

                            layout.cell_vertices.push_back(tet);
                            layout.cell_piece.push_back(pieceIndex);
                            layout.cell_tet_cell.push_back(ci);
                            layout.cell_hex_cell.push_back(-1);
                        }
                    }

                    layout.pieces[pieceIndex].n_cells =
                            (int)layout.cell_vertices.size() - layout.pieces[pieceIndex].first_cell;
                    // slivers can have a relatively large error but it is
                    // absolutely negligible compared to the input tet
                    if (std::fabs(tetrahedralizedVolume - pieceVolume) > 1e-6 * pieceVolume &&
                        std::fabs(tetrahedralizedVolume - pieceVolume) > 1e-9 * tetVolume3D)
                        ++layout.n_pieces_with_bad_tetrahedralization;
                }
    }

    // -----------------------------------------------------------------------
    // 4. groups: all pieces of one unit cube of the integer grid
    //
    //    Pieces are connected across the facets that lie on a tet face. They
    //    belong to the same cube iff the neighbouring piece's cube equals the
    //    image of this piece's cube under the face's transition function. This
    //    uses the transition functions only *per face*, so it does not require
    //    a globally consistent integer grid frame.
    // -----------------------------------------------------------------------
    const int nPieces = (int)layout.pieces.size();

    struct Link
    {
        int piece;
        int face; ///< input face index
    };
    std::vector<std::vector<Link> > links(nPieces);
    for (size_t f = 0; f < layout.facet_vertices.size(); ++f)
    {
        const int a = layout.facet_piece_a[f];
        const int b = layout.facet_piece_b[f];
        if (a < 0 || b < 0 || a == b) continue;
        if (facetCubePlaneAxis[f] >= 0) continue; // a hex face, not a tet face
        const int face = layout.facet_tet_face[f];
        if (face < 0) continue;
        Link la; la.piece = b; la.face = face;
        Link lb; lb.piece = a; lb.face = face;
        links[a].push_back(la);
        links[b].push_back(lb);
    }

    std::vector<GridIsomorphism> pieceToCanonical(nPieces);
    std::vector<unsigned char> grouped(nPieces, 0);
    size_t nLinksTried = 0, nLinksRejected = 0;

    for (int seed = 0; seed < nPieces; ++seed)
    {
        if (grouped[seed]) continue;
        if (layout.pieces[seed].group == -2) { grouped[seed] = 1; continue; }
        const int groupIndex = (int)layout.groups.size();
        HexTrajectoryLayout::Group group;
        group.corners.fill(-1);
        layout.groups.push_back(group);

        std::vector<int> queue;
        queue.push_back(seed);
        grouped[seed] = 1;
        pieceToCanonical[seed] = GridIsomorphism(0);
        layout.pieces[seed].group = groupIndex;

        for (size_t qi = 0; qi < queue.size(); ++qi)
        {
            const int p = queue[qi];
            layout.groups[groupIndex].pieces.push_back(p);
            const int ch = layout.pieces[p].tet_cell;

            for (size_t l = 0; l < links[p].size(); ++l)
            {
                const int q = links[p][l].piece;
                if (grouped[q]) continue;

                const int face = links[p][l].face;
                HalfFaceHandle hfh = inputMesh.halfface_handle(FaceHandle(face), 0);
                if (!hfh.is_valid() || inputMesh.incident_cell(hfh) != CellHandle(ch))
                    hfh = inputMesh.opposite_halfface_handle(hfh);
                if (!hfh.is_valid() || inputMesh.incident_cell(hfh) != CellHandle(ch))
                {
                    if (debugDiagnostics())
                        std::cerr << "[trajectory] could not find the halfface of tet " << ch
                                  << " on face " << face << std::endl;
                    continue;
                }

                ++nLinksTried;
                GridIsomorphism transition = getTransitionFunction(hfh); // piece frame -> neighbour frame
                GridIsomorphism inverse = transition.inverted();

                // the neighbouring piece's cube is the image of this piece's
                // cube under the transition
                Vec3d expectedMin(1e300, 1e300, 1e300);
                for (int di = 0; di < 2; ++di)
                    for (int dj = 0; dj < 2; ++dj)
                        for (int dk = 0; dk < 2; ++dk)
                        {
                            const Vec3d corner = transition.transform_point(Vec3d(
                                    (double)(layout.pieces[p].cube[0] + di),
                                    (double)(layout.pieces[p].cube[1] + dj),
                                    (double)(layout.pieces[p].cube[2] + dk)));
                            for (int k = 0; k < 3; ++k) expectedMin[k] = std::min(expectedMin[k], corner[k]);
                        }
                const std::array<int, 3> expectedCube =
                        std::array<int, 3>{{ (int)std::llround(expectedMin[0]),
                                             (int)std::llround(expectedMin[1]),
                                             (int)std::llround(expectedMin[2]) }};
                if (expectedCube != layout.pieces[q].cube)
                {
                    // the two pieces share a tet face but not a cube; do not link
                    ++nLinksRejected;
                    if (debugDiagnostics() && nLinksRejected <= 4)
                    {
                        // report the area of the shared facet
                        double area = 0;
                        const int ff = links[p][l].face;
                        for (size_t k = 0; k < layout.facet_vertices.size(); ++k)
                            if (layout.facet_tet_face[k] == ff &&
                                (layout.facet_piece_a[k] == p || layout.facet_piece_b[k] == p) &&
                                (layout.facet_piece_a[k] == q || layout.facet_piece_b[k] == q))
                            {
                                const Vec3d& A = layout.vertex_position[layout.facet_vertices[k][0]];
                                const Vec3d& B = layout.vertex_position[layout.facet_vertices[k][1]];
                                const Vec3d& C = layout.vertex_position[layout.facet_vertices[k][2]];
                                area = std::max(area, ((B - A) % (C - A)).length());
                            }
                        const HalfFaceHandle h0 = inputMesh.halfface_handle(FaceHandle(ff), 0);
                        const HalfFaceHandle h1 = inputMesh.opposite_halfface_handle(h0);
                        std::cerr << "[trajectory] rejected link area " << area
                                  << " p.tet " << layout.pieces[p].tet_cell
                                  << " q.tet " << layout.pieces[q].tet_cell
                                  << " face " << ff
                                  << " face cells " << inputMesh.incident_cell(h0).idx()
                                  << "/" << inputMesh.incident_cell(h1).idx()
                                  << " ch " << ch << std::endl;
                    }
                    if (debugDiagnostics())
                        std::cerr << "[trajectory] pieces " << p << " (cube " << layout.pieces[p].cube[0]
                                  << "," << layout.pieces[p].cube[1] << "," << layout.pieces[p].cube[2]
                                  << ") and " << q << " (cube " << layout.pieces[q].cube[0]
                                  << "," << layout.pieces[q].cube[1] << "," << layout.pieces[q].cube[2]
                                  << ") share a tet face but not a cube; expected ("
                                  << expectedCube[0] << "," << expectedCube[1] << "," << expectedCube[2] << ")"
                                  << std::endl;
                    continue;
                }

                GridIsomorphism canonicalToP = pieceToCanonical[p];
                pieceToCanonical[q] = canonicalToP * inverse;
                grouped[q] = 1;
                layout.pieces[q].group = groupIndex;
                queue.push_back(q);
            }
        }
        layout.n_groups = layout.groups.size();
    }
    if (debugDiagnostics())
        std::cerr << "[trajectory] group links: tried " << nLinksTried
                  << ", rejected (different cubes) " << nLinksRejected << std::endl;

    // -----------------------------------------------------------------------
    // 5. the hex vertex of every grid vertex
    // -----------------------------------------------------------------------
    for (size_t v = 0; v < layout.vertex_position.size(); ++v)
    {
        if (v >= vertexHasThreeGridPlanes.size() || !vertexHasThreeGridPlanes[v]) continue;
        const int hexVertex = matchHexVertex(layout.vertex_position[v]);
        if (hexVertex < 0)
        {
            ++layout.n_unmatched_grid_vertices;
            continue;
        }
        layout.vertex_hex_vertex[v] = hexVertex;
    }

    // -----------------------------------------------------------------------
    // 6. the eight corner hex vertices of every group
    // -----------------------------------------------------------------------
    std::set<int> claimedHexCells;
    for (size_t g = 0; g < layout.groups.size(); ++g)
    {
        HexTrajectoryLayout::Group& group = layout.groups[g];
        if (group.pieces.empty()) continue;
        const std::array<int, 3> rootCube = layout.pieces[group.pieces[0]].cube;

        for (size_t pi = 0; pi < group.pieces.size(); ++pi)
        {
            const int p = group.pieces[pi];
            const HexTrajectoryLayout::Piece& piece = layout.pieces[p];
            const GridIsomorphism& toCanonical = pieceToCanonical[p];

            for (size_t i = 0; i < pieceVertices[p].size(); ++i)
            {
                const LocalVertexInfo& info = pieceVertices[p][i];
                if (info.nPlanes != 3 || info.vertex < 0) continue;

                int local[3];
                bool isCorner = true;
                for (int k = 0; k < 3; ++k)
                {
                    bool found = false;
                    for (int j = 0; j < info.nPlanes; ++j)
                        if (info.planes[j][0] == k)
                        {
                            const int c = info.planes[j][1] - piece.cube[k];
                            if (c != 0 && c != 1) { isCorner = false; break; }
                            local[k] = c;
                            found = true;
                            break;
                        }
                    if (!found || !isCorner) { isCorner = false; break; }
                }
                if (!isCorner) continue;

                const Vec3d canonical = toCanonical.transform_point(Vec3d(
                        (double)(piece.cube[0] + local[0]),
                        (double)(piece.cube[1] + local[1]),
                        (double)(piece.cube[2] + local[2])));
                const int cx = (int)std::llround(canonical[0] - (double)rootCube[0]);
                const int cy = (int)std::llround(canonical[1] - (double)rootCube[1]);
                const int cz = (int)std::llround(canonical[2] - (double)rootCube[2]);
                if (cx < 0 || cx > 1 || cy < 0 || cy > 1 || cz < 0 || cz > 1)
                {
                    ++layout.n_corner_conflicts;
                    continue;
                }
                const int corner = cx + 2 * cy + 4 * cz;
                const int hexVertex = layout.vertex_hex_vertex[info.vertex];
                if (hexVertex < 0) continue;
                if (group.corners[corner] < 0)
                    group.corners[corner] = hexVertex;
                else if (group.corners[corner] != hexVertex)
                    ++layout.n_corner_conflicts;
            }
        }

        bool complete = true;
        for (int c = 0; c < 8; ++c)
            if (group.corners[c] < 0) complete = false;
        if (!complete)
        {
            ++layout.n_groups_without_all_corners;
            if (debugDiagnostics() && layout.n_groups_without_all_corners <= 20)
            {
                int n = 0;
                for (int c = 0; c < 8; ++c) if (group.corners[c] >= 0) ++n;
                std::cerr << "[trajectory] group " << g << " of " << group.pieces.size()
                          << " pieces has " << n << "/8 corners" << std::endl;
            }
            continue;
        }

        std::vector<int> vertexSet(group.corners.begin(), group.corners.end());
        std::sort(vertexSet.begin(), vertexSet.end());
        const std::map<std::vector<int>, int>::const_iterator it = finalCellByVertices.find(vertexSet);
        if (it == finalCellByVertices.end())
        {
            ++layout.n_groups_without_hex_cell;
            continue;
        }
        group.hex_cell = it->second;
        if (claimedHexCells.find(group.hex_cell) != claimedHexCells.end())
            ++layout.n_duplicate_hex_cells;
        claimedHexCells.insert(group.hex_cell);
        for (size_t pi = 0; pi < group.pieces.size(); ++pi)
            layout.pieces[group.pieces[pi]].hex_cell = group.hex_cell;
    }

    // -----------------------------------------------------------------------
    // 7. hex faces: a refined facet on a cube plane belongs to a hex face
    // -----------------------------------------------------------------------
    for (size_t f = 0; f < layout.facet_vertices.size(); ++f)
    {
        const int ax = facetCubePlaneAxis[f];
        if (ax < 0) continue;
        const int side = facetCubePlaneSide[f];
        const int p = layout.facet_piece_a[f];
        if (p < 0) continue;
        const HexTrajectoryLayout::Piece& piece = layout.pieces[p];
        const int g = piece.group;
        if (g < 0) continue;
        HexTrajectoryLayout::Group& group = layout.groups[g];
        if (group.hex_cell < 0) continue;

        const std::array<int, 3> rootCube = layout.pieces[group.pieces[0]].cube;
        const GridIsomorphism& toCanonical = pieceToCanonical[p];
        const int a0 = (ax + 1) % 3;
        const int a1 = (ax + 2) % 3;

        int cornerHex[4];
        bool complete = true;
        int n = 0;
        for (int u = 0; u < 2 && complete; ++u)
            for (int w = 0; w < 2 && complete; ++w)
            {
                int local[3];
                local[ax] = side;
                local[a0] = u;
                local[a1] = w;
                const Vec3d canonical = toCanonical.transform_point(Vec3d(
                        (double)(piece.cube[0] + local[0]),
                        (double)(piece.cube[1] + local[1]),
                        (double)(piece.cube[2] + local[2])));
                const int cx = (int)std::llround(canonical[0] - (double)rootCube[0]);
                const int cy = (int)std::llround(canonical[1] - (double)rootCube[1]);
                const int cz = (int)std::llround(canonical[2] - (double)rootCube[2]);
                if (cx < 0 || cx > 1 || cy < 0 || cy > 1 || cz < 0 || cz > 1)
                {
                    complete = false;
                    break;
                }
                cornerHex[n] = group.corners[cx + 2 * cy + 4 * cz];
                if (cornerHex[n] < 0) complete = false;
                ++n;
            }
        if (!complete) { ++layout.n_unresolved_grid_facets; continue; }

        std::vector<int> faceVertices(cornerHex, cornerHex + 4);
        std::sort(faceVertices.begin(), faceVertices.end());
        const std::map<std::vector<int>, int>::const_iterator it = finalFaceByVertices.find(faceVertices);
        if (it == finalFaceByVertices.end())
        {
            ++layout.n_unresolved_grid_facets;
            continue;
        }
        layout.facet_hex_face[f] = it->second;
        layout.facet_kind[f] = 1;
    }

    // -----------------------------------------------------------------------
    // 8. hex edges: a refined edge on a grid line belongs to a hex edge
    // -----------------------------------------------------------------------
    auto addEdge = [&](int a, int b) -> int
    {
        if (a == b) return -1;
        if (a > b) std::swap(a, b);
        const std::pair<int, int> key(a, b);
        const std::unordered_map<std::pair<int, int>, int, PairHash>::const_iterator it = edgeIndex.find(key);
        if (it != edgeIndex.end()) return it->second;
        const int id = (int)edgeVertices.size();
        edgeVertices.push_back(std::array<int, 2>{{ a, b }});
        edgePiece.push_back(-1);
        edgeFacetA.push_back(-1);
        edgeFacetB.push_back(-1);
        edgeIndex.insert(std::make_pair(key, id));
        return id;
    };

    for (size_t f = 0; f < layout.facet_vertices.size(); ++f)
    {
        const std::array<int, 3>& tri = layout.facet_vertices[f];
        for (int i = 0; i < 3; ++i)
        {
            const int e = addEdge(tri[i], tri[(i + 1) % 3]);
            if (e < 0) continue;
            if (edgeFacetA[e] < 0) edgeFacetA[e] = (int)f;
            else if (edgeFacetB[e] < 0) edgeFacetB[e] = (int)f;
            if (edgePiece[e] < 0) edgePiece[e] = layout.facet_piece_a[f];
        }
    }

    layout.edge_vertices = edgeVertices;
    layout.edge_hex_edge.assign(edgeVertices.size(), -1);
    layout.edge_hex_edge_step.assign(edgeVertices.size(), -1);
    layout.edge_tet_edge.assign(edgeVertices.size(), -1);
    layout.edge_piece = edgePiece;
    layout.edge_on_grid.assign(edgeVertices.size(), 0);

    for (size_t e = 0; e < edgeVertices.size(); ++e)
    {
        const int a = edgeVertices[e][0];
        const int b = edgeVertices[e][1];
        const std::pair<int, int> pairKey(std::min(a, b), std::max(a, b));
        const std::unordered_map<std::pair<int, int>, EdgeLocalInfo, PairHash>::const_iterator localIt =
                edgeLocalByPair.find(pairKey);
        if (localIt == edgeLocalByPair.end()) continue;

        const EdgeLocalInfo& info = localIt->second;
        layout.edge_on_grid[e] = 1;
        layout.edge_piece[e] = info.piece;

        const HexTrajectoryLayout::Piece& piece = layout.pieces[info.piece];
        const int g = piece.group;
        if (g < 0) { ++layout.n_unresolved_grid_edges; continue; }
        HexTrajectoryLayout::Group& group = layout.groups[g];
        if (group.hex_cell < 0) { ++layout.n_unresolved_grid_edges; continue; }

        int local[3];
        bool ok = true;
        for (int k = 0; k < 3; ++k)
        {
            if (k == info.firstAxis) local[k] = info.firstCoord - piece.cube[k];
            else if (k == info.secondAxis) local[k] = info.secondCoord - piece.cube[k];
            else local[k] = 0;
            if (local[k] != 0 && local[k] != 1) ok = false;
        }
        if (!ok) { ++layout.n_unresolved_grid_edges; continue; }

        const std::array<int, 3> rootCube = layout.pieces[group.pieces[0]].cube;
        const GridIsomorphism& toCanonical = pieceToCanonical[info.piece];

        int cornerHex[2];
        bool complete = true;
        for (int u = 0; u < 2; ++u)
        {
            int l[3] = { local[0], local[1], local[2] };
            l[info.varying] = u;
            const Vec3d canonical = toCanonical.transform_point(Vec3d(
                    (double)(piece.cube[0] + l[0]),
                    (double)(piece.cube[1] + l[1]),
                    (double)(piece.cube[2] + l[2])));
            const int cx = (int)std::llround(canonical[0] - (double)rootCube[0]);
            const int cy = (int)std::llround(canonical[1] - (double)rootCube[1]);
            const int cz = (int)std::llround(canonical[2] - (double)rootCube[2]);
            if (cx < 0 || cx > 1 || cy < 0 || cy > 1 || cz < 0 || cz > 1)
            {
                complete = false;
                break;
            }
            cornerHex[u] = group.corners[cx + 2 * cy + 4 * cz];
            if (cornerHex[u] < 0) complete = false;
        }
        if (!complete) { ++layout.n_unresolved_grid_edges; continue; }

        int ha = cornerHex[0];
        int hb = cornerHex[1];
        if (ha > hb) std::swap(ha, hb);
        const std::unordered_map<std::pair<int, int>, int, PairHash>::const_iterator it =
                finalEdgeByVertices.find(std::make_pair(ha, hb));
        if (it == finalEdgeByVertices.end())
        {
            ++layout.n_unresolved_grid_edges;
            continue;
        }
        layout.edge_hex_edge[e] = it->second;

        if (layout.vertex_tet_edge[a] >= 0 && layout.vertex_tet_edge[a] == layout.vertex_tet_edge[b])
            layout.edge_tet_edge[e] = layout.vertex_tet_edge[a];
    }

    // -----------------------------------------------------------------------
    // 9. aggregate the correspondence per hex element
    // -----------------------------------------------------------------------
    layout.hex_vertex_vertices.assign(nHexVertices, std::vector<int>());
    layout.hex_edge_edges.assign(nHexEdges, std::vector<int>());
    layout.hex_edge_polyline.assign(nHexEdges, std::vector<int>());
    layout.hex_face_facets.assign(nHexFaces, std::vector<int>());
    layout.hex_cell_pieces.assign(nHexCells, std::vector<int>());
    layout.hex_cell_cells.assign(nHexCells, std::vector<int>());
    layout.hex_vertex_n_refined.assign(nHexVertices, 0);
    layout.hex_edge_n_refined.assign(nHexEdges, 0);
    layout.hex_face_n_refined.assign(nHexFaces, 0);
    layout.hex_cell_n_refined.assign(nHexCells, 0);
    layout.hex_cell_n_tets.assign(nHexCells, 0);

    std::vector<std::vector<int> > edgesOfHexEdge(nHexEdges);
    for (size_t e = 0; e < edgeVertices.size(); ++e)
        if (layout.edge_hex_edge[e] >= 0)
            edgesOfHexEdge[layout.edge_hex_edge[e]].push_back((int)e);

    for (size_t h = 0; h < edgesOfHexEdge.size(); ++h)
    {
        std::vector<int>& edges = edgesOfHexEdge[h];
        if (edges.empty()) continue;

        std::map<int, std::vector<int> > adjacency;
        for (size_t i = 0; i < edges.size(); ++i)
        {
            const int a = edgeVertices[edges[i]][0];
            const int b = edgeVertices[edges[i]][1];
            adjacency[a].push_back(b);
            adjacency[b].push_back(a);
        }

        // start at the endpoint whose hex vertex is the smaller one of the two
        // hex edge endpoints (this gives a deterministic direction)
        const int ha = std::min(view.edges[h][0], view.edges[h][1]);
        int startVertex = -1;
        bool branched = false;
        for (std::map<int, std::vector<int> >::const_iterator it = adjacency.begin(); it != adjacency.end(); ++it)
        {
            if (it->second.size() > 2) branched = true;
            if (it->second.size() == 1 && layout.vertex_hex_vertex[it->first] == ha)
                startVertex = it->first;
        }

        std::vector<int> polyline;
        std::vector<int> orderedEdges;
        if (!branched && startVertex >= 0)
        {
            int previous = -1;
            int current = startVertex;
            while (current >= 0)
            {
                polyline.push_back(current);
                int next = -1;
                for (size_t i = 0; i < adjacency[current].size(); ++i)
                    if (adjacency[current][i] != previous) { next = adjacency[current][i]; break; }
                if (next >= 0)
                {
                    int found = -1;
                    for (size_t i = 0; i < edges.size(); ++i)
                    {
                        const int a = edgeVertices[edges[i]][0];
                        const int b = edgeVertices[edges[i]][1];
                        if ((a == current && b == next) || (a == next && b == current)) { found = edges[i]; break; }
                    }
                    if (found >= 0) orderedEdges.push_back(found);
                }
                previous = current;
                current = next;
            }
        }
        else
        {
            ++layout.n_hex_edges_with_branched_polyline;
            orderedEdges = edges;
            for (size_t i = 0; i < edges.size(); ++i)
            {
                polyline.push_back(edgeVertices[edges[i]][0]);
                polyline.push_back(edgeVertices[edges[i]][1]);
            }
            std::sort(polyline.begin(), polyline.end());
            polyline.erase(std::unique(polyline.begin(), polyline.end()), polyline.end());
        }

        for (size_t i = 0; i < orderedEdges.size(); ++i)
            layout.edge_hex_edge_step[orderedEdges[i]] = (int)i;

        layout.hex_edge_edges[h] = orderedEdges;
        layout.hex_edge_polyline[h] = polyline;
    }

    for (size_t v = 0; v < layout.vertex_position.size(); ++v)
    {
        const int h = layout.vertex_hex_vertex[v];
        if (h < 0 || (size_t)h >= layout.hex_vertex_vertices.size()) continue;
        layout.hex_vertex_vertices[h].push_back((int)v);
        ++layout.hex_vertex_n_refined[h];
    }

    for (size_t e = 0; e < edgeVertices.size(); ++e)
    {
        const int h = layout.edge_hex_edge[e];
        if (h < 0 || (size_t)h >= layout.hex_edge_edges.size()) continue;
        ++layout.hex_edge_n_refined[h];
    }

    for (size_t f = 0; f < layout.facet_vertices.size(); ++f)
    {
        const int h = layout.facet_hex_face[f];
        if (h < 0 || (size_t)h >= layout.hex_face_facets.size()) continue;
        layout.hex_face_facets[h].push_back((int)f);
        ++layout.hex_face_n_refined[h];
    }

    for (size_t g = 0; g < layout.groups.size(); ++g)
    {
        const HexTrajectoryLayout::Group& group = layout.groups[g];
        if (group.hex_cell < 0 || (size_t)group.hex_cell >= layout.hex_cell_pieces.size()) continue;
        for (size_t pi = 0; pi < group.pieces.size(); ++pi)
        {
            const int p = group.pieces[pi];
            layout.hex_cell_pieces[group.hex_cell].push_back(p);
            layout.hex_cell_n_tets[group.hex_cell] += layout.pieces[p].n_cells;
        }
    }

    for (size_t c = 0; c < layout.cell_vertices.size(); ++c)
    {
        const int pieceIndex = layout.cell_piece[c];
        int h = -1;
        if (pieceIndex >= 0 && (size_t)pieceIndex < layout.pieces.size())
            h = layout.pieces[pieceIndex].hex_cell;
        layout.cell_hex_cell[c] = h;
        if (h >= 0 && (size_t)h < layout.hex_cell_cells.size())
        {
            layout.hex_cell_cells[h].push_back((int)c);
            ++layout.hex_cell_n_refined[h];
        }
    }

    // -----------------------------------------------------------------------
    // 10. diagnostics
    // -----------------------------------------------------------------------
    layout.n_refined_vertices = layout.vertex_position.size();
    layout.n_refined_edges = edgeVertices.size();
    layout.n_refined_facets = layout.facet_vertices.size();
    layout.n_refined_cells = layout.cell_vertices.size();
    layout.n_pieces = layout.pieces.size();

    for (size_t h = 0; h < layout.hex_vertex_n_refined.size(); ++h)
        if (layout.hex_vertex_n_refined[h] == 0) ++layout.n_hex_vertices_without_trajectory;
    for (size_t h = 0; h < layout.hex_edge_n_refined.size(); ++h)
        if (layout.hex_edge_n_refined[h] == 0) ++layout.n_hex_edges_without_trajectory;
    for (size_t h = 0; h < layout.hex_face_n_refined.size(); ++h)
        if (layout.hex_face_n_refined[h] == 0) ++layout.n_hex_faces_without_trajectory;
    for (size_t h = 0; h < layout.hex_cell_n_refined.size(); ++h)
        if (layout.hex_cell_n_refined[h] == 0) ++layout.n_hex_cells_without_trajectory;

    for (size_t e = 0; e < edgeVertices.size(); ++e)
    {
        if (edgeFacetB[e] >= 0) continue;
        ++layout.n_open_edges;
        const int f = edgeFacetA[e];
        bool onDomainBoundary = false;
        if (f >= 0)
        {
            const int tf = layout.facet_tet_face[f];
            onDomainBoundary = (tf >= 0) && inputMesh.is_boundary(FaceHandle(tf));
        }
        if (!onDomainBoundary) ++layout.n_crack_edges;
    }
}

} // namespace HexEx
