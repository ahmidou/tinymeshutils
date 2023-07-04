/*
MIT License

Copyright (c) 2019 Syoyo Fujita

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#ifndef TINYMESHUTILS_HALF_EDGE_HH_
#define TINYMESHUTILS_HALF_EDGE_HH_

// TODO(syoyo): Support uint64_t for Edge.

#include <algorithm>
#include <cstdint>
#include <vector>
#include "tsl/robin_map.h"
//#include <iostream>  // DBG

// Using robin-map may give better performance

#if defined(TINYMESHUTILS_USE_ROBINMAP)
#if __has_include("tsl/robin_map.h")
//#error TODO
#include <tsl/robin_map.h>
template<class K, class V, class _Hasher = std::hash<K>>
using umap = tsl::robin_map<K, V>;
#endif
#else
#include <unordered_map>
template<class K, class V, class _Hasher = std::hash<K>>
using umap = std::unordered_map<K, V, _Hasher>;
#endif

namespace tinymeshutils {

struct Edge {
  Edge(const uint32_t _v0, const uint32_t _v1) : v0(_v0), v1(_v1) {}

  // create an 64 bit identifier that is unique for the not oriented edge
  operator uint64_t() const {
    uint32_t p0 = v0, p1 = v1;
    if (p0 < p1) std::swap(p0, p1);
    return (uint64_t(p0) << 32) | uint64_t(p1);
  }

  uint32_t v0, v1;  // start, end
};

struct HalfEdge {
  // invalid or undefined = MAXUINT32
  uint32_t twin{MAXUINT32};  // index to the halfedges array.
  uint32_t next{MAXUINT32};      // index to the halfedges array.
  uint32_t vertex_index{MAXUINT32}; // vertex index at the start of the edge
  uint32_t face_index{MAXUINT32};  // index to face indices
  uint32_t index{MAXUINT32};  // index to edge indices
};

// Functor object for Edge
struct EdgeHash {
  size_t operator()(const std::pair<uint32_t, uint32_t> &k) const {
    if (sizeof(size_t) == 8) {
      // We can create unique value
      return (uint64_t(k.first) << 32) | uint64_t(k.second);
    } else {
      // 32bit environment. we cannot create unique hash value.
      return std::hash<uint32_t>()(k.first) ^ std::hash<uint32_t>()(k.second);
    }
  }
};

///
/// Simple half-edge construction for polygons.
///
/// @param[in] face_vert_indices Array of face vertex indices.
/// @param[in] face_vert_counts Array of the number of vertices for each face(3
/// = triangle, 4 = quad, ...).
/// @param[out] edges Array of edges constructed
/// @param[out] halfedges Array of half-edges comstructed. length =
/// face_vert_indices.size()
/// @param[out] vertex_starting_halfedge_indices Index to starting half-edge in
/// `halfedges` for each vertex. MAXUINT32 when no half-edge assigned to its vertex.
/// length = the highest value in `face_vert_indices`(= the number of vertices)
/// face_vert_indices.size()
/// @param[out] err Optional error message.
///
/// @return true upon success. Return false when input mesh has invalid
/// topology.
///
/// face_vert_indices.size() is equal to the sum of each elements in
/// `face_vert_counts`. Assume edge is defined in v0 -> v1, v1 -> v2, ... v(-1)
/// -> v0 order.
///
bool BuildHalfEdge(const std::vector<uint32_t> &face_vert_indices,
                   const std::vector<uint32_t> &face_vert_counts,
                   std::vector<Edge> *edges, std::vector<HalfEdge> *halfedges,
                   std::vector<uint32_t> *vertex_starting_halfedge_indices,std::vector<uint32_t> *faceData,
                   std::string *err = nullptr);

}  // namespace tinymeshutils

#endif  // TINYMESHUTILS_HALF_EDGE_HH_

#if defined(TINYMESHUTILS_HALF_EDGE_IMPLEMENTATION)
namespace tinymeshutils {

bool BuildHalfEdge(const std::vector<uint32_t> &face_vert_indices,
                   const std::vector<uint32_t> &face_vert_counts,
                   std::vector<Edge> *edges, std::vector<HalfEdge> *halfedges,
                   std::vector<uint32_t> *vertex_starting_halfedge_indices, std::vector<uint32_t> *faceData,
                   std::string *err) {
  // Based on documents at Internet and an implementation
  // https://github.com/yig/halfedge

  std::vector<uint32_t> faceData_buf(face_vert_counts);

  size_t num_indices = 0;
  for (size_t i = 0; i < face_vert_counts.size(); i++) {
    if (face_vert_counts[i] < 3) {
      // invalid # of vertices for a face.
      return false;
    }
    num_indices += face_vert_counts[i];
  }

  // Find larget number = the number of vertices in input mesh.
  uint32_t num_vertices = (*std::max_element(face_vert_indices.begin(), face_vert_indices.end())) +1;
  std::vector<uint32_t> vertex_starting_halfedge_indices_buf(num_vertices, MAXUINT32);

  // allocate buffer for half-edge.
  std::vector<HalfEdge> halfedge_buf(num_indices);

  umap<std::pair<int32_t, int32_t>, size_t, EdgeHash>
      halfedge_table;  // <<v0, v1>, index to `half_edges`>

  //
  // 1. Build edges.
  //
  std::vector<Edge> edge_buf;  // linear array of edges.
  umap<uint64_t, size_t>
      edge_map;  // <un oriented edge index, index to edge array>

  {
    std::vector<Edge>
        tmp_edges;  // linear array of edges. may have duplicated edges

    // list up and register edges.
    size_t f_offset = 0;
    for (size_t f = 0; f < face_vert_counts.size(); f++) {
      uint32_t nv = face_vert_counts[f];
      if (nv < 3) {
        if (err) {
          (*err) = "Face " + std::to_string(f) + " has invalid # of vertices " +
                   std::to_string(nv) + "\n";
        }
        return false;
      }

      uint32_t ne = nv;

      // for each edge
      for (size_t e = 0; e < size_t(ne); e++) {
        // std::cout << "e = " << e << ", " << (e + 1) % nv << "\n";
        uint32_t v0 = face_vert_indices[f_offset + e];
        uint32_t v1 = face_vert_indices[f_offset + (e + 1) % nv];
        // std::cout << "v0 = " << v0 << ", v1 = " << v1 << "\n";

        tmp_edges.push_back({v0, v1});
      }

      f_offset += nv;
    }

    // create edge_map and unique array of edges.
    for (const auto &edge : tmp_edges) {
      uint64_t key = edge;
      if (!edge_map.count(key)) {
        size_t edge_idx = edge_buf.size();
        edge_map[edge] = edge_idx;

        edge_buf.push_back(edge);
      }
    }
  }

#if 0
  // dbg
  for (size_t i = 0; i < edge_buf.size(); i++) {
    std::cout << "edge[" << i << "] = " << edge_buf[i].v0 << ", "
              << edge_buf[i].v1 << "\n";
  }
#endif

  //
  // 2. Register half edges
  //
  {
    size_t f_offset = 0;
    for (size_t f = 0; f < face_vert_counts.size(); f++) {
      // for each edge
      uint32_t nv = face_vert_counts[f];
      if (nv < 3) {
        if (err) {
          (*err) = "Face " + std::to_string(f) + " has invalid # of vertices " +
                   std::to_string(nv) + "\n";
        }
        return false;
      }

      uint32_t ne = nv;

      // register
      for (size_t e = 0; e < size_t(ne); e++) {
        uint32_t v0 = face_vert_indices[f_offset + e];
        uint32_t v1 = face_vert_indices[f_offset + (e + 1) % nv];

        Edge edge(v0, v1);

        // vertex pair must be unique over the input mesh
        if (halfedge_table.count(std::make_pair(v0, v1))) {
          if (err) {
            (*err) = "(Register half edges). Invalid topology. Edge (v0: " +
                     std::to_string(v0) + ", v1: " + std::to_string(v1) +
                     ") must be unique but duplicated one exists for Face " +
                     std::to_string(f) + "\n";
          }

          return false;
        }

        uint64_t eid = edge;  // non oriented

        if (!edge_map.count(eid)) {
          if (err) {
            (*err) = "??? Edge (v0: " + std::to_string(edge.v0) +
                     ", v1: " + std::to_string(edge.v1) +
                     ") does not found for Face " + std::to_string(f) +
                     ". This should not be happen.\n";
          }
          return false;
        }

        size_t index = edge_map[eid];

        size_t halfedge_offset = f_offset + e;

        halfedge_table[std::make_pair(v0, v1)] = halfedge_offset;

        halfedge_buf[halfedge_offset].index = uint32_t(index);
        halfedge_buf[halfedge_offset].vertex_index = uint32_t(v0);
        halfedge_buf[halfedge_offset].face_index = uint32_t(f);
        halfedge_buf[halfedge_offset].next = uint32_t(f_offset + (e + 1) % nv);

        if(e==0)
          faceData_buf[f] = uint32_t(halfedge_offset);

        if (size_t(v0) >= vertex_starting_halfedge_indices_buf.size()) {
          if (err) {
            (*err) =
                "Out-of-bounds access. v0 " + std::to_string(v0) +
                " must be less than " +
                std::to_string(vertex_starting_halfedge_indices_buf.size()) +
                "\n";
          }
          return false;
        }

        if (vertex_starting_halfedge_indices_buf[size_t(v0)] == MAXUINT32) {
          // Set as starting half-edge
          vertex_starting_halfedge_indices_buf[size_t(v0)] = uint32_t(halfedge_offset);
        }
      }

      f_offset += nv;
    }
  }

  // dbg
  // for (size_t i = 0; i < halfedge_buf.size(); i++) {
  //  std::cout << "halfedge_buf[" << i << "].index = " <<
  //  halfedge_buf[i].index << "\n";
  //}

  //
  // 3. Find opposite half edges
  //
  for (size_t i = 0; i < halfedge_buf.size(); i++) {
    HalfEdge &halfedge = halfedge_buf[i];
    if ((halfedge.index == MAXUINT32) ||
        (halfedge.index >= uint32_t(edge_buf.size()))) {
      if (err) {
        (*err) = "Invalid index " + std::to_string(halfedge.index) +
                 ". Must be >= 0 and < " + std::to_string(edge_buf.size()) +
                 "\n";
      }

      return false;
    }

    const Edge &edge = edge_buf[size_t(halfedge.index)];

    if (halfedge_table.count(std::make_pair(edge.v1, edge.v0)) &&
        halfedge_table.count(std::make_pair(edge.v0, edge.v1))) {
      // Opposite halfedge exists. Make a link.

      size_t halfedge_index0 =
          halfedge_table.at(std::make_pair(edge.v0, edge.v1));

      size_t halfedge_index1 =
          halfedge_table.at(std::make_pair(edge.v1, edge.v0));

      if (halfedge_index0 == halfedge_index1) {
        if (err) {
          (*err) = "Invalid halfedge_index. Both indices has same value.\n";
        }
        return false;
      }

      // Check if self-referencing. Choose different index compared to current
      // index(`i`).
      size_t opposite_halfedge_index =
          (halfedge_index0 == i) ? halfedge_index1 : halfedge_index0;

      if (opposite_halfedge_index >= halfedge_buf.size()) {
        if (err) {
          (*err) = "Invalid halfedge_index " +
                   std::to_string(opposite_halfedge_index) + ". Must be < " +
                   std::to_string(halfedge_buf.size()) + "\n";
        }
        return false;
      }

      HalfEdge &twin = halfedge_buf[opposite_halfedge_index];

      if (twin.index != halfedge.index) {
        if (err) {
          (*err) = "Edge id mismatch. twin.index " +
                   std::to_string(twin.index) +
                   " must be equal to halfedge.index " +
                   std::to_string(halfedge.index) + "\n";
        }
        return false;
      }

      // Make a link.
      halfedge.twin = uint32_t(opposite_halfedge_index);
    }
  }

  (*faceData) = faceData_buf;
  (*edges) = edge_buf;
  (*halfedges) = halfedge_buf;
  (*vertex_starting_halfedge_indices) = vertex_starting_halfedge_indices_buf;

#if 0
  // dbg
  std::cout << "halfedge_buf.size = " << halfedge_buf.size() << "\n";

  for (size_t i = 0; i < halfedge_buf.size(); i++) {
    std::cout << "halfedge[" << i << "]. face = " << halfedge_buf[i].face_index
              << ", edge = " << halfedge_buf[i].index
              << ", opposite he = " << halfedge_buf[i].twin
              << ", next he = " << halfedge_buf[i].next << "\n";
  }

  for (size_t i = 0; i < vertex_starting_halfedge_indices_buf.size(); i++) {
    std::cout << "v[" << i << "].halfedge_index = " << vertex_starting_halfedge_indices_buf[i] << "\n";
  }
#endif

  return true;
}

}  // namespace tinymeshutils
#endif  // TINYMESHUTILS_HALF_EDGE_IMPLEMENTATION
