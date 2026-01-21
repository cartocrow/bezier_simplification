#pragma once

#include <concepts>
#include "graph_curve_traits_2.h"

namespace cartocrow {
template<class G>
concept Graph2Like =
requires(G g) {
    typename G::Vertex;
    typename G::Vertex_data;
    typename G::Vertex_handle;
    typename G::Vertex_const_handle;
    typename G::Vertex_iterator;
    typename G::Vertex_const_iterator;
    typename G::Edge;
    typename G::Edge_data;
    typename G::Edge_handle;
    typename G::Edge_const_handle;
    typename G::Edge_iterator;
    typename G::Edge_const_iterator;
    typename G::Kernel;
    typename G::Point_2;
    typename G::Curve_2;

    requires GraphCurveTraits_2<typename G::Curve_traits>;

    /// Retrieves the vertices begin iterator
    {
    g.vertices_begin()
    } -> std::same_as<typename G::Vertex_iterator>;

    /// Retrieves the vertices end iterator
    {
    g.vertices_end()
    } -> std::same_as<typename G::Vertex_iterator>;

    /// Retrieves the edges begin iterator
    {
    g.edges_begin()
    } -> std::same_as<typename G::Edge_iterator>;

    /// Retrieves the edges end iterator
    {
    g.edges_end()
    } -> std::same_as<typename G::Edge_iterator>;

    /// Retrieves the number of edges
    {
    g.number_of_edges()
    } -> std::same_as<std::size_t>;

    /// Retrieves the number of vertices
    {
    g.number_of_vertices()
    } -> std::same_as<std::size_t>;
};
}