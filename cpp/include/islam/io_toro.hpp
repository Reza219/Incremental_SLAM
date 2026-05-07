#pragma once

#include "islam/graph.hpp"

#include <filesystem>

namespace islam {

Graph read_graph_toro(const std::filesystem::path& path);

} // namespace islam
