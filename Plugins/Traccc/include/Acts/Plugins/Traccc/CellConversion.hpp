// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

// Plugin include(s)
#include "Acts/Plugins/Traccc/BarcodeMap.hpp"
#include "Acts/Plugins/Traccc/Detail/Module.hpp"

// Detray include(s)
#include "detray/core/detector.hpp"

// Traccc include(s)
#include "traccc/edm/cell.hpp"
#include "traccc/geometry/geometry.hpp"
#include "traccc/io/digitization_config.hpp"
#include "traccc/io/reader_edm.hpp"
#include "traccc/io/read_geometry.hpp"

// System include(s)
#include <cstdint>
#include <cstdlib>
#include <map>
#include <tuple>
#include <memory>
#include <utility>

// temp
#include "traccc/io/read_cells.hpp"
#include "traccc/io/data_format.hpp"

// This code is borrowed from traccc/io/src/csv/read_cells.cpp with minor modifications
// so that a cells map rather than a file path is needed.
namespace {

/// Comparator used for sorting cells. This sorting is one of the assumptions
/// made in the clusterization algorithm
struct cell_order {
    bool operator()(const traccc::cell& lhs, const traccc::cell& rhs) const {
        if (lhs.module_link != rhs.module_link) {
            return lhs.module_link < rhs.module_link;
        } else if (lhs.channel1 != rhs.channel1) {
            return (lhs.channel1 < rhs.channel1);
        } else {
            return (lhs.channel0 < rhs.channel0);
        }
    }
};  // struct cell_order

}

namespace Acts::TracccPlugin{

/// @brief Converts a geometry ID -> traccc cells map to traccc cells and modules.
/// @param mr The memory resource to use.
/// @param cellsMap A map from Acts geometry ID (value) to traccc cell.
/// @param geom The traccc geometry 
/// @param dconfig The traccc digitization configuration.
/// @param barcode_map A map from Acts geometry IDs (value) to detray barcodes.
/// @return A tuple containing the traccc cells (first item) and traccc modules (second item).
inline auto createCellsAndModules(
    vecmem::memory_resource* mr,
    std::map<std::uint64_t, std::vector<traccc::cell>> cellsMap,
    const traccc::geometry* geom,
    const traccc::digitization_config* dconfig,
    const std::map<std::uint64_t, detray::geometry::barcode>* barcode_map) {

    traccc::io::cell_reader_output out(mr);

    // Sort the cells. Deduplication or not, they do need to be sorted.
    for (auto& [_, cells] : cellsMap) {
        std::sort(cells.begin(), cells.end(), cell_order());
    }

    // Fill the output containers with the ordered cells and modules.
    for (const auto& [original_geometry_id, cells] : cellsMap) {
        // Modify the geometry ID of the module if a barcode map is
        // provided.
        std::uint64_t geometry_id = original_geometry_id;
        if (barcode_map != nullptr) {
            const auto it = barcode_map->find(geometry_id);
            if (it != barcode_map->end()) {
                geometry_id = it->second.value();
            } else {
                throw std::runtime_error(
                    "Could not find barcode for geometry ID " +
                    std::to_string(geometry_id));
            }
        }

        // Add the module and its cells to the output.
        out.modules.push_back(
            Detail::get_module(geometry_id, geom, dconfig, original_geometry_id));
        for (auto& cell : cells) {
            out.cells.push_back(cell);
            // Set the module link.
            out.cells.back().module_link = out.modules.size() - 1;
        }
    }
    return std::make_tuple(std::move(out.cells), std::move(out.modules));
}

}
