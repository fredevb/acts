// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Geometry/GeometryIdentifier.hpp"

// Detray include(s)
#include "detray/core/detector.hpp"

// System include(s)
#include <cstdint>
#include <cstdlib>
#include <map>
#include <tuple>
#include <memory>
#include <utility>

namespace Acts::TracccPlugin {

/// @brief Creates an map from Acts geometry ID (value) to detray barcode
/// by only using the information contained in the detray detector.
/// @param detector the detray detector.
/// @return A map: geometry ID value type -> detray geometry barcode.
template <typename metadata_t, typename container_t>
inline std::map<std::uint64_t, detray::geometry::barcode> createBarcodeMap(const detray::detector<metadata_t, container_t>& detector){
    // Construct a map from Acts surface identifiers to Detray barcodes.
    std::map<std::uint64_t, detray::geometry::barcode> barcode_map;
    for (const auto& surface : detector.surfaces()) {
        barcode_map[surface.source] = surface.barcode();
    }
    return barcode_map;
}

std::uint64_t barcodeMapTry(const std::map<std::uint64_t, detray::geometry::barcode>* barcode_map, std::uint64_t geometry_id){
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
    return geometry_id;
}
}
