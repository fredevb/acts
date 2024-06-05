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

// This code is borrowed from traccc/io/src/csv/read_cells.cpp with minor modifications
// so that a cells map rather than a file path is needed.
namespace Acts::TracccPlugin::Detail {

/// Helper function which finds module from csv::cell in the geometry and
/// digitization config, and initializes the modules limits with the cell's
/// properties
inline traccc::cell_module get_module(const std::uint64_t geometry_id,
                               const traccc::geometry* geom,
                               const traccc::digitization_config* dconfig,
                               const std::uint64_t original_geometry_id) {

    traccc::cell_module result;
    result.surface_link = detray::geometry::barcode{geometry_id};

    // Find/set the 3D position of the detector module.
    if (geom != nullptr) {

        // Check if the module ID is known.
        if (!geom->contains(result.surface_link.value())) {
            throw std::runtime_error(
                "Could not find placement for geometry ID " +
                std::to_string(result.surface_link.value()));
        }

        // Set the value on the module description.
        result.placement = (*geom)[result.surface_link.value()];
    }

    // Find/set the digitization configuration of the detector module.
    if (dconfig != nullptr) {

        // Check if the module ID is known.
        const traccc::digitization_config::Iterator geo_it =
            dconfig->find(original_geometry_id);
        if (geo_it == dconfig->end()) {
            throw std::runtime_error(
                "Could not find digitization config for geometry ID " +
                std::to_string(original_geometry_id));
        }

        // Set the value on the module description.
        const auto& binning_data = geo_it->segmentation.binningData();
        assert(binning_data.size() >= 2);
        result.pixel = {binning_data[0].min, binning_data[1].min,
                        binning_data[0].step, binning_data[1].step};
    }

    return result;
}

}
