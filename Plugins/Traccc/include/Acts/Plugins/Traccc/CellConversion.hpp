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

// VecMem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

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

/// Helper function which finds module from csv::cell in the geometry and
/// digitization config, and initializes the modules limits with the cell's
/// properties
traccc::cell_module get_module(const std::uint64_t geometry_id,
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

void read_cells(
    traccc::io::cell_reader_output& out,
    std::map<std::uint64_t, std::vector<traccc::cell>> cellsMap,
    const traccc::geometry* geom,
    const traccc::digitization_config* dconfig,
    const std::map<std::uint64_t, detray::geometry::barcode>* barcode_map) {

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
            get_module(geometry_id, geom, dconfig, original_geometry_id));
        for (auto& cell : cells) {
            out.cells.push_back(cell);
            // Set the module link.
            out.cells.back().module_link = out.modules.size() - 1;
        }
    }
}

}

namespace Acts::TracccPlugin::Internal {
    template <typename, typename T = void>
    struct cellTypeHasRequiredFunctions : std::false_type {};

    template <typename T>
    struct cellTypeHasRequiredFunctions<
        T,
        std::void_t<decltype(getCellRow(std::declval<T>())),
                    decltype(getCellColumn(std::declval<T>())),
                    decltype(getCellActivation(std::declval<T>())),
                    decltype(getCellTime(std::declval<T>()))>> : std::true_type {
    };

    template <typename T>
    constexpr void staticCheckCellType() {
    constexpr bool hasFns = cellTypeHasRequiredFunctions<T>();
    static_assert(hasFns,
                    "Cell type should have the following functions: "
                    "'int getCellRow(const Cell&)', "
                    "'int getCellColumn(const Cell&)', "
                    "'int getCellActivation(const Cell&)', "
                    "'int getCellTime(const Cell&)'");
    }
}

namespace Acts::TracccPlugin::CellConversion{

template <typename CellCollection>
std::map<std::uint64_t, std::vector<traccc::cell>> newCellsMap(const std::map<Acts::GeometryIdentifier, CellCollection>& map){
    using Cell = typename CellCollection::value_type;
    Internal::staticCheckCellType<Cell>();

    std::map<std::uint64_t, std::vector<traccc::cell>> tracccCellMap;
    for (const auto& [geometryID, cells] : map){
        std::vector<traccc::cell> tracccCells;
        for (const auto& cell : cells){
            traccc::cell tracccCell{
                getCellRow(cell),
                getCellColumn(cell),        
                getCellActivation(cell), //cell.activation
                getCellTime(cell)
                //tracccModules.size() - 1
            };
            tracccCells.push_back(tracccCell);
        }
        tracccCellMap.insert({geometryID.value(), tracccCells});
    }
    return tracccCellMap;
}

template <typename metadata_t, typename container_t>
traccc::geometry getSurfaceTransforms(const detray::detector<metadata_t, container_t>& detector) {
    return traccc::io::alt_read_geometry(detector);
}

template <typename metadata_t, typename container_t>
std::map<std::uint64_t, detray::geometry::barcode> getBarcodeMap(const detray::detector<metadata_t, container_t>& detector){
    // Construct a map from Acts surface identifiers to Detray barcodes.
    std::map<std::uint64_t, detray::geometry::barcode> barcode_map;
    for (const auto& surface : detector.surfaces()) {
        barcode_map[surface.source] = surface.barcode();
    }
    return barcode_map;
}

traccc::digitization_config getConfig(const Acts::GeometryHierarchyMap<Acts::BinUtility>& segmentations){
    using elem_t = std::pair<Acts::GeometryIdentifier, traccc::module_digitization_config>;
    std::vector<elem_t> vec;
    for (auto& e : segmentations.getElements()){
        vec.push_back({e.first, traccc::module_digitization_config{e.second}});
    }
    return traccc::digitization_config(vec);
}

template <typename metadata_t, typename container_t>
class CellDataConverter{
    public:
    using detector_type = detray::detector<metadata_t, container_t>;
    CellDataConverter(vecmem::memory_resource& mr, const detector_type& det, Acts::GeometryHierarchyMap<Acts::BinUtility> segs) :
    memoryResource(mr),
    detector(det),
    digitizationConfig(getConfig(segs)),
    surface_transforms(getSurfaceTransforms(det)),
    barcodeMap(getBarcodeMap(det)){}

    template <typename CellsCollection>
    auto operator()(std::map<Acts::GeometryIdentifier, CellsCollection> map) const{
        const auto cellsMap = newCellsMap(map);
        traccc::io::cell_reader_output readOut(&memoryResource);
        Detail::read_cells(readOut, cellsMap, &surface_transforms, &digitizationConfig, &barcodeMap);
        return std::make_tuple(std::move(readOut.cells), std::move(readOut.modules));
    }
    
    private:

    vecmem::memory_resource& memoryResource;
    const detector_type& detector;
    const traccc::digitization_config digitizationConfig;
    const traccc::geometry surface_transforms;
    const std::map<std::uint64_t, detray::geometry::barcode> barcodeMap;
};

}