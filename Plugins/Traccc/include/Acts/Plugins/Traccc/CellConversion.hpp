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

// temp
#include "traccc/io/read_cells.hpp"
#include "traccc/io/data_format.hpp"

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

namespace Acts::TracccPlugin{

template <typename Cell>
struct CellDataGetter{
    virtual double getTime(const Cell& cell) const = 0;
    virtual double getActivation(const Cell& cell) const = 0;
    virtual int getRow(const Cell& cell) const = 0;
    virtual int getColumn(const Cell& cell) const = 0;
};

template <typename CellCollection>
std::map<std::uint64_t, std::vector<traccc::cell>> newTracccCellsMap(const std::map<Acts::GeometryIdentifier, CellCollection>& map, const CellDataGetter<typename CellCollection::value_type>& getter){
    std::map<std::uint64_t, std::vector<traccc::cell>> tracccCellMap;
    for (const auto& [geometryID, cells] : map){
        std::vector<traccc::cell> tracccCells;
        for (const auto& cell : cells){
            traccc::cell tracccCell{
                static_cast<traccc::channel_id>(getter.getRow(cell)),
                static_cast<traccc::channel_id>(getter.getColumn(cell)),        
                static_cast<traccc::scalar>(getter.getActivation(cell)),
                static_cast<traccc::scalar>(getter.getTime(cell)),
                0//tracccModules.size() - 1
            };
            tracccCells.push_back(std::move(tracccCell));
        }
        tracccCellMap.insert({geometryID.value(), std::move(tracccCells)});
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

template <typename metadata_t, typename container_t, typename cell_data_getter_t>
class CellDataConverter{
    public:
    using detector_t = detray::detector<metadata_t, container_t>;
    CellDataConverter(const detector_t& det, const Acts::GeometryHierarchyMap<Acts::BinUtility>& segs, const cell_data_getter_t& g) :
    detector(det),
    digitizationConfig(getConfig(segs)),
    getter(g),
    surfaceTransforms(getSurfaceTransforms(det)),
    barcodeMap(getBarcodeMap(det)){}

    template <typename host_mr_t, typename CellCollection>
    auto operator()(host_mr_t& mr, std::map<Acts::GeometryIdentifier, CellCollection> map) const{
        std::cout << "Running Conversion" << std::endl;
        const auto cellsMap = newTracccCellsMap(map, getter);
        std::cout << "Created Cells Map" << std::endl;
        traccc::io::cell_reader_output readOut(&mr);
        Detail::read_cells(readOut, cellsMap, &surfaceTransforms, &digitizationConfig, &barcodeMap);
        //                traccc::io::read_cells(readOut, "/home/frederik/Desktop/CERN-TECH-OTHER/traccc/data/tml_pixels/event000000000-cells.csv",
        //                               traccc::data_format::csv, &surfaceTransforms,
        //                               &digitizationConfig, &barcodeMap);
        std::cout << "Read Cells" << std::endl;
        return std::make_tuple(std::move(readOut.cells), std::move(readOut.modules));
    }
    
    private:
    const detector_t& detector;
    const traccc::digitization_config digitizationConfig;
    const cell_data_getter_t getter;
    const traccc::geometry surfaceTransforms;
    const std::map<std::uint64_t, detray::geometry::barcode> barcodeMap;
};

}