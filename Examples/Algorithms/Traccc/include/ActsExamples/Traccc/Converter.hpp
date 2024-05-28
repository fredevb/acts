// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Plugin include(s)
#include "Acts/Plugins/Traccc/CellConversion.hpp"
#include "Acts/Plugins/Traccc/TrackConversion.hpp"
#include "Acts/Plugins/Traccc/MeasurementConversion.hpp"
#include "Acts/Plugins/Covfie/CovfieConversion.hpp"

// Acts include(s)
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "ActsExamples/EventData/Track.hpp"

// Acts examples include(s)
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"

// Traccc include(s)
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation_algorithm.hpp"
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/edm/cell.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/io/digitization_config.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/propagator/rk_stepper.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <map>

namespace ActsExamples::TracccConversion {

auto readDetector(vecmem::memory_resource* mr, const std::string& detectorFilePath, const std::string& materialFilePath = "", const std::string& gridFilePath = "")
{
    using detector_type = detray::detector<detray::default_metadata,
                                        detray::host_container_types>;
    // Set up the detector reader configuration.
    detray::io::detector_reader_config cfg;
    cfg.add_file(detectorFilePath);
    if (!materialFilePath.empty()) {
        cfg.add_file(materialFilePath);
    }
    if (!gridFilePath.empty()) {
        cfg.add_file(gridFilePath);
    }

    // Read the detector.
    auto [det, names] = detray::io::read_detector<detector_type>(*mr, cfg);
    return std::move(det);
}

/// @brief Gets the time of the cell.
/// @note Currently, it always returns 0.
float getTime(const Cluster::Cell& /*cell*/){
    return 0.f;
}

/// @brief Gets the activation of the cell.
float getActivation(const Cluster::Cell& cell){
    return static_cast<float>(cell.activation);
}

/// @brief Gets the row of the cell.
unsigned int getRow(const Cluster::Cell& cell){
    if (cell.bin[0] > UINT_MAX) {
        throw std::runtime_error("Overflow will occur when casting to unsigned int.");
    }
    return static_cast<unsigned int>(cell.bin[0]);
}

/// @brief Gets the column of the cell.
unsigned int getColumn(const Cluster::Cell& cell){
        if (cell.bin[0] > UINT_MAX) {
        throw std::runtime_error("Overflow will occur when casting to unsigned int.");
    }
    return static_cast<unsigned int>(cell.bin[1]);
}

Acts::BinUtility getSegmentation(const DigiComponentsConfig& dcc){
    return dcc.geometricDigiConfig.segmentation;
}

/// @brief Converts a "geometry ID -> generic cell collection type" map to a "geometry ID -> traccc cell collection" map.
/// @note The function sets the module link of the cells in the output to 0.
/// @return Map from geometry ID to its cell data (as a vector of traccc cell data)
template <typename CellCollection>
std::map<std::uint64_t, std::vector<traccc::cell>> tracccCellsMap(
    const std::map<Acts::GeometryIdentifier, CellCollection>& map)
    {
    std::map<std::uint64_t, std::vector<traccc::cell>> tracccCellMap;
    for (const auto& [geometryID, cells] : map){
        std::vector<traccc::cell> tracccCells;
        for (const auto& cell : cells){
            tracccCells.push_back(
                traccc::cell{
                    getRow(cell),
                    getColumn(cell),        
                    getActivation(cell),
                    getTime(cell),
                    0
                }
            );
        }
        tracccCellMap.insert({geometryID.value(), std::move(tracccCells)});
    }
    return tracccCellMap;
}

/// @brief Creates a traccc digitalization config from an Acts geometry hierarchy map
/// that contains the digitization configuration.
/// @param config the Acts geometry hierarchy map that contains the digitization configuration.
/// @param getSegmentation a function that gets the segmentation from an item in the geometry hierarchy map.
/// @return a traccc digitization config.
template <typename data_t, typename get_segmentation_fn_t>
traccc::digitization_config tracccConfig(
    const Acts::GeometryHierarchyMap<data_t>& config,
    const get_segmentation_fn_t& getSegmentation){
    using elem_t = std::pair<Acts::GeometryIdentifier, traccc::module_digitization_config>;
    std::vector<elem_t> vec;
    for (auto& e : config.getElements()){
        vec.push_back({e.first, traccc::module_digitization_config{getSegmentation(e.second)}});
    }
    return traccc::digitization_config(vec);
}

/// @brief Creates a source links and measurements. 
/// Creates the measurements by copying the data in the traccc measurements.
/// @param detector The detray detector
/// @param measurements The traccc measurements
/// @return A tuple containing a vector of source links and a vector of Acts measurements.
/// @note The indices of the measurements and source links corresponds to the index of their
/// respective traccc measurement in the traccc measurement vector.
template <typename detector_t, typename allocator_t>
auto getIndexSourceLinkAndMeasurements(const detector_t& detector, const std::vector<traccc::measurement, allocator_t>& measurements){
    std::vector<Acts::SourceLink>sourceLinks;
    std::vector<Acts::Measurement<Acts::BoundIndices, 2>> measurementContainer;

    for (const traccc::measurement& m : measurements) 
    {
        Acts::GeometryIdentifier moduleGeoId(detector.surface(m.surface_link).source);
        Index measurementIdx = measurementContainer.size();
        IndexSourceLink idxSourceLink{moduleGeoId, measurementIdx};
        sourceLinks.insert(sourceLinks.end(), Acts::SourceLink{idxSourceLink});
        measurementContainer.push_back(Acts::TracccPlugin::measurement<2UL>(m, Acts::SourceLink{idxSourceLink}));
    }
    return std::make_tuple(std::move(sourceLinks), std::move(measurementContainer));
}

/// @brief This class provides functions for converting input and output data using the data structures in Acts Examples.
/// @tparam detector_t The type of the detray detector this converter expects.
template <typename detector_t>
class Converter{

    public:
    
    Converter(
        std::shared_ptr<const Acts::TrackingGeometry> tg,
        std::shared_ptr<const detector_t> det,
        const Acts::GeometryHierarchyMap<DigiComponentsConfig>& dccMap):
        trackingGeometry(tg),
        detector(det),
        digitizationConfig(tracccConfig(dccMap, getSegmentation)),
        surfaceTransforms(traccc::io::alt_read_geometry(*det)),
        barcodeMap(Acts::TracccPlugin::createBarcodeMap(*det))
    {}

    /// @brief Converts a map of cells to traccc input
    /// @param map The geometry id -> cell collection map which corresponds to each geometry's cells.
    /// @param mr The memory resource to use.
    /// @returns The converted input as a tuple containing traccc input data, i.e. (cells, modules).
    auto convertInput(const std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>& map, vecmem::memory_resource* mr) const {
        auto tcm = tracccCellsMap(map);
        return Acts::TracccPlugin::createCellsAndModules(
            mr, tcm, &surfaceTransforms, &digitizationConfig, &barcodeMap);
    }

    /// @brief Converts a container of traccc tracks to a container of Acts tracks.
    /// @param ms The traccc measurements 
    /// (this is needed to set the source links and calibrated data of the newly converted track states).
    /// @param ts The (indexable) container of traccc tracks.
    /// @return An Acts const track container containing the same track data as the traccc tracks container.
    /// Furthermore, the track states sourcelinks and measurements are also set.
    template <typename traccc_track_container_t, typename allocator_t>
    auto convertOutput(const std::vector<traccc::measurement, allocator_t>& ms, const traccc_track_container_t& ts) const {
        auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
        auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
        TrackContainer tracks(trackContainer, trackStateContainer);

        auto [sourceLinks, measurements] = getIndexSourceLinkAndMeasurements(*detector, ms);

        Acts::TracccPlugin::makeTracks(ts, tracks, *detector, *trackingGeometry, measurements, sourceLinks);

        ConstTrackContainer constTracks{
            std::make_shared<Acts::ConstVectorTrackContainer>(std::move(*trackContainer)),
            std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*trackStateContainer))
            };

        return constTracks;
    }
    
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const detector_t> detector;

    // Cache the converted digitalization configuration, the surface transforms, and the barcode map.
    const traccc::digitization_config digitizationConfig;
    const traccc::geometry surfaceTransforms;
    const std::map<std::uint64_t, detray::geometry::barcode> barcodeMap;
};

}