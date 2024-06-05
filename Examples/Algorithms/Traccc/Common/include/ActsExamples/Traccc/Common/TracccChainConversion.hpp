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

// Boost include(s)
#include <boost/range/combine.hpp>

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <map>

namespace ActsExamples::Traccc::Common {

/// @brief Gets the time of the cell.
/// @note Currently, it always returns 0.
inline float getTime(const Cluster::Cell& /*cell*/){
    return 0.f;
}

/// @brief Gets the activation of the cell.
inline float getActivation(const Cluster::Cell& cell){
    return static_cast<float>(cell.activation);
}

/// @brief Gets the row of the cell.
inline unsigned int getRow(const Cluster::Cell& cell){
    if (cell.bin[0] > UINT_MAX) {
        throw std::runtime_error("Overflow will occur when casting to unsigned int.");
    }
    return static_cast<unsigned int>(cell.bin[0]);
}

/// @brief Gets the column of the cell.
inline unsigned int getColumn(const Cluster::Cell& cell){
        if (cell.bin[0] > UINT_MAX) {
        throw std::runtime_error("Overflow will occur when casting to unsigned int.");
    }
    return static_cast<unsigned int>(cell.bin[1]);
}

inline Acts::BinUtility getSegmentation(const DigiComponentsConfig& dcc){
    return dcc.geometricDigiConfig.segmentation;
}

/// @brief Converts a "geometry ID -> generic cell collection type" map to a "geometry ID -> traccc cell collection" map.
/// @note The function sets the module link of the cells in the output to 0.
/// @return Map from geometry ID to its cell data (as a vector of traccc cell data)
template <typename cell_collection_t>
inline std::map<std::uint64_t, std::vector<traccc::cell>> tracccCellsMap(
    const std::map<Acts::GeometryIdentifier, cell_collection_t>& map)
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
inline traccc::digitization_config tracccConfig(
    const Acts::GeometryHierarchyMap<data_t>& config,
    const get_segmentation_fn_t& getSegmentation){
    using ElementType = std::pair<Acts::GeometryIdentifier, traccc::module_digitization_config>;
    std::vector<ElementType> vec;
    for (auto& e : config.getElements()){
        vec.push_back({e.first, traccc::module_digitization_config{getSegmentation(e.second)}});
    }
    return traccc::digitization_config(vec);
}

/// @brief This class provides functions for converting input and output data using the data structures in Acts Examples.
/// @tparam detector_t The type of the detray detector this converter expects.
template <typename detector_t>
class TracccChainDataConverter{

    public:
    
    TracccChainDataConverter(
        const Acts::TrackingGeometry& tg,
        const detector_t& det,
        const Acts::GeometryHierarchyMap<DigiComponentsConfig>& dccMap):
        trackingGeometry(tg),
        detector(det),
        digitizationConfig(tracccConfig(dccMap, getSegmentation)),
        surfaceTransforms(traccc::io::alt_read_geometry(det)),
        barcodeMap(Acts::TracccPlugin::createBarcodeMap(det))
    {}

    /// @brief Converts a map of cells to traccc input
    /// @param map The geometry id -> cell collection map which corresponds to each geometry's cells.
    /// @param mr The memory resource to use.
    /// @returns The converted input as a tuple containing traccc input data, i.e. (cells, modules).
    auto createTracccInput(const std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>& map, vecmem::memory_resource* mr) const {
        auto tcm = tracccCellsMap(map);
        return Acts::TracccPlugin::createCellsAndModules(
            mr, tcm, &surfaceTransforms, &digitizationConfig, &barcodeMap);
    }

    /// @brief Converts a container of traccc tracks to a container of Acts tracks.
    /// (this is needed to set the source links and calibrated data of the newly converted track states).
    /// @param tracccTrackContainer The (indexable) container of traccc tracks.
    /// @param sourceLinks the source links used to set the uncalibrated source links.
    /// The source link is set by indexing the container using the traccc track state measurement id.
    /// @param measurements the measurements used to set the calibrated data.
    /// The measurement is set by indexing the container using the traccc track state measurement id.
    /// @return An Acts const track container containing the same track data as the traccc tracks container.
    /// Furthermore, the track states sourcelinks and measurements are also set.
    template <typename traccc_track_container_t>
    auto createActsTracks(
        const traccc_track_container_t& tracccTrackContainer,
        const std::map<traccc::measurement, Acts::BoundVariantMeasurement>& measurementConversionMap
    ) const {
        auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
        auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
        TrackContainer tracks(trackContainer, trackStateContainer);

        for (std::size_t i = 0; i < tracccTrackContainer.size(); i++) {
            auto ttrack = tracccTrackContainer[i];
            auto atrack = Acts::TracccPlugin::makeTrack(ttrack, tracks, detector, trackingGeometry);

            Acts::TracccPlugin::setSourceAndMeasurements(ttrack, atrack, measurementConversionMap);
        }

        ConstTrackContainer constTracks{
            std::make_shared<Acts::ConstVectorTrackContainer>(std::move(*trackContainer)),
            std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*trackStateContainer))
        };

        return constTracks;
    }
    
    const Acts::TrackingGeometry& trackingGeometry;
    const detector_t& detector;

    // Cache the converted digitalization configuration, the surface transforms, and the barcode map.
    const traccc::digitization_config digitizationConfig;
    const traccc::geometry surfaceTransforms;
    const std::map<std::uint64_t, detray::geometry::barcode> barcodeMap;
};

}