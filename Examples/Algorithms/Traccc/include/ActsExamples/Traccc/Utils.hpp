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

namespace ActsExamples::TracccPluginUtils {

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

float getTime(const Cluster::Cell& /*cell*/){
    return 0.f;
}

float getActivation(const Cluster::Cell& cell){
    return static_cast<float>(cell.activation);
}

unsigned int getRow(const Cluster::Cell& cell){
    if (cell.bin[0] > UINT_MAX) {
        throw std::runtime_error("Overflow will occur when casting to unsigned int.");
    }
    return static_cast<unsigned int>(cell.bin[0]);
}

unsigned int getColumn(const Cluster::Cell& cell){
        if (cell.bin[0] > UINT_MAX) {
        throw std::runtime_error("Overflow will occur when casting to unsigned int.");
    }
    return static_cast<unsigned int>(cell.bin[1]);
}

Acts::BinUtility getSegmentation(const DigiComponentsConfig& dcc){
    return dcc.geometricDigiConfig.segmentation;
}


template <typename detector_t, typename allocator_t>
auto getIndexSourceLinkAndMeasurements(const detector_t& detector, const std::vector<traccc::measurement, allocator_t>& measurements){
    IndexSourceLinkContainer sourceLinks;
    MeasurementContainer measurementContainer;

    for (const traccc::measurement& m : measurements) 
    {
        Acts::GeometryIdentifier moduleGeoId(detector.surface(m.surface_link).source);
        Index measurementIdx = measurementContainer.size();
        IndexSourceLink sourceLink{moduleGeoId, measurementIdx};
        sourceLinks.insert(sourceLinks.end(), sourceLink);
        measurementContainer.push_back(Acts::TracccPlugin::boundVariantMeasurement(m, sourceLink));
    }
    return std::make_tuple(std::move(sourceLinks), std::move(measurementContainer));
}

template <typename track_container_t, typename trajectory_t, template <typename> class holder_t>
void setUncalibratedSourceLinks(
    Acts::TrackContainer<track_container_t, trajectory_t, holder_t>& trackContainer,
    const IndexSourceLinkContainer& sourceLinks){
    
    using IndexType = typename Acts::TrackContainer<track_container_t, trajectory_t, holder_t>::IndexType;

    auto& tsc = trackContainer.trackStateContainer();
    for (IndexType i = 0; i < tsc.size(); i++){
        const auto geoID = tsc.getTrackState(i).referenceSurface().geometryId();
        const auto& it = sourceLinks.find(geoID); // also needs to be the one that also has the correct measurement ID
        if (it != sourceLinks.end()){
            const IndexSourceLink& isl = *it;
            Acts::SourceLink sl{isl};
            tsc.getTrackState(i).setUncalibratedSourceLink(std::move(sl));
            tsc.getTrackState(i).setCalibrated
        }
        else{
            auto msg = "Source link not found for geometry ID " + std::to_string(geoID.value());
            throw std::runtime_error(msg.c_str());
        }
    }
}

template <typename track_container_t, typename trajectory_t, template <typename> class holder_t>
void setCalibrated(
    Acts::TrackContainer<track_container_t, trajectory_t, holder_t>& trackContainer,
    const IndexSourceLinkContainer& sourceLinks){
    
    using IndexType = typename Acts::TrackContainer<track_container_t, trajectory_t, holder_t>::IndexType;

    auto& tsc = trackContainer.trackStateContainer();
    for (IndexType i = 0; i < tsc.size(); i++){
        const auto geoID = tsc.getTrackState(i).referenceSurface().geometryId();
        const auto& it = sourceLinks.find(geoID); // the one that also has the correct measurement ID
        if (it != sourceLinks.end()){
            const IndexSourceLink& isl = *it;
            Acts::SourceLink sl{isl};
            tsc.getTrackState(i).setUncalibratedSourceLink(std::move(sl));
            tsc.getTrackState(i).setCalibrated
        }
        else{
            auto msg = "Source link not found for geometry ID " + std::to_string(geoID.value());
            throw std::runtime_error(msg.c_str());
        }
    }
}


template <typename detector_t>
class Converter{

    public:
    
    Converter(
        std::shared_ptr<const Acts::TrackingGeometry> tg,
        std::shared_ptr<const detector_t> det,
        const Acts::GeometryHierarchyMap<DigiComponentsConfig>& dccMap):
        trackingGeometry(tg),
        detector(det),
        cellDataConverter(*detector, Acts::TracccPlugin::tracccConfig(dccMap, getSegmentation))
    {}

    auto convertInput(const std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>> map, vecmem::memory_resource* mr) const {
        auto tcm = Acts::TracccPlugin::tracccCellsMap(
            map, 
            getRow, 
            getColumn, 
            getActivation, 
            getTime
        );
        return cellDataConverter(mr, tcm);
    }

    template <typename tracks_t, typename allocator_t>
    auto convertOutput(const std::vector<traccc::measurement, allocator_t>& ms, const tracks_t& ts) const {
        auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
        auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
        TrackContainer tracks(trackContainer, trackStateContainer);

        Acts::TracccPlugin::copyTrackContainer(ts, tracks, *detector, *trackingGeometry);

        auto [sourceLinks, measurementContainer] = getIndexSourceLinkAndMeasurements(*detector, ms);
        setSourceLinks(tracks, sourceLinks);

        ConstTrackContainer constTracks{
            std::make_shared<Acts::ConstVectorTrackContainer>(std::move(*trackContainer)),
            std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*trackStateContainer))
            };
        return constTracks;
    }

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const detector_t> detector;
    const Acts::TracccPlugin::CellDataConverter cellDataConverter;
};

}