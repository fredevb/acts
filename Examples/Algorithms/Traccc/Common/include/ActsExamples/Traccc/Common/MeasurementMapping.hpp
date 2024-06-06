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
#include "ActsExamples/Traccc/Common/LSH.hpp"

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

namespace{

inline Acts::GeometryIdentifier getGeometryID(const Acts::BoundVariantMeasurement& measurement){
    return std::visit([](auto& m){
        return m.sourceLink().template get<ActsExamples::IndexSourceLink>().geometryId();
        }, measurement);
}

inline bool mEq(const Acts::BoundVariantMeasurement& measurement1, const Acts::BoundVariantMeasurement& measurement2, const double radius = .1){
    auto gidEq = getGeometryID(measurement1) == getGeometryID(measurement2);
    
    auto sqNorm = (Acts::TracccPlugin::getLocal(measurement1) - Acts::TracccPlugin::getLocal(measurement2)).squaredNorm();
    auto locEq = sqNorm <= radius * radius;

    return gidEq && locEq;
}

inline std::size_t mHash(const Acts::BoundVariantMeasurement& measurement){
    // The hash function can be optimized to reduce collisions.
    return static_cast<std::size_t>(getGeometryID(measurement).value());
}

const auto wrappedHash = std::function<std::size_t(const Acts::BoundVariantMeasurement&)>(mHash);
const auto wrappedEq = std::function<bool(const Acts::BoundVariantMeasurement&, const Acts::BoundVariantMeasurement&)>(
    [](const Acts::BoundVariantMeasurement& m1, const Acts::BoundVariantMeasurement& m2) { return mEq(m1, m2); });

}

namespace ActsExamples::Traccc::Common {

inline auto getMeasurementMatchingMap(const std::vector<Acts::BoundVariantMeasurement>& from, const std::vector<Acts::BoundVariantMeasurement>& to){
    return MatchingMap(from, to, wrappedHash, wrappedEq);
}

/// @brief Creates a source links and measurements. 
/// Creates the measurements by copying the data in the traccc measurements.
/// @param detector The detray detector
/// @param measurements The traccc measurements
/// @return A vector of Acts bound variant measurements.
template <typename detector_t, typename allocator_t>
inline auto createActsMeasurements(const detector_t& detector, const std::vector<traccc::measurement, allocator_t>& measurements){
    std::vector<Acts::BoundVariantMeasurement> measurementContainer;
    for (const traccc::measurement& m : measurements) 
    {
        Acts::GeometryIdentifier moduleGeoId(detector.surface(m.surface_link).source);
        Index measurementIdx = measurementContainer.size();
        IndexSourceLink idxSourceLink{moduleGeoId, measurementIdx};
        measurementContainer.push_back(Acts::TracccPlugin::boundVariantMeasurement(m, Acts::SourceLink{idxSourceLink}));
    }
    return measurementContainer;
}

template <typename detector_t, typename allocator_t>
std::map<traccc::measurement, Acts::BoundVariantMeasurement> measurementConversionMap(const detector_t& detector, const std::vector<traccc::measurement, allocator_t>& tracccMeasurements, const std::vector<Acts::BoundVariantMeasurement>& actsMeasurements){
    auto convertedMeasurements = createActsMeasurements(detector, tracccMeasurements);
    auto indexMap = getMeasurementMatchingMap(convertedMeasurements, actsMeasurements);
    return referenceMap(tracccMeasurements, actsMeasurements, indexMap);
}

}