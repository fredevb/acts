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
#include <sstream>

namespace ActsExamples::Traccc::Common::Debug {

namespace{

std::string toString(const Acts::ActsVector<2>& vec){
    std::stringstream ss;
    ss << "(" << vec[0] << ", " << vec[1] << ")";
    return ss.str();
}

// Structure to hold table data
struct MeasurementMatchRow {
    std::size_t idx1;
    Acts::ActsVector<2> local1;
    Acts::ActsVector<2> variance1;
    std::size_t idx2;
    Acts::ActsVector<2> local2;
    Acts::ActsVector<2> variance2;
    Acts::ActsScalar distanceLocal;
};

auto getTable(const std::vector<Acts::BoundVariantMeasurement>& measurements1, const std::vector<Acts::BoundVariantMeasurement>& measurements2, const std::map<std::size_t, std::size_t>& indexMap){
    std::vector<MeasurementMatchRow> table;
    for (std::size_t idx1 = 0; idx1 < measurements1.size(); ++idx1){
        MeasurementMatchRow row;
        auto measurement1 = measurements1[idx1];
        row.idx1 = idx1;
        row.local1 = Acts::TracccPlugin::getLocal(measurement1);
        row.variance1 = Acts::TracccPlugin::getVariance(measurement1);

        auto idx2 = indexMap.at(idx1);
        auto measurement2 = measurements2[idx2];
        row.idx2 = idx2;
        row.local2 = Acts::TracccPlugin::getLocal(measurement2);
        row.variance2 = Acts::TracccPlugin::getVariance(measurement2);

        row.distanceLocal = (row.local1 - row.local2).norm();
        table.push_back(row);
    }   
    return table;
}

}

// Function to print a formatted table
auto pairingStatistics(const std::vector<Acts::BoundVariantMeasurement>& measurements1, const std::vector<Acts::BoundVariantMeasurement>& measurements2, const std::map<std::size_t, std::size_t>& indexMap) {
    auto table = getTable(measurements1, measurements2, indexMap);
    
    std::stringstream ss;
    // Column headers
    ss  << std::setw(6) << "Idx1"
        << std::setw(25) << "Local1"
        << std::setw(35) << "Variance1"
        << std::setw(20) << "Idx2"
        << std::setw(25) << "Local2"
        << std::setw(35) << "Variance2"
        << std::setw(35) << "Distance Local"
        << std::endl;

    // Line separator
    ss << std::string(182, '-') << std::endl;

    // Print each row
    for (const auto& row : table) {
        ss  << std::setw(6) << row.idx1
            << std::setw(25) << toString(row.local1)
            << std::setw(35) << toString(row.variance1)
            << std::setw(20) << row.idx2
            << std::setw(25) << toString(row.local2)
            << std::setw(35) << toString(row.variance2)
            << std::setw(35) << std::fixed << std::setprecision(2) << row.distanceLocal
            << std::endl;
    }
    return ss.str();
}
}