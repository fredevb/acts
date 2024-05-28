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

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/propagator/rk_stepper.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <map>

namespace ActsExamples::Chain {

struct Host {
    using MemoryResourceType = vecmem::host_memory_resource;
    using DetectorType = detray::detector<detray::default_metadata, detray::host_container_types>;
    using StepperType = detray::rk_stepper<typename detray::bfield::const_field_t::view_t, typename DetectorType::algebra_type, detray::constrained_step<>>;
    using NavigatorType = detray::navigator<const DetectorType>;
    using FitterType = traccc::kalman_fitter<StepperType, NavigatorType>;

    using ClusterizationAlgorithmType = traccc::host::clusterization_algorithm;
    using SpacepointFormationAlgorithmType = traccc::host::spacepoint_formation_algorithm;
    using SeedingAlgorithmType = traccc::seeding_algorithm;
    using TrackParametersEstimationAlgorithmType = traccc::track_params_estimation;
    using FindingAlgorithmType = traccc::finding_algorithm<StepperType, NavigatorType>;
    using FittingAlgorithmType = traccc::fitting_algorithm<FitterType>;
    using AmbiguityResolutionAlgorithmType = traccc::greedy_ambiguity_resolution_algorithm;
};

template <typename platform_t>
struct Chain {
    using PlatformType = platform_t;

    typename PlatformType::ClusterizationAlgorithmType clusterizationAlgorithm;
    typename PlatformType::SpacepointFormationAlgorithmType spacepointFormationAlgorithm;
    typename PlatformType::SeedingAlgorithmType seedingAlgorithm;
    typename PlatformType::TrackParametersEstimationAlgorithmType trackParametersEstimationAlgorithm;
    typename PlatformType::FindingAlgorithmType findingAlgorithm;
    typename PlatformType::FittingAlgorithmType fittingAlgorithm;
    typename PlatformType::AmbiguityResolutionAlgorithmType ambiguityResolutionAlgorithm;

    struct Config{
        const traccc::seedfinder_config seedfinderConfig;
        const traccc::spacepoint_grid_config spacepointGridConfig{seedfinderConfig};
        const traccc::seedfilter_config seedfilterConfig;
        const typename PlatformType::FindingAlgorithmType::config_type findingConfig;
        const typename PlatformType::FittingAlgorithmType::config_type fittingConfig;
        const typename PlatformType::AmbiguityResolutionAlgorithmType::config_t ambiguityResolutionConfig;
    };

    Chain(typename PlatformType::MemoryResourceType* mr, const Config& config):
        clusterizationAlgorithm(*mr),
        spacepointFormationAlgorithm(*mr),
        seedingAlgorithm(config.seedfinderConfig, config.spacepointGridConfig, config.seedfilterConfig, *mr),
        trackParametersEstimationAlgorithm(*mr),
        findingAlgorithm(config.findingConfig),
        fittingAlgorithm(config.fittingConfig),
        ambiguityResolutionAlgorithm(config.ambiguityResolutionConfig)
    {}
};

}
