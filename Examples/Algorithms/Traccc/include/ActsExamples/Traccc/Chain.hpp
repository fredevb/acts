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
#include "ActsExamples/Traccc/StandardChain.hpp"

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
    using memory_resource_type = vecmem::host_memory_resource;
    using detector_type = detray::detector<detray::default_metadata, detray::host_container_types>;
    using stepper_type = detray::rk_stepper<typename detray::bfield::const_field_t::view_t, typename detector_type::algebra_type, detray::constrained_step<>>;
    using navigator_type = detray::navigator<const detector_type>;
    using finding_algorithm_type = traccc::finding_algorithm<stepper_type, navigator_type>;
    using fitter_type = traccc::kalman_fitter<stepper_type, navigator_type>;
    using fitting_algorithm_type = traccc::fitting_algorithm<fitter_type>;
    using clusterization_algorithm_type = traccc::host::clusterization_algorithm;
    using spacepoint_formation_algorithm_type = traccc::host::spacepoint_formation_algorithm;
    using seeding_algorithm_type = traccc::seeding_algorithm;
    using track_parameters_estimation_algorithm_type = traccc::track_params_estimation;
    using ambiguity_resolution_algorithm_type = traccc::greedy_ambiguity_resolution_algorithm;
};

template <typename platform_t>
struct Chain {
    using platform_type = platform_t;

    typename platform_type::clusterization_algorithm_type clusterizationAlgorithm;
    typename platform_type::spacepoint_formation_algorithm_type spacepointFormationAlgorithm;
    typename platform_type::seeding_algorithm_type seedingAlgorithm;
    typename platform_type::track_parameters_estimation_algorithm_type trackParametersEstimationAlgorithm;
    typename platform_type::finding_algorithm_type findingAlgorithm;
    typename platform_type::fitting_algorithm_type fittingAlgorithm;
    typename platform_type::ambiguity_resolution_algorithm_type ambiguityResolutionAlgorithm;

    struct Config{
        const traccc::seedfinder_config seedfinderConfig;
        const traccc::spacepoint_grid_config spacepointGridConfig{seedfinderConfig};
        const traccc::seedfilter_config seedfilterConfig;
        const typename platform_type::finding_algorithm_type::config_type findingConfig;
        const typename platform_type::fitting_algorithm_type::config_type fittingConfig;
        const typename platform_type::ambiguity_resolution_algorithm_type::config_t ambiguityResolutionConfig;
    };

    Chain(typename platform_type::memory_resource_type* mr, const Config& config):
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
