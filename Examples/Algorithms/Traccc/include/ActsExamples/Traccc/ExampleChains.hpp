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

namespace ActsExamples::ExampleChains {

auto buildHost(vecmem::host_memory_resource& mr){
    using detector_t = detray::detector<detray::default_metadata, detray::host_container_types>;
    using stepper_t = detray::rk_stepper<typename detray::bfield::const_field_t::view_t, typename detector_t::algebra_type, detray::constrained_step<>>;
    using navigator_t = detray::navigator<const detector_t>;
    using finding_algorithm_t = traccc::finding_algorithm<stepper_t, navigator_t>;
    using fitter_t = traccc::kalman_fitter<stepper_t, navigator_t>;
    using fitting_algorithm_t = traccc::fitting_algorithm<fitter_t>;
    using clusterization_algorithm_t = traccc::host::clusterization_algorithm;
    using spacepoint_formation_algorithm_t = traccc::host::spacepoint_formation_algorithm;
    using seeding_algorithm_t = traccc::seeding_algorithm;
    using track_params_estimation_algorithm_t = traccc::track_params_estimation;
    using resolution_algorithm_t = traccc::greedy_ambiguity_resolution_algorithm;
    using chain_t = StandardChain<clusterization_algorithm_t, spacepoint_formation_algorithm_t, seeding_algorithm_t, track_params_estimation_algorithm_t, finding_algorithm_t, fitting_algorithm_t, resolution_algorithm_t>;

    const traccc::seedfinder_config finderConfig;
    const traccc::spacepoint_grid_config gridConfig{finderConfig};
    const traccc::seedfilter_config filterConfig;
    const typename finding_algorithm_t::config_type findingConfig;
    const typename fitting_algorithm_t::config_type fittingConfig;

    // Algorithms
    traccc::host::clusterization_algorithm ca(mr);
    traccc::host::spacepoint_formation_algorithm sf(mr);
    traccc::seeding_algorithm sa(finderConfig, gridConfig, filterConfig, mr);
    traccc::track_params_estimation tp(mr);
    finding_algorithm_t findAlg(findingConfig);
    fitting_algorithm_t fitAlg(fittingConfig);
    traccc::greedy_ambiguity_resolution_algorithm res;

    const chain_t chain{
        std::move(ca),
        std::move(sf),
        std::move(sa),
        std::move(tp),
        std::move(findAlg),
        std::move(fitAlg),
        std::move(res)
    };

    return chain;
}
}