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
#include "traccc/fitting/fitting_config.hpp"
#include "traccc/edm/cell.hpp"
#include "traccc/definitions/primitives.hpp"

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

namespace ActsExamples::Traccc::Cuda {

template <typename field_view_t>
struct Types {
    using DetectorType = detray::detector<detray::default_metadata, detray::device_container_types>;
    using StepperType = detray::rk_stepper<field_view_t, typename DetectorType::algebra_type, detray::constrained_step<>>;
    using NavigatorType = detray::navigator<const DetectorType>;
    using FitterType = traccc::kalman_fitter<StepperType, NavigatorType>;

    using ClusterizationAlgorithmType = traccc::cuda::clusterization_algorithm;
    using SpacepointFormationAlgorithmType = traccc::cuda::spacepoint_formation_algorithm;
    using SeedingAlgorithmType = traccc::cuda::seeding_algorithm;
    using TrackParametersEstimationAlgorithmType = traccc::cuda::track_params_estimation;
    using FindingAlgorithmType = traccc::finding_algorithm<StepperType, NavigatorType>;
    using FittingAlgorithmType = traccc::fitting_algorithm<FitterType>;
};

}
