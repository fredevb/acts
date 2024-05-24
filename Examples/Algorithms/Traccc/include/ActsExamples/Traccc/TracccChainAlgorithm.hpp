// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Traccc plugin include(s)
#include "ActsExamples/Traccc/StandardChain.hpp"

// Acts include(s)
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Definitions/Algebra.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <memory>
#include <string>

namespace ActsExamples {

/// Construct a traccc algorithm
class TracccChainAlgorithm final : public IAlgorithm {
public:
struct Config {
    std::string inputCells;
    std::string outputTracks;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    std::shared_ptr<const Acts::ConstantBField> field;
    Acts::GeometryHierarchyMap<DigiComponentsConfig> digitizationConfigs;
};

using field_t = Acts::CovfieConversion::constant_field_t;

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
//using chain_t = StandardChainHost<clusterization_algorithm_t, spacepoint_formation_algorithm_t, seeding_algorithm_t, track_params_estimation_algorithm_t, finding_algorithm_t, fitting_algorithm_t, resolution_algorithm_t>;

/// Construct the traccc algorithm.
///
/// @param cfg is the algorithm configuration
/// @param lvl is the logging level
TracccChainAlgorithm(Config cfg, Acts::Logging::Level lvl);

/// Run the algorithm.
///
/// @param ctx is the algorithm context with event information
/// @return a process code indication success or failure
ProcessCode execute(const AlgorithmContext& ctx) const final override;

/// Const access to the config
const Config& config() const { return m_cfg; }

private:

using CellsMap = std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>;

Config m_cfg;

ReadDataHandle<CellsMap> m_inputCells{this, "InputCells"};

WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};

// memory resources should be declared before detector to ensure order of destructor call.
std::shared_ptr<vecmem::memory_resource> detector_mr;
// Same goes for dectector. It should come before the chain runner.
std::shared_ptr<const detector_t> detector;

std::shared_ptr<const StandardChainRunner> chainRunner;

std::shared_ptr<const Acts::CovfieConversion::constant_field_t> field;
};

}  // namespace ActsExamples
