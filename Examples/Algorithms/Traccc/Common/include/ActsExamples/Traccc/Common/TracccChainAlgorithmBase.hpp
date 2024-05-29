// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Traccc plugin include(s)
#include "ActsExamples/Traccc/Common/TracccChainConfig.hpp"
#include "ActsExamples/Traccc/Common/TracccChainConversion.hpp"

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

namespace ActsExamples::Traccc::Common {


class TracccChainAlgorithmBase : public IAlgorithm {
public:

using DetectorHostType = detray::detector<detray::default_metadata, detray::host_container_types>;
using FieldType = Acts::CovfieConversion::constant_field_t;
using CellsMapType = std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>;

struct Config {
    std::string inputCells;
    std::string outputTracks;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    std::shared_ptr<const Acts::ConstantBField> field;
    Acts::GeometryHierarchyMap<DigiComponentsConfig> digitizationConfigs;
    std::shared_ptr<const TracccChainConfig> chainConfig;
};

/// Construct the traccc algorithm.
///
/// @param cfg is the algorithm configuration
/// @param lvl is the logging level
TracccChainAlgorithmBase(Config cfg, Acts::Logging::Level lvl);

/// Const access to the config
const Config& config() const { return m_cfg; }

protected:

Config m_cfg;

ReadDataHandle<CellsMapType> m_inputCells{this, "InputCells"};
WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};

// Ensure order of destructor call.
vecmem::host_memory_resource hostMemoryResource;
const DetectorHostType detector;
const FieldType field;
const TracccChainDataConverter<DetectorHostType> dataConverter;
//std::shared_ptr<const DetectorHostType> detector;
//std::shared_ptr<const FieldType> field;
//std::shared_ptr<const TracccChainDataConverter<DetectorHostType>> dataConverter;

private:

void TestValidConfig(){
  if (m_cfg.inputCells.empty()) {
    throw std::invalid_argument("Missing input cells");
  }

  if (m_cfg.field == nullptr) {
    throw std::invalid_argument("Missing field");
  }

  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing track geometry");
  }

  if (m_cfg.digitizationConfigs.empty()) {
    throw std::invalid_argument("Missing digitization configuration");
  }
}

};

}  // namespace ActsExamples