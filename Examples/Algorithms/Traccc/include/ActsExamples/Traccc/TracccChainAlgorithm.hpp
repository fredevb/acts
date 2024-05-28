// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Traccc plugin include(s)
#include "ActsExamples/Traccc/Chain.hpp"
#include "ActsExamples/Traccc/Converter.hpp"

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
template <typename platform_t>
class TracccChainAlgorithm final : public IAlgorithm {
private:

// Temporarily used to get the corresponding detray detector for the Acts geometry.
// Will be replaced when the detray plugin containing geometry conversion is complete.
inline auto readDetector(vecmem::memory_resource* mr, const std::string& detectorFilePath, const std::string& materialFilePath = "", const std::string& gridFilePath = "")
{
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
    auto [det, names] = detray::io::read_detector<DetectorType>(*mr, cfg);
    return std::move(det);
}

public:

using PlatformType = platform_t;
using ChainType = Chain::Chain<PlatformType>;
using DetectorType = typename PlatformType::DetectorType;
using FieldType = Acts::CovfieConversion::constant_field_t;

using CellsMapType = std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>>;

struct Config {
    std::string inputCells;
    std::string outputTracks;
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    std::shared_ptr<const Acts::ConstantBField> field;
    Acts::GeometryHierarchyMap<DigiComponentsConfig> digitizationConfigs;
    std::shared_ptr<const typename ChainType::Config> chainConfig;
};

/// Construct the traccc algorithm.
///
/// @param cfg is the algorithm configuration
/// @param lvl is the logging level
TracccChainAlgorithm(
    Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::IAlgorithm("TracccChainAlgorithm", lvl),
      m_cfg(std::move(cfg)) {

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

  m_inputCells.initialize(m_cfg.inputCells);
  m_outputTracks.initialize(m_cfg.outputTracks);

  // Read detector from file as temporary solution until detray plugin is complete.
  detector = std::make_shared<const DetectorType>(readDetector(&platformMemoryResource, "/home/frederik/Desktop/CERN-TECH/input/odd-detray_geometry_detray.json"));
  field = std::make_shared<const FieldType>(Acts::CovfieConversion::covfieField(*m_cfg.field));
  chain = std::make_shared<const ChainType>(&platformMemoryResource, *m_cfg.chainConfig);
  dataConverter = std::make_shared<const ActsExamples::TracccConversion::Converter<DetectorType>>(m_cfg.trackingGeometry, detector, m_cfg.digitizationConfigs);
}

/// Run the algorithm.
///
/// @param ctx is the algorithm context with event information
/// @return a process code indication success or failure
ProcessCode execute(const AlgorithmContext& ctx) const final override;

/// Const access to the config
const Config& config() const { return m_cfg; }

private:

Config m_cfg;

ReadDataHandle<CellsMapType> m_inputCells{this, "InputCells"};

WriteDataHandle<ConstTrackContainer> m_outputTracks{this, "OutputTracks"};

// Ensure order of destructor call.
typename PlatformType::MemoryResourceType platformMemoryResource;
std::shared_ptr<const DetectorType> detector;
std::shared_ptr<const ChainType> chain;
std::shared_ptr<const FieldType> field;
std::shared_ptr<const TracccConversion::Converter<DetectorType>> dataConverter;

};

}  // namespace ActsExamples
