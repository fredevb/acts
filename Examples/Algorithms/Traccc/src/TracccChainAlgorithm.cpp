// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Acts examples include(s)
#include "ActsExamples/Traccc/TracccChainAlgorithm.hpp"
#include "ActsExamples/Traccc/ExampleChains.hpp"

// Acts include (s)
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/Utilities/Result.hpp"

// Traccc include (s)
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation_algorithm.hpp"
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/finding/finding_config.hpp"
// #include "traccc/io/data_format.hpp"
// #include "traccc/io/read_digitization_config.hpp"

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <type_traits>
#include <tuple>

#include <stdexcept>

using namespace Acts::TracccPlugin;

ActsExamples::TracccChainAlgorithm::TracccChainAlgorithm(
    ActsExamples::TracccChainAlgorithm::Config cfg, Acts::Logging::Level lvl)
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

  detector_mr = std::make_shared<vecmem::host_memory_resource>();
  detector = std::make_shared<const detector_t>(ActsExamples::TracccConversion::readDetector(detector_mr.get(), "/home/frederik/Desktop/CERN-TECH/input/odd-detray_geometry_detray.json"));
  field = std::make_shared<const Acts::CovfieConversion::constant_field_t>(Acts::CovfieConversion::covfieField(*m_cfg.field));

  chainRunner = std::make_shared<const StandardChainRunner>(m_cfg.trackingGeometry, detector, m_cfg.digitizationConfigs);
}


ActsExamples::ProcessCode ActsExamples::TracccChainAlgorithm::execute(
const AlgorithmContext& ctx) const {

  const auto cellsMap = m_inputCells(ctx);

  vecmem::host_memory_resource mr;

  auto chain = ExampleChains::buildHost(mr);

  auto result = chainRunner->run(cellsMap, *field, chain); 

  ACTS_INFO("Ran the traccc algorithm!");
  
  std::stringstream ss;
  result.trackStateContainer().statistics().toStream(ss);
  ACTS_INFO(ss.str());

  m_outputTracks(ctx, std::move(result));
  std::cout << "Wrote Container" << std::endl;
  return ActsExamples::ProcessCode::SUCCESS;
  
}