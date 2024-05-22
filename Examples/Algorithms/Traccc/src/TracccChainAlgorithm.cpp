// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Acts examples include(s)
#include "ActsExamples/Traccc/TracccChainAlgorithm.hpp"

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
  //if (m_cfg.field == nullptr) {
  //  throw std::invalid_argument("Missing field");
  //}
  if (m_cfg.trackingGeometry == nullptr) {
    throw std::invalid_argument("Missing track geometry");
  }

  if (m_cfg.digitizationConfigs.empty()) {
    throw std::invalid_argument("Missing digitization configuration");
  }

  m_inputCells.initialize(m_cfg.inputCells);
  m_outputTracks.initialize(m_cfg.outputTracks);

  chainAdapter = std::make_shared<const ChainAdapter>(m_cfg.trackingGeometry, m_cfg.digitizationConfigs, *m_cfg.field);
}


ActsExamples::ProcessCode ActsExamples::TracccChainAlgorithm::execute(
const AlgorithmContext& ctx) const {
  // Read input data
  const auto cellsMap = m_inputCells(ctx);


  //TODO: Compile traccc chain here!
  //Remember to update the hard coded filepath that loads the detray detector
  //Make getCellTime and getCellActivation
  //Convert to correct Acts::GeometryHeirarchy input
  //TODO: Load input cellsMap
  //SETUP WRAPPED CHAIN

  // Get the geometry.
  //auto [surface_transforms, barcode_map] = getGeometry(detector);

  //const std::string inputDirectory = "/home/frederik/Desktop/CERN-TECH/input/";
 // const std::string digitalizationFile = "/home/frederik/Desktop/CERN-TECH-OTHER/traccc/data/tml_detector/default-geometric-config-generic.json";
 // const std::string eventFile = inputDirectory + "tml_pixels/event000000000-cells.csv";
  //traccc::data_format format = traccc::data_format::csv;

  // Read the digitization configuration file

  //Error
  //auto digitizationConfiguration = traccc::io::read_digitization_config(digitalizationFile);


  //traccc::io::cell_reader_output readOut(&host_mr);

  // Read the cells from the relevant event file
  //traccc::io::read_cells(readOut, eventFile, format, &surface_transforms, &digitizationConfiguration, barcode_map.get());


  vecmem::host_memory_resource mr;

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

  const chain_t chain(
    std::move(ca),
    std::move(sf),
    std::move(sa),
    std::move(tp),
    std::move(findAlg),
    std::move(fitAlg),
    std::move(res));

  auto result = chainAdapter->runChain(cellsMap, chain, &mr); 

  ACTS_INFO("Ran the traccc algorithm!");
  
  std::stringstream ss;
  result.trackStateContainer().statistics().toStream(ss);
  ACTS_INFO(ss.str());

  m_outputTracks(ctx, std::move(result));
  std::cout << "Wrote Container" << std::endl;
  return ActsExamples::ProcessCode::SUCCESS;
}