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

  using elem_t = std::pair<Acts::GeometryIdentifier, Acts::BinUtility>;
  
  std::vector<elem_t> vec;
  for (auto& e : m_cfg.digitizationConfigs.getElements()){
    auto geoID = e.first;
    auto segs = e.second.geometricDigiConfig.segmentation;
    vec.push_back({geoID, segs});
  }
  std::shared_ptr<const Acts::GeometryHierarchyMap<Acts::BinUtility>> segmentations = std::make_shared<const Acts::GeometryHierarchyMap<Acts::BinUtility>>(vec);

  wrappedChain = std::make_shared<WrappedChain<chain_t>>(m_cfg.trackingGeometry, m_cfg.digitizationConfigs, *m_cfg.field, &mr, chain);
  // Load geosegs
}


ActsExamples::ProcessCode ActsExamples::TracccChainAlgorithm::execute(
const AlgorithmContext& ctx) const {
  // Read input data
  const auto& cellsMap = m_inputCells(ctx);


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

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  (*wrappedChain)(tracks, cellsMap); 

  ACTS_INFO("Ran the traccc algorithm!");

  std::stringstream ss;
  trackStateContainer->statistics().toStream(ss);
  ACTS_INFO(ss.str());

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};
  
  std::cout << "Created Container" << std::endl;

  m_outputTracks(ctx, std::move(constTracks));
  std::cout << "Wrote Container" << std::endl;
  return ActsExamples::ProcessCode::SUCCESS;
}