// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccChainAlgorithm.hpp"
#include "Acts/Plugins/Traccc/TracccConversion.hpp"

// acts include (s)
#include "Acts/Utilities/Result.hpp"

// traccc plugin
#include "Acts/Plugins/Traccc/TracccConversion.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackProxy.hpp"

// tracc include (s)
// algorithms
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation_algorithm.hpp"
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"

// configs
#include "traccc/finding/finding_config.hpp"

// io
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_digitization_config.hpp"
#include "traccc/io/read_geometry.hpp"
#include "traccc/io/utils.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/navigation/navigator.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <type_traits>
#include <tuple>

#include <stdexcept>

using namespace Acts::TracccConversion;


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

  m_inputCells.initialize(m_cfg.inputCells);
  m_outputTracks.initialize(m_cfg.outputTracks);

}

ActsExamples::ProcessCode ActsExamples::TracccChainAlgorithm::execute(
const AlgorithmContext& ctx) const {
  // Read input data
  const auto& cell = m_inputCells(ctx);

  // Memory resource used by the application.
  vecmem::host_memory_resource host_mr;

  const std::string detectorFile = "/home/frederik/Downloads/traccc-data-v6/tml_detector/trackml-detector.csv";
  const std::string digitalizationFile = "/home/frederik/Downloads/traccc-data-v6/tml_detector/default-geometric-config-generic.json";
  const std::string inputDirectory = "/home/frederik/Downloads/traccc-data-v6/tml_pixels/";
  const std::string eventFile = inputDirectory + "event000000000-cells.csv";

  traccc::data_format format = traccc::data_format::csv;

  auto detector = readDetector(host_mr, detectorFile);

  // Get the geometry.
  auto [surface_transforms, barcode_map] = getGeometry(detector);

  using field_t = detray::bfield::const_field_t;
  const traccc::vector3 field_vec = {0.f, 0.f, 1.0f};
  const field_t field = detray::bfield::create_const_field(field_vec);

  // Read the digitization configuration file
  auto digitizationConfiguration = traccc::io::read_digitization_config(digitalizationFile);

  TracccChainFactory<decltype(detector)> factory;

  const traccc::seedfinder_config finderConfig{};
  const traccc::spacepoint_grid_config gridConfig{finderConfig};
  const traccc::seedfilter_config filterConfig{};
  const typename decltype(factory)::finding_algorithm_t::config_type findingConfig{};
  const typename decltype(factory)::fitting_algorithm_t::config_type fittingConfig{};

  auto chain = factory.buildChainHost(host_mr, detector, field, finderConfig, gridConfig, filterConfig, findingConfig, fittingConfig);

  traccc::io::cell_reader_output readOut(&host_mr);

  // Read the cells from the relevant event file
  traccc::io::read_cells(readOut, eventFile, format, &surface_transforms, &digitizationConfiguration, barcode_map.get());

  auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
  auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
  TrackContainer tracks(trackContainer, trackStateContainer);

  chain.run(readOut.cells, readOut.modules, tracks);

  ACTS_INFO("Ran the traccc algorithm!");

  std::stringstream ss;
  trackStateContainer->statistics().toStream(ss);
  ACTS_INFO(ss.str());

  ConstTrackContainer constTracks{
      std::make_shared<Acts::ConstVectorTrackContainer>(
          std::move(*trackContainer)),
      std::make_shared<Acts::ConstVectorMultiTrajectory>(
          std::move(*trackStateContainer))};

  m_outputTracks(ctx, std::move(constTracks));
  return ActsExamples::ProcessCode::SUCCESS;
}