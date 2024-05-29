// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Acts examples include(s)
#include "ActsExamples/Traccc/Host/TracccChainAlgorithm.hpp"

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <memory>
#include <type_traits>
#include <tuple>

#include <stdexcept>

ActsExamples::Traccc::Host::TracccChainAlgorithm::TracccChainAlgorithm(
    Config cfg, Acts::Logging::Level lvl)
    : ActsExamples::Traccc::Common::TracccChainAlgorithmBase(std::move(cfg), std::move(lvl)),
    clusterizationAlgorithm(hostMemoryResource),
    spacepointFormationAlgorithm(hostMemoryResource),
    seedingAlgorithm(m_cfg.chainConfig->seedfinderConfig, m_cfg.chainConfig->spacepointGridConfig, m_cfg.chainConfig->seedfilterConfig, hostMemoryResource),
    trackParametersEstimationAlgorithm(hostMemoryResource),
    findingAlgorithm(m_cfg.chainConfig->findingConfig),
    fittingAlgorithm(m_cfg.chainConfig->fittingConfig),
    ambiguityResolutionAlgorithm(m_cfg.chainConfig->ambiguityResolutionConfig)
{}

ActsExamples::ProcessCode ActsExamples::Traccc::Host::TracccChainAlgorithm::execute(
const ActsExamples::AlgorithmContext& ctx) const {

  vecmem::host_memory_resource mr;

  typename HostTypes::ClusterizationAlgorithmType::output_type measurements{&mr};
  typename HostTypes::SpacepointFormationAlgorithmType::output_type spacepoints{&mr};
  typename HostTypes::SeedingAlgorithmType::output_type seeds{&mr};
  typename HostTypes::TrackParametersEstimationAlgorithmType::output_type params{&mr};
  typename HostTypes::FindingAlgorithmType::output_type trackCandidates{&mr};
  typename HostTypes::FittingAlgorithmType::output_type trackStates{&mr};
  typename HostTypes::AmbiguityResolutionAlgorithmType::output_type resolvedTrackStates{&mr};

  const auto cellsMap = m_inputCells(ctx);

  auto [cells, modules] = dataConverter.convertInput(cellsMap, &mr);

  ACTS_VERBOSE("Converted the Acts input data to traccc input data");

  measurements = clusterizationAlgorithm(vecmem::get_data(cells), vecmem::get_data(modules)); //HERE, note gcda error in python script at the start. (different behaviour on full clean build?)
  
  ACTS_INFO("Ran the clusterization algorithm");

  spacepoints = spacepointFormationAlgorithm(vecmem::get_data(measurements), vecmem::get_data(modules));
  
  ACTS_INFO("Ran the spacepoint formation algorithm");

  seeds = seedingAlgorithm(spacepoints);
  
  ACTS_INFO("Ran the seeding algorithm");

  const typename FieldType::view_t fieldView(field);

  static_assert(std::is_same<FieldType, typename detray::bfield::const_field_t>::value, "Currently, traccc expects a constant field.");
  params = trackParametersEstimationAlgorithm(spacepoints, seeds, fieldView.at(0.f,0.f,0.f));
  
  ACTS_INFO("Ran the parameters estimation algorithm");

  trackCandidates = findingAlgorithm(detector, field, measurements, params);
  
  ACTS_INFO("Ran the finding algorithm");

  trackStates = fittingAlgorithm(detector, field, trackCandidates);

  ACTS_INFO("Ran the fitting algorithm");

  resolvedTrackStates = ambiguityResolutionAlgorithm(trackStates);
  
  ACTS_INFO("Ran the ambiguity resolution algorithm");

  auto result = dataConverter.convertOutput(measurements, resolvedTrackStates); 

  ACTS_VERBOSE("Converted the traccc track data to Acts track data");

  m_outputTracks(ctx, std::move(result));
  return ActsExamples::ProcessCode::SUCCESS;
}
