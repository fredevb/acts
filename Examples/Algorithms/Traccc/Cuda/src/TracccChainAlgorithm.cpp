// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Acts examples include(s)
#include "ActsExamples/Traccc/TracccChainAlgorithm.hpp"

// VecMem include(s).
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>


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

template <>
ActsExamples::ProcessCode ActsExamples::TracccChainAlgorithm<ActsExamples::Chain::Cuda>::execute(
const ActsExamples::AlgorithmContext& ctx) const {

      using DetectorType = detray::detector<detray::default_metadata, detray::device_container_types>;
    using StepperType = detray::rk_stepper<typename detray::bfield::const_field_t::view_t, typename DetectorType::algebra_type, detray::constrained_step<>>;
    using NavigatorType = detray::navigator<const DetectorType>;
    using FitterType = traccc::kalman_fitter<StepperType, NavigatorType>;

    using ClusterizationAlgorithmType = traccc::cuda::clusterization_algorithm;
    using SpacepointFormationAlgorithmType = traccc::cuda::spacepoint_formation_algorithm;
    using SeedingAlgorithmType = traccc::cuda::seeding_algorithm;
    using TrackParametersEstimationAlgorithmType = traccc::cuda::track_params_estimation;
    using FindingAlgorithmType = traccc::finding_algorithm<StepperType, NavigatorType>;
    using FittingAlgorithmType = traccc::fitting_algorithm<FitterType>;
    using AmbiguityResolutionAlgorithmType = traccc::greedy_ambiguity_resolution_algorithm;

  const auto cellsMap = m_inputCells(ctx);

  auto [cells, modules] = dataConverter->convertInput(cellsMap, &mr);

  ACTS_VERBOSE("Converted the Acts input data to traccc input data");

  // Create device copy of input collections
  cell_collection_types::buffer cells_buffer(cells.size(),
                                              *m_cached_device_mr);
  m_copy(vecmem::get_data(cells), cells_buffer)->ignore();
  cell_module_collection_types::buffer modules_buffer(modules.size(),
                                                      *m_cached_device_mr);
  m_copy(vecmem::get_data(modules), modules_buffer)->ignore();

  // Run the clusterization (asynchronously).
  const clusterization_algorithm::output_type measurements =
      m_clusterization(cells_buffer, modules_buffer);
  m_measurement_sorting(measurements);

  // Run the seed-finding (asynchronously).
  const spacepoint_formation_algorithm::output_type spacepoints =
      m_spacepoint_formation(measurements, modules_buffer);
  const track_params_estimation::output_type track_params =
      m_track_parameter_estimation(spacepoints, m_seeding(spacepoints),
                                    m_field_vec);

  // Create the buffer needed by track finding and fitting.
  auto navigation_buffer = detray::create_candidates_buffer(
      *m_detector,
      m_finding_config.navigation_buffer_size_scaler *
          m_copy.get_size(track_params),
      *m_cached_device_mr, &m_host_mr);

  // Run the track finding (asynchronously).
  const finding_algorithm::output_type track_candidates =
      m_finding(m_device_detector_view, m_field, navigation_buffer,
                measurements, track_params);

  // Run the track fitting (asynchronously).
  const fitting_algorithm::output_type track_states =
      m_fitting(m_device_detector_view, m_field, navigation_buffer,
                track_candidates);

  // Copy a limited amount of result data back to the host.
  output_type trackStates{&m_host_mr};
  m_copy(track_states.headers, trackStates)->wait();

  auto result = dataConverter->convertOutput(measurements, trackStates); 

  ACTS_VERBOSE("Converted the traccc track data to Acts track data");

  m_outputTracks(ctx, std::move(result));
  return ActsExamples::ProcessCode::SUCCESS;
}
