// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
/*
// Boost.Test include(s).
#include <boost/test/unit_test.hpp>

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

using namespace Acts::TracccConversion;

BOOST_AUTO_TEST_CASE(Traccc_Conversion_Chain) {
    // Memory resource used by the application.
    vecmem::host_memory_resource host_mr;

    const std::string detectorFile = "/home/frederik/Downloads/traccc-data-v6/tml_detector/trackml-detector.csv";
    const std::string digitalizationFile = "/home/frederik/Downloads/traccc-data-v6/tml_detector/default-geometric-config-generic.json";
    const std::string inputDirectory = "/home/frederik/Downloads/traccc-data-v6/tml_pixels/";
    const std::string eventFile = inputDirectory + "event000000000-cells.csv";

    traccc::data_format format = traccc::data_format::csv;

    // Read in the geometry.
    //auto [surface_transforms, barcode_map] = traccc::io::read_geometry(
    //        detectorFile, traccc::data_format::json);
    
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

    //auto chain = factory.buildChainHost(host_mr, detector, field, finderConfig, gridConfig, filterConfig, findingConfig, fittingConfig);

    //traccc::io::cell_reader_output readOut(&host_mr);

    // Read the cells from the relevant event file
    //traccc::io::read_cells(readOut, eventFile, format, &surface_transforms, &digitizationConfiguration, barcode_map.get());

    // container_types<fitting_result<transform3>, track_state<transform3>>
    //traccc::track_state_container_types::host result = chain.run(readOut.cells, readOut.modules);
    //for (std::size_t i = 0; i < result.size(); i++){
    //    auto c = result[i];
    //    auto fittingResult = c.header;
    //    auto trackStates = c.items;
    //}
    
}*/

// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/utils/algorithm.hpp
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/container.hpp#L154
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/details/host_container.hpp#L25
// https://github.com/acts-project/traccc/blob/f51a1f8c4031fae664c752695e177b8315fcf22b/core/include/traccc/edm/details/container_base.hpp#L42
// https://github.com/acts-project/traccc/blob/f51a1f8c4031fae664c752695e177b8315fcf22b/core/include/traccc/edm/details/container_element.hpp#L28

// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/track_state.hpp