// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost.Test include(s).
#include <boost/test/unit_test.hpp>

// traccc plugin
#include "Acts/Plugins/Traccc/TracccConversion.hpp"

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



void convertData(){

}

template <typename detector_t,
          typename field_t,
          typename clusterization_func_t, 
          typename spacepoint_formation_func_t, 
          typename seeding_func_t, 
          typename track_params_estimation_func_t, 
          typename finding_func_t, 
          typename fitting_func_t,
          typename resolution_func_t>
class TracccChain{

    public:

    TracccChain(
        vecmem::memory_resource& hostmr,
        detector_t& det,
        field_t& f,
        clusterization_func_t& cFunc,
                spacepoint_formation_func_t& sfFunc,
                seeding_func_t& sFunc,
                track_params_estimation_func_t& tpeFunc,
                finding_func_t& findFunc,
                fitting_func_t& fitFunc,
                resolution_func_t& resFunc) 
            : 
            hostMemoryResource(hostmr),
            detector(det),
            field(f),
            clusterizationFunc(cFunc),
            spacepointFormationFunc(sfFunc),
            seedingFunc(sFunc),
            trackParameterEstimationFunc(tpeFunc),
            findingFunc(findFunc),
            fittingFunc(fitFunc),
            resolutionFunc(resFunc){}


    void operator()(const traccc::cell_collection_types::host& cells,
        const traccc::cell_module_collection_types::host& modules) const {

        typename clusterization_func_t::output_type measurements{&hostMemoryResource};
        typename spacepoint_formation_func_t::output_type spacepoints{&hostMemoryResource};
        typename seeding_func_t::output_type seeds{&hostMemoryResource};
        typename track_params_estimation_func_t::output_type params{&hostMemoryResource};
        typename finding_func_t::output_type track_candidates{&hostMemoryResource};
        typename finding_func_t::output_type track_states{&hostMemoryResource};
        typename resolution_func_t::output_type resolved_track_states{&hostMemoryResource};

        measurements = vecmem::get_data(clusterizationFunc(vecmem::get_data(cells), vecmem::get_data(modules)));
        spacepoints = spacepointFormationFunc(measurements, vecmem::get_data(modules));
        seeds = seedingFunc(spacepoints);
        
        const typename field_t::view_t fieldView(field);
        static_assert(std::is_same<field_t, typename detray::bfield::const_field_t>::value, "Currently, traccc expects a constant field.");
        params = trackParameterEstimationFunc(spacepoints, seeds, fieldView.at(0,0,0));

        track_candidates = findingFunc(detector, field, measurements, params);
        track_states = fittingFunc(detector, field, track_candidates);

        if (ambiguityResolution) {
            resolved_track_states = resolutionFunc(track_states);
        }
    }

    private:

    bool ambiguityResolution = true;

    vecmem::memory_resource& hostMemoryResource;
    detector_t& detector;
    field_t& field;
    clusterization_func_t clusterizationFunc;
    spacepoint_formation_func_t spacepointFormationFunc;
    seeding_func_t seedingFunc;
    track_params_estimation_func_t trackParameterEstimationFunc;
    finding_func_t findingFunc;
    fitting_func_t fittingFunc;
    resolution_func_t resolutionFunc;
};

template <typename detector_t>
class TracccChainFactory{

    TracccChainFactory(detector_t& det) : detector(det) {}
    
    using stepper_t = detray::rk_stepper<typename detray::bfield::const_field_t::view_t, typename detector_t::transform3, detray::constrained_step<>>;
    
    using navigator_t = detray::navigator<const detector_t>;
    
    using finding_algorithm_t = traccc::finding_algorithm<stepper_t, navigator_t>;
    
    using fitter_t = traccc::kalman_fitter<stepper_t, navigator_t>;
    using fitting_algorithm_t = traccc::fitting_algorithm<fitter_t>;

    template <typename field_t>
    auto buildChainHost(
                vecmem::host_memory_resource& host_mr,
                detector_t& detector,
                field_t& field,
                const traccc::seedfinder_config& finderConfig,
                const traccc::spacepoint_grid_config& gridConfig,
                const traccc::seedfilter_config& filterConfig,
                const typename finding_algorithm_t::config_type& findingConfig,
                const typename fitting_algorithm_t::config_type& fittingConfig){
        
        // Algorithms
        traccc::host::clusterization_algorithm ca(host_mr);
        traccc::host::spacepoint_formation_algorithm sf(host_mr);
        traccc::seeding_algorithm sa(finderConfig,
                                    gridConfig,
                                    filterConfig,
                                    host_mr);

        traccc::track_params_estimation tp(host_mr);
        finding_algorithm_t findFunc(findingConfig);
        fitting_algorithm_t fitFunc(fittingConfig);
        traccc::greedy_ambiguity_resolution_algorithm res;

        TracccChain tc(host_mr, detector, field, ca, sf, sa, tp, findFunc, fitFunc, res);
        return tc;
    }
};


/*void x(){
    traccc::geometry surface_transforms;
    std::unique_ptr<std::map<std::uint64_t, detray::geometry::barcode>>
        barcode_map;     
}*/

auto readDetector(vecmem::host_memory_resource& host_mr, const std::string& detectorFile, const std::string& materialFile = "", const std::string& gridFile = ""){
    using detector_t = detray::detector<detray::default_metadata,
                                        detray::host_container_types>;
    // Set up the detector reader configuration.
    detray::io::detector_reader_config cfg;
    cfg.add_file(traccc::io::data_directory() +
                detectorFile);
    if (materialFile.empty() == false) {
        cfg.add_file(traccc::io::data_directory() +
                    materialFile);
    }
    if (gridFile.empty() == false) {
        cfg.add_file(traccc::io::data_directory() +
                    gridFile);
    }

    // Read the detector.
    auto [det, names] = detray::io::read_detector<detector_t>(host_mr, cfg);
    return std::move(det);
}

#include <iostream>

BOOST_AUTO_TEST_CASE(Traccc_Conversion_Chain) {
    // Memory resource used by the application.
    vecmem::host_memory_resource host_mr;

    const auto detectorFile = "/home/frederik/Downloads/traccc-data-v6/tml_detector/trackml-detector.csv";
    const auto digitalizationFile = "/home/frederik/Downloads/traccc-data-v6/tml_detector/default-geometric-config-generic.json";
    const auto inputDirectory = "/home/frederik/Downloads/traccc-data-v6/tml_detector/tml_pixels/";

    // Read in the geometry.
    //auto [surface_transforms, barcode_map] = traccc::io::read_geometry(
    //        detectorFile, traccc::data_format::json);
    
    auto detector = readDetector(host_mr, detectorFile);

    const traccc::vector3 field_vec = {0.f, 0.f, 1.0f};
    const detray::bfield::const_field_t field =
        detray::bfield::create_const_field(field_vec);

    // Read the digitization configuration file
    auto digi_cfg =
        traccc::io::read_digitization_config(digitalizationFile);

    TracccChainFactory factory(detector);

    const traccc::seedfinder_config finderConfig{};
    const traccc::spacepoint_grid_config gridConfig{finderConfig};
    const traccc::seedfilter_config filterConfig{};
    const typename decltype(factory)::finding_algorithm_t::config_type findingConfig{};
    const typename decltype(factory)::fitting_algorithm_t::config_type fittingConfig{};

    auto chain = factory.buildChainHost(host_mr, detector, field, finderConfig, gridConfig, filterConfig, findingConfig, fittingConfig);
    
}