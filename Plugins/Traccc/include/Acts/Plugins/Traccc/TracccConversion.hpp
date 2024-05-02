// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Boost.Test include(s).
#include <boost/test/unit_test.hpp>

// acts include (s)
#include "Acts/Utilities/Result.hpp"

// traccc plugin include (s)
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

namespace Acts::TracccConversion{

// https://github.com/acts-project/traccc/blob/436424f777b45c583754718c470f1c70b87ad11e/core/include/traccc/edm/cell.hpp#L27
// https://github.com/acts-project/traccc/blob/436424f777b45c583754718c470f1c70b87ad11e/core/include/traccc/edm/cell.hpp#L27
// https://github.com/acts-project/traccc/blob/436424f777b45c583754718c470f1c70b87ad11e/io/src/csv/read_cells.cpp#L165
// m_channelizer.channelize in digitalization algorithm.


auto convertData(std::vector<std::vector<std::tuple<double, int>>> grid,
                const traccc::geometry* geom,
                const std::map<std::uint64_t, detray::geometry::barcode>* barcode_map){
    
    traccc::cell_collection_types::host cells{};
    traccc::cell_module_collection_types::host modules{};

    return std::tie(cells, modules);

}

// Code modified based off of traccc/io/src/csv/read_cells.cpp 
// (keep the conversion code while removing the requirement that the input is a file)
namespace {
    
/// Helper function which finds module from csv::cell in the geometry and
/// digitization config, and initializes the modules limits with the cell's
/// properties
traccc::cell_module get_module(const std::uint64_t geometry_id,
                               const traccc::geometry* geom,
                               const traccc::digitization_config* dconfig,
                               const std::uint64_t original_geometry_id) {

    traccc::cell_module result;
    result.surface_link = detray::geometry::barcode{geometry_id};

    // Find/set the 3D position of the detector module.
    if (geom != nullptr) {

        // Check if the module ID is known.
        if (!geom->contains(result.surface_link.value())) {
            throw std::runtime_error(
                "Could not find placement for geometry ID " +
                std::to_string(result.surface_link.value()));
        }

        // Set the value on the module description.
        result.placement = (*geom)[result.surface_link.value()];
    }

    // Find/set the digitization configuration of the detector module.
    if (dconfig != nullptr) {

        // Check if the module ID is known.
        const traccc::digitization_config::Iterator geo_it =
            dconfig->find(original_geometry_id);
        if (geo_it == dconfig->end()) {
            throw std::runtime_error(
                "Could not find digitization config for geometry ID " +
                std::to_string(original_geometry_id));
        }

        // Set the value on the module description.
        const auto& binning_data = geo_it->segmentation.binningData();
        assert(binning_data.size() >= 2);
        result.pixel = {binning_data[0].min, binning_data[1].min,
                        binning_data[0].step, binning_data[1].step};
    }

    return result;
}

auto convertCellsMap(
    const std::map<std::uint64_t, std::vector<traccc::cell>> cellsMap,
    const traccc::geometry* geom,
    const traccc::digitization_config* dconfig,
    const std::map<std::uint64_t,
    detray::geometry::barcode>* barcode_map) {

    traccc::cell_collection_types::host resCells{};
    traccc::cell_module_collection_types::host resModules{};

    // Fill the output containers with the ordered cells and modules.
    for (const auto& [original_geometry_id, cells] : cellsMap) {
        // Modify the geometry ID of the module if a barcode map is
        // provided.
        std::uint64_t geometry_id = original_geometry_id;
        if (barcode_map != nullptr) {
            const auto it = barcode_map->find(geometry_id);
            if (it != barcode_map->end()) {
                geometry_id = it->second.value();
            } else {
                throw std::runtime_error(
                    "Could not find barcode for geometry ID " +
                    std::to_string(geometry_id));
            }
        }

        // Add the module and its cells to the output.
        resModules.push_back(
            get_module(geometry_id, geom, dconfig, original_geometry_id));
        for (auto& cell : cells) {
            resCells.push_back(cell);
            // Set the module link.
            resCells.back().module_link = resModules.size() - 1;
        }
    }
    return std::tie(resCells, resModules);
}

};

template <typename track_container_t, typename traj_t, template <typename> class holder_t>
//Acts::Result<typename Acts::TrackContainer<track_container_t, traj_t, holder_t>::TrackProxy>
void convert(const traccc::track_state_container_types::host& data, Acts::TrackContainer<track_container_t, traj_t, holder_t>& trackContainer){

    //Acts::Result<typename Acts::TrackContainer<track_container_t, traj_t, holder_t>::TrackProxy> result{};

    for (std::size_t i = 0; i < data.size(); i++){
        
        auto c = data[i];
        auto fittingResult = c.header;
        auto trackStates = c.items;

        auto track = trackContainer.makeTrack();
        


    }

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

    template <typename track_container_t, typename traj_t, template <typename> class holder_t>
    auto run(traccc::cell_collection_types::host& cells,
        traccc::cell_module_collection_types::host& modules,
        Acts::TrackContainer<track_container_t, traj_t, holder_t>& trackContainer){
        auto res = run(cells, modules);
        convert(res, trackContainer);
        }

    private:

    auto run(traccc::cell_collection_types::host& cells,
        traccc::cell_module_collection_types::host& modules) const {

        typename clusterization_func_t::output_type measurements{&hostMemoryResource};
        typename spacepoint_formation_func_t::output_type spacepoints{&hostMemoryResource};
        typename seeding_func_t::output_type seeds{&hostMemoryResource};
        typename track_params_estimation_func_t::output_type params{&hostMemoryResource};
        typename finding_func_t::output_type trackCandidates{&hostMemoryResource};
        typename fitting_func_t::output_type trackStates{&hostMemoryResource};
        typename resolution_func_t::output_type resolvedTrackStates{&hostMemoryResource};

        measurements = clusterizationFunc(vecmem::get_data(cells), vecmem::get_data(modules));
        spacepoints = spacepointFormationFunc(vecmem::get_data(measurements), vecmem::get_data(modules));
        seeds = seedingFunc(spacepoints);

        const typename field_t::view_t fieldView(field);

        //static_assert(std::is_same<field_t, typename detray::bfield::const_field_t>::value, "Currently, traccc expects a constant field.");
        params = trackParameterEstimationFunc(spacepoints, seeds, fieldView.at(0,0,0));

        trackCandidates = findingFunc(detector, field, measurements, params);
        trackStates = fittingFunc(detector, field, trackCandidates);
        resolvedTrackStates = resolutionFunc(trackStates);
        return resolvedTrackStates;
    }


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

    public:

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


template <typename detector_t> 
std::pair<traccc::geometry, std::unique_ptr<std::map<std::uint64_t, detray::geometry::barcode>>>
getGeometry(const detector_t& detector) {

    traccc::geometry surface_transforms = traccc::io::alt_read_geometry(detector);

    // Construct a map from Acts surface identifiers to Detray barcodes.
    std::unique_ptr<std::map<std::uint64_t, detray::geometry::barcode>> barcode_map = std::make_unique<std::map<std::uint64_t, detray::geometry::barcode>>();
    for (const auto& surface : detector.surfaces()) {
        (*barcode_map)[surface.source] = surface.barcode();
    }

    // Return the created objects.
    return {surface_transforms, std::move(barcode_map)};
}

auto readDetector(vecmem::host_memory_resource& host_mr, const std::string& detectorFile, const std::string& materialFile = "", const std::string& gridFile = ""){
    using detector_type = detray::detector<detray::default_metadata,
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
    auto [det, names] = detray::io::read_detector<detector_type>(host_mr, cfg);
    return std::move(det);
}

// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/utils/algorithm.hpp
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/container.hpp#L154
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/details/host_container.hpp#L25
// https://github.com/acts-project/traccc/blob/f51a1f8c4031fae664c752695e177b8315fcf22b/core/include/traccc/edm/details/container_base.hpp#L42
// https://github.com/acts-project/traccc/blob/f51a1f8c4031fae664c752695e177b8315fcf22b/core/include/traccc/edm/details/container_element.hpp#L28

// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/track_state.hpp

};
