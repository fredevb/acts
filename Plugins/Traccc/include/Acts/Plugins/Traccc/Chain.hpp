// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Plugin include(s)
#include "Acts/Plugins/Traccc/CellConversion.hpp"
#include "Acts/Plugins/Traccc/TrackConversion.hpp"
#include "Acts/Plugins/Covfie/CovfieConversion.hpp"

// Acts include(s)
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Utilities/BinUtility.hpp"

// Traccc include(s)
#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clusterization_algorithm.hpp"
#include "traccc/clusterization/spacepoint_formation_algorithm.hpp"
#include "traccc/finding/finding_algorithm.hpp"
#include "traccc/fitting/fitting_algorithm.hpp"
#include "traccc/seeding/seeding_algorithm.hpp"
#include "traccc/seeding/track_params_estimation.hpp"
#include "traccc/finding/finding_config.hpp"
#include "traccc/edm/cell.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
#include "detray/io/frontend/detector_reader.hpp"
#include "detray/propagator/rk_stepper.hpp"

// VecMem include(s).
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s).
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <map>

// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/utils/algorithm.hpp
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/container.hpp#L154
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/details/host_container.hpp#L25
// https://github.com/acts-project/traccc/blob/f51a1f8c4031fae664c752695e177b8315fcf22b/core/include/traccc/edm/details/container_base.hpp#L42
// https://github.com/acts-project/traccc/blob/f51a1f8c4031fae664c752695e177b8315fcf22b/core/include/traccc/edm/details/container_element.hpp#L28

// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/track_state.hpp

// https://github.com/acts-project/traccc/blob/436424f777b45c583754718c470f1c70b87ad11e/core/include/traccc/edm/cell.hpp#L27
// https://github.com/acts-project/traccc/blob/436424f777b45c583754718c470f1c70b87ad11e/core/include/traccc/edm/cell.hpp#L27
// https://github.com/acts-project/traccc/blob/436424f777b45c583754718c470f1c70b87ad11e/io/src/csv/read_cells.cpp#L165
// m_channelizer.channelize in digitalization algorithm.


namespace{

auto readDetector(vecmem::memory_resource& mr, const std::string& detectorFilePath, const std::string& materialFilePath = "", const std::string& gridFilePath = ""){
    using detector_type = detray::detector<detray::default_metadata,
                                        detray::host_container_types>;
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
    auto [det, names] = detray::io::read_detector<detector_type>(mr, cfg);
    return std::move(det);
}

}

namespace Acts::TracccPlugin::Chain {
    
template <typename field_t,
        typename clusterization_func_t, 
        typename spacepoint_formation_func_t, 
        typename seeding_func_t, 
        typename track_params_estimation_func_t, 
        typename finding_func_t, 
        typename fitting_func_t,
        typename resolution_func_t>
class TracccChain{

    public:
    using detray_detector_type = detray::detector<detray::default_metadata,
                                        detray::host_container_types>;
    TracccChain(
        std::shared_ptr<const Acts::TrackingGeometry> trackingGeom,
        vecmem::memory_resource& mr,
        std::shared_ptr<const GeometryHierarchyMap<Acts::BinUtility>> geoSegs,
        const field_t&& f,
        clusterization_func_t& cFunc,
                spacepoint_formation_func_t& sfFunc,
                seeding_func_t& sFunc,
                track_params_estimation_func_t& tpeFunc,
                finding_func_t& findFunc,
                fitting_func_t& fitFunc,
                resolution_func_t& resFunc) 
            : 
            trackingGeometry(trackingGeom),
            memoryResource(mr),
            detector(readDetector(memoryResource, "/home/frederik/Desktop/CERN-TECH/input/detray_detector.json")), //replace with detector conversion fucnction when the work on it is complete
            segmentations(geoSegs),
            field(f),
            clusterizationFunc(cFunc),
            spacepointFormationFunc(sfFunc),
            seedingFunc(sFunc),
            trackParameterEstimationFunc(tpeFunc),
            findingFunc(findFunc),
            fittingFunc(fitFunc),
            resolutionFunc(resFunc),
            cellDataConverter(CellConversion::CellDataConverter(memoryResource, detector, *segmentations)){}

    template <typename track_container_t, typename traj_t, template <typename> class holder_t, typename CellCollection>
    void operator()(const std::map<Acts::GeometryIdentifier, CellCollection> map,
        Acts::TrackContainer<track_container_t, traj_t, holder_t>& out) const{
        auto [cells, modules] = cellDataConverter(map);
        auto res = run(cells, modules);
        copyTrackContainer(res, out, detector, trackingGeometry);
        }

    private:

    auto run(const traccc::cell_collection_types::host& cells,
        const traccc::cell_module_collection_types::host& modules) const {

        typename clusterization_func_t::output_type measurements{&memoryResource};
        typename spacepoint_formation_func_t::output_type spacepoints{&memoryResource};
        typename seeding_func_t::output_type seeds{&memoryResource};
        typename track_params_estimation_func_t::output_type params{&memoryResource};
        typename finding_func_t::output_type trackCandidates{&memoryResource};
        typename fitting_func_t::output_type trackStates{&memoryResource};
        typename resolution_func_t::output_type resolvedTrackStates{&memoryResource};

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


    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    vecmem::memory_resource& memoryResource;
    detray_detector_type detector;
    std::shared_ptr<const Acts::GeometryHierarchyMap<Acts::BinUtility>> segmentations;
    const field_t field;
    clusterization_func_t clusterizationFunc;
    spacepoint_formation_func_t spacepointFormationFunc;
    seeding_func_t seedingFunc;
    track_params_estimation_func_t trackParameterEstimationFunc;
    finding_func_t findingFunc;
    fitting_func_t fittingFunc;
    resolution_func_t resolutionFunc;
    const CellConversion::CellDataConverter<detray::default_metadata, detray::host_container_types> cellDataConverter;
};

class TracccChainFactory{

    public:

    using detector_t = detray::detector<detray::default_metadata, detray::host_container_types>;

    using stepper_t = detray::rk_stepper<typename detray::bfield::const_field_t::view_t, typename detector_t::transform3, detray::constrained_step<>>;
    
    using navigator_t = detray::navigator<const detector_t>;
    
    using finding_algorithm_t = traccc::finding_algorithm<stepper_t, navigator_t>;
    
    using fitter_t = traccc::kalman_fitter<stepper_t, navigator_t>;
    using fitting_algorithm_t = traccc::fitting_algorithm<fitter_t>;

    using clusterization_func_t = traccc::host::clusterization_algorithm;
    using spacepoint_formation_func_t = traccc::host::spacepoint_formation_algorithm;
    using seeding_func_t = traccc::seeding_algorithm;
    using track_params_estimation_func_t = traccc::track_params_estimation;
    using resolution_func_t = traccc::greedy_ambiguity_resolution_algorithm;

    using field_t = Acts::CovfieConversion::constant_field_t;

    using chain_t = TracccChain<field_t, clusterization_func_t, spacepoint_formation_func_t, seeding_func_t, track_params_estimation_func_t, finding_algorithm_t, fitting_algorithm_t, resolution_func_t>;

    auto buildChainHost(
                std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
                vecmem::host_memory_resource& host_mr,
                std::shared_ptr<const GeometryHierarchyMap<Acts::BinUtility>> geoSegs,
                //std::shared_ptr<const Acts::MagneticFieldProvider> fieldProvider,
                std::shared_ptr<const Acts::ConstantBField> fieldProvider,
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

        const field_t field = Acts::CovfieConversion::covfieField(*fieldProvider);

        TracccChain tc(trackingGeometry, host_mr, geoSegs, std::move(field), ca, sf, sa, tp, findFunc, fitFunc, res);
        return tc;
    }
};

}