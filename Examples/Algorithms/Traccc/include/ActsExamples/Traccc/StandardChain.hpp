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
#include "Acts/Plugins/Traccc/MeasurementConversion.hpp"
#include "Acts/Plugins/Covfie/CovfieConversion.hpp"

// Acts include(s)
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "ActsExamples/EventData/Track.hpp"

// Acts examples include(s)
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"
#include "ActsExamples/Traccc/Converter.hpp"

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

namespace ActsExamples {

template <typename ca_t, typename sfa_t, typename sa_t, typename tpea_t, typename fnda_t, typename fta_t, typename ra_t>
struct StandardChain{
    using clusterization_algorithm_t = ca_t;
    using spacepoint_formation_algorithm_t = sfa_t;
    using seeding_algorithm_t = sa_t;
    using track_params_estimation_algorithm_t = tpea_t;
    using finding_algorithm_t = fnda_t;
    using fitting_algorithm_t = fta_t;
    using resolution_algorithm_t = ra_t;

    const clusterization_algorithm_t clusterizationAlgorithm;
    const spacepoint_formation_algorithm_t spacepointFormationAlgorithm;
    const seeding_algorithm_t seedingAlgorithm;
    const track_params_estimation_algorithm_t trackParameterEstimationAlgorithm;
    const finding_algorithm_t findingAlgorithm;
    const fitting_algorithm_t fittingAlgorithm;
    const resolution_algorithm_t resolutionAlgorithm;
};


class StandardChainRunner {

    public:

    using detector_t = detray::detector<detray::default_metadata, detray::host_container_types>;

    StandardChainRunner(
        std::shared_ptr<const Acts::TrackingGeometry> tg,
        std::shared_ptr<const detector_t> det,
        const Acts::GeometryHierarchyMap<DigiComponentsConfig>& dccMap)
        :
        trackingGeometry(tg),
        detector(det),
        converter(trackingGeometry, detector, dccMap)
        {}

    template <typename field_t, typename chain_t>
    ConstTrackContainer run(const std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>> map, const field_t& field, const chain_t& chain) const{
        vecmem::host_memory_resource mr;
        auto [cells, modules] = converter.convertInput(map, &mr);
        std::cout << "Converted." << std::endl;

        typename chain_t::clusterization_algorithm_t::output_type measurements{&mr};
        typename chain_t::spacepoint_formation_algorithm_t::output_type spacepoints{&mr};
        typename chain_t::seeding_algorithm_t::output_type seeds{&mr};
        typename chain_t::track_params_estimation_algorithm_t::output_type params{&mr};
        typename chain_t::finding_algorithm_t::output_type trackCandidates{&mr};
        typename chain_t::fitting_algorithm_t::output_type trackStates{&mr};
        typename chain_t::resolution_algorithm_t::output_type resolvedTrackStates{&mr};

        measurements = chain.clusterizationAlgorithm(vecmem::get_data(cells), vecmem::get_data(modules)); //HERE, note gcda error in python script at the start. (different behaviour on full clean build?)
        std::cout << "Clustered." << std::endl;
        spacepoints = chain.spacepointFormationAlgorithm(vecmem::get_data(measurements), vecmem::get_data(modules));
        std::cout << "Spacepointed." << std::endl;
        seeds = chain.seedingAlgorithm(spacepoints);
        std::cout << "Seeded." << std::endl;

        const typename field_t::view_t fieldView(field);

        //static_assert(std::is_same<field_t, typename detray::bfield::const_field_t>::value, "Currently, traccc expects a constant field.");
        params = chain.trackParameterEstimationAlgorithm(spacepoints, seeds, fieldView.at(0.f,0.f,0.f));//traccc::vector3{0.f, 0.f, 1.f}); //fieldView.at(0.f,0.f,0.f
        std::cout << "Params." << std::endl;

        trackCandidates = chain.findingAlgorithm(*detector, field, measurements, params); //here really
        std::cout << "Found." << std::endl;
        trackStates = chain.fittingAlgorithm(*detector, field, trackCandidates);
        std::cout << "Fitted." << std::endl;
        resolvedTrackStates = chain.resolutionAlgorithm(trackStates);
        std::cout << "Resolved." << std::endl;

        return converter.convertOutput(measurements, resolvedTrackStates);
    }

    private:
    
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const detector_t> detector;
    ActsExamples::TracccConversion::Converter<detector_t> converter;
};

}