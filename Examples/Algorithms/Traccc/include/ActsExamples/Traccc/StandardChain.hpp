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
#include "Acts/EventData/TrackContainer.hpp"
#include "ActsExamples/EventData/Track.hpp"

// Acts examples include(s)
#include "ActsExamples/Digitization/DigitizationAlgorithm.hpp"

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

float getTime(const Cluster::Cell& /*cell*/){
    return 0.f;
}

float getActivation(const Cluster::Cell& cell){
    return static_cast<float>(cell.activation);
}

unsigned int getRow(const Cluster::Cell& cell){
    if (cell.bin[0] > UINT_MAX) {
        throw std::runtime_error("Overflow will occur when casting to unsigned int.");
    }
    return static_cast<unsigned int>(cell.bin[0]);
}

unsigned int getColumn(const Cluster::Cell& cell){
        if (cell.bin[0] > UINT_MAX) {
        throw std::runtime_error("Overflow will occur when casting to unsigned int.");
    }
    return static_cast<unsigned int>(cell.bin[1]);
}



Acts::BinUtility getSegmentation(const DigiComponentsConfig& dcc){
    return dcc.geometricDigiConfig.segmentation;
}

template <typename detector_t, typename allocator_t>
auto getIndexSourceLinkAndMeasurements(const detector_t& detector, const std::vector<traccc::measurement, allocator_t>& measurements){
  IndexSourceLinkContainer sourceLinks;
  MeasurementContainer measurementContainer;

  for (traccc::measurement& m : measurements) {
    Acts::GeometryIdentifier moduleGeoId(detector.surface(m.surface_link).source);
    Index measurementIdx = measurementContainer.size();
    IndexSourceLink sourceLink{moduleGeoId, measurementIdx};
    sourceLinks.insert(sourceLinks.end(), sourceLink);
    measurements.push_back(Acts::TracccPlugin::boundVariantMeasurement(m, sourceLink));
  }
  return std::make_tuple(std::move(sourceLinks), std::move(measurementContainer));
}

class ChainAdapter{
    public:

    using detector_t = detray::detector<detray::default_metadata, detray::host_container_types>;
    using field_t = Acts::CovfieConversion::constant_field_t;

    ChainAdapter(
        std::shared_ptr<const Acts::TrackingGeometry> tg,
        const Acts::GeometryHierarchyMap<DigiComponentsConfig>& dccMap,
        const Acts::ConstantBField& f)
        :
        trackingGeometry(tg),
        detector(readDetector(resource, "/home/frederik/Desktop/CERN-TECH/input/odd-detray_geometry_detray.json")),
        field(Acts::CovfieConversion::covfieField(f)),
        cellDataConverter(detector, Acts::TracccPlugin::tracccConfig(dccMap, getSegmentation))
        {}

    template <typename chain_t>
    auto runChain(const std::map<Acts::GeometryIdentifier, std::vector<Cluster::Cell>> map, const chain_t& chain, vecmem::host_memory_resource* mr) const{
        // Convert input data
        auto tcm = Acts::TracccPlugin::tracccCellsMap(map, getRow, getColumn, getActivation, getTime);
        auto [cells, modules] = cellDataConverter(mr, tcm);
        std::cout << "Input converted." << std::endl;

        // Run chain
        auto res = chain(cells, modules, detector, field, mr);

        auto trackContainer = std::make_shared<Acts::VectorTrackContainer>();
        auto trackStateContainer = std::make_shared<Acts::VectorMultiTrajectory>();
        TrackContainer tracks(trackContainer, trackStateContainer);
        Acts::TracccPlugin::copyTrackContainer(res, tracks, detector, *trackingGeometry);
        std::cout << "Output converted." << std::endl;

        ConstTrackContainer constTracks{
            std::make_shared<Acts::ConstVectorTrackContainer>(std::move(*trackContainer)),
            std::make_shared<Acts::ConstVectorMultiTrajectory>(std::move(*trackStateContainer))
            };
        return constTracks;
    }

    private:

    vecmem::host_memory_resource resource;

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    const traccc::digitization_config digitalizationConfig;
    const detector_t detector;
    const field_t field;
    const Acts::TracccPlugin::CellDataConverter cellDataConverter;
};

template <typename clusterization_func_t, 
          typename spacepoint_formation_func_t, 
          typename seeding_func_t, 
          typename track_params_estimation_func_t, 
          typename finding_func_t, 
          typename fitting_func_t,
          typename resolution_func_t>
class StandardChainHost{
    public:
    StandardChainHost(
        const clusterization_func_t&& cFunc,
        const spacepoint_formation_func_t&& sfFunc,
        const seeding_func_t&& sFunc,
        const track_params_estimation_func_t&& tpeFunc,
        const finding_func_t&& findFunc,
        const fitting_func_t&& fitFunc,
        const resolution_func_t&& resFunc
        ): 
        clusterizationFunc(cFunc),
        spacepointFormationFunc(sfFunc),
        seedingFunc(sFunc),
        trackParameterEstimationFunc(tpeFunc),
        findingFunc(findFunc),
        fittingFunc(fitFunc),
        resolutionFunc(resFunc){}

    template <typename detector_t, typename field_t>
    auto operator()(
        const traccc::cell_collection_types::host& cells,
        const traccc::cell_module_collection_types::host& modules,
        const detector_t& detector,
        const field_t& field,
        vecmem::host_memory_resource* mr) 
        const {

        typename clusterization_func_t::output_type measurements{mr};
        typename spacepoint_formation_func_t::output_type spacepoints{mr};
        typename seeding_func_t::output_type seeds{mr};
        typename track_params_estimation_func_t::output_type params{mr};
        typename finding_func_t::output_type trackCandidates{mr};
        typename fitting_func_t::output_type trackStates{mr};
        typename resolution_func_t::output_type resolvedTrackStates{mr};

        measurements = clusterizationFunc(vecmem::get_data(cells), vecmem::get_data(modules)); //HERE, note gcda error in python script at the start. (different behaviour on full clean build?)
        std::cout << "Clustered." << std::endl;
        spacepoints = spacepointFormationFunc(vecmem::get_data(measurements), vecmem::get_data(modules));
        std::cout << "Spacepointed." << std::endl;
        seeds = seedingFunc(spacepoints);
        std::cout << "Seeded." << std::endl;

        const typename field_t::view_t fieldView(field);

        //static_assert(std::is_same<field_t, typename detray::bfield::const_field_t>::value, "Currently, traccc expects a constant field.");
        params = trackParameterEstimationFunc(spacepoints, seeds, fieldView.at(0.f,0.f,0.f));//traccc::vector3{0.f, 0.f, 1.f}); //fieldView.at(0.f,0.f,0.f
        std::cout << "Params." << std::endl;

        trackCandidates = findingFunc(detector, field, measurements, params); //here really
        std::cout << "Found." << std::endl;
        trackStates = fittingFunc(detector, field, trackCandidates);
        std::cout << "Fitted." << std::endl;
        resolvedTrackStates = resolutionFunc(trackStates);
        std::cout << "Resolved." << std::endl;
        return resolvedTrackStates;
    }

    const clusterization_func_t clusterizationFunc;
    const spacepoint_formation_func_t spacepointFormationFunc;
    const seeding_func_t seedingFunc;
    const track_params_estimation_func_t trackParameterEstimationFunc;
    const finding_func_t findingFunc;
    const fitting_func_t fittingFunc;
    const resolution_func_t resolutionFunc;
};

}