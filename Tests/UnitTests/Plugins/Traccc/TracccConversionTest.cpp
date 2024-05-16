// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Boost.Test include(s).
#include <boost/test/unit_test.hpp>

// Acts Traccc plugin include(s).
#include "Acts/Plugins/Traccc/CellConversion.hpp"
#include "Acts/Plugins/Traccc/TrackConversion.hpp"

// Detray include(s).
#include "detray/core/detector.hpp"
#include "detray/detectors/bfield.hpp"
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

using namespace Acts::TracccPlugin;

/*namespace{

using detector_t = detray::detector<detray::default_metadata, detray::host_container_types>;

Acts::GeometryHierarchyMap<Acts::BinUtility> createSegmentations(){
    return Acts::GeometryHierarchyMap<Acts::BinUtility>();
}

struct Cell{
    int bin0, bin1;
    double time, activation;
};

struct Getter : CellDataGetter<Cell>{
    double getTime(const Cell& cell){
        return cell.time;
    }
    double getActivation(const Cell& cell){
        return cell.activation;
    }
    int getRow(const Cell& cell){
        cell.bin0;
    }
    int getColumn(const Cell& cell){
        cell.bin1;
    }
};

detector_t createDetrayDetector(){
    auto det = detector_t();
    det.add_detector(detray::bfield_detector::create());
    return det;
}

std::shared_ptr<const Acts::TrackingGeometry> createTrackingGeometery(){
    return std::make_shared<const Acts::TrackingGeometry>();
}

}

BOOST_AUTO_TEST_CASE(Traccc_Conversion_Cell) {
    auto det = createDetrayDetector();
    auto segs = createSegmentations();
    Getter getter;
    CellDataConverter converter(det, segs, getter);

    auto tg = createTrackingGeometery();
    auto inputData = createActsCellData(tg);
    vecmem::host_memory_resource mr;
    auto outputData = createTracccCellData(tg, mr);
    auto result = converter(inputData, mr);
    BOOST_TEST(CellDataEqual(result, outputData));
}*/

/*BOOST_AUTO_TEST_CASE(Traccc_Conversion_Track) {
    auto det = createDetrayDetector();
    auto tg = createTrackingGeometery();
    auto inputData = createTracccTracks();
    auto outputData = createActsTracks();
    auto result;
    copyTrackContainer(inputData, result, det, tg);
    BOOST_TEST(TrackDataEqual(outTrackContainer, outputData));
}*/

// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/utils/algorithm.hpp
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/container.hpp#L154
// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/details/host_container.hpp#L25
// https://github.com/acts-project/traccc/blob/f51a1f8c4031fae664c752695e177b8315fcf22b/core/include/traccc/edm/details/container_base.hpp#L42
// https://github.com/acts-project/traccc/blob/f51a1f8c4031fae664c752695e177b8315fcf22b/core/include/traccc/edm/details/container_element.hpp#L28

// https://github.com/acts-project/traccc/blob/710b62cac11c0dd9e139bd82d7cbafa4bc863b6c/core/include/traccc/edm/track_state.hpp