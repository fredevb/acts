// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Plugin include(s)
#include "Acts/Plugins/Traccc/Utils.hpp"

// Acts include(s)
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackProxy.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Definitions/Algebra.hpp"

// Detray include(s)
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/core/detector.hpp"

// Traccc include(s)
#include "traccc/edm/track_state.hpp"

// System include(s)
#include <memory>

namespace Acts::TracccPlugin {

template <typename algebra_t, typename metadata_t, typename container_t>
auto newParams(const detray::bound_track_parameters<algebra_t>& dparams, const detray::detector<metadata_t, container_t>& detector, const Acts::TrackingGeometry& trackingGeometry){
    Acts::ActsVector<6U> parameterVector = Utils::newVector<6U>(dparams.vector()[0]);
    typename Acts::BoundTrackParameters::CovarianceMatrix cov = Utils::newSqaureMatrix<6U>(dparams.covariance());
    Acts::ParticleHypothesis particleHypothesis = Acts::ParticleHypothesis::pion();

    auto geoID = Acts::GeometryIdentifier(detector.surface(dparams.surface_link()).source);

    auto surface = trackingGeometry.findSurface(geoID);

    Acts::BoundTrackParameters params( //Here really
        surface->getSharedPtr(),
        parameterVector,
        std::make_optional(std::move(cov)),
        particleHypothesis
    );

    return params;
}

template <typename algebra_t, typename track_container_t, typename trajectory_t, template <typename> class holder_t, typename metadata_t, typename container_t>
void copyFittingResult(const traccc::fitting_result<algebra_t>& source, Acts::TrackProxy<track_container_t, trajectory_t, holder_t, false>& destination, const detray::detector<metadata_t, container_t>& detector, const Acts::TrackingGeometry& trackingGeometry){
    const auto params = newParams(source.fit_params, detector, trackingGeometry);
    //track.tipIndex() = kalmanResult.lastMeasurementIndex;
    destination.parameters() = params.parameters();
    destination.covariance() = params.covariance().value();
    destination.setReferenceSurface(params.referenceSurface().getSharedPtr());
}

/// @note will not contain a source link
template <typename algebra_t, typename metadata_t, typename container_t, typename trajectory_t, std::size_t M>
void copyTrackState(const traccc::track_state<algebra_t>& source, Acts::TrackStateProxy<trajectory_t, M, false>& destination, const detray::detector<metadata_t, container_t>& detector, const Acts::TrackingGeometry& trackingGeometry) {
    auto geoID = Acts::GeometryIdentifier(detector.surface(source.surface_link()).source);
    auto surface = trackingGeometry.findSurface(geoID)->getSharedPtr();
    destination.setReferenceSurface(surface);

    using Parameters = typename Acts::TrackStateProxy<trajectory_t, M, false>::Parameters;
    using Covariance = typename Acts::TrackStateProxy<trajectory_t, M, false>::Covariance;

    destination.predicted() = Parameters(Utils::newVector<6U>(source.predicted().vector()[0]).data());
    destination.predictedCovariance() = Covariance(Utils::newSqaureMatrix<6U>(source.predicted().covariance()).data());

    destination.smoothed() = Parameters(Utils::newVector<6U>(source.smoothed().vector()[0]).data());
    destination.smoothedCovariance() = Covariance(Utils::newSqaureMatrix<6U>(source.smoothed().covariance()).data());

    destination.filtered() = Parameters(Utils::newVector<6U>(source.filtered().vector()[0]).data());
    destination.filteredCovariance() = Covariance(Utils::newSqaureMatrix<6U>(source.filtered().covariance()).data());

    destination.jacobian() = Covariance(Utils::newSqaureMatrix<6U>(source.jacobian()).data());

    destination.chi2() = source.smoothed_chi2();

    auto typeFlags = destination.typeFlags();
    typeFlags.set(TrackStateFlag::ParameterFlag);
    if (surface->surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
    }
    if (source.is_hole){
        typeFlags.set(TrackStateFlag::HoleFlag);
    }
    typeFlags.set(TrackStateFlag::MeasurementFlag);
}

template <typename algebra_t>
const traccc::fitting_result<algebra_t>& getFittingResult(const traccc::container_element<const traccc::fitting_result<algebra_t> &, const vecmem::vector<traccc::track_state<algebra_t>> &>& e){
    return e.header;
}

template <typename algebra_t>
const vecmem::vector<traccc::track_state<algebra_t>>& getTrackStates(const traccc::container_element<const traccc::fitting_result<algebra_t> &, const vecmem::vector<traccc::track_state<algebra_t>> &>& e){
    return e.items;
}

template <typename algebra_t, typename track_container_t, typename trajectory_t, template <typename> class holder_t, typename metadata_t, typename container_t>
void makeTrack(
    const traccc::container_element<const traccc::fitting_result<algebra_t> &, const vecmem::vector<traccc::track_state<algebra_t>> &>& tracccContainerElement, 
    Acts::TrackContainer<track_container_t, trajectory_t, holder_t>& trackContainer, 
    const detray::detector<metadata_t, container_t>& detector, 
    const Acts::TrackingGeometry& trackingGeometry, 
    const std::vector<Acts::Measurement<Acts::BoundIndices, 2>>& measurements,
    const std::vector<Acts::SourceLink>& sourceLinks) {

    auto fittingResult = getFittingResult(tracccContainerElement);
    auto trackStates = getTrackStates(tracccContainerElement);

    auto track = trackContainer.makeTrack(); //CHeck the chi2 is not filled. Check where tracks get lost with printing

    copyFittingResult(fittingResult, track, detector, trackingGeometry);

    std::cout << trackContainer.size() << std::endl;


    // set the track states
    for (const auto& tstate : trackStates){
        auto astate = track.appendTrackState();

        copyTrackState(tstate, astate, detector, trackingGeometry);

        auto idx = tstate.get_measurement().measurement_id;

        auto& measurement = measurements[idx];
        astate.setCalibrated(measurement);

        auto sourceLink = sourceLinks[idx];
        astate.setUncalibratedSourceLink(sourceLink);
    }
    std::cout << trackContainer.size() << std::endl;

}

template <typename traccc_container_t, typename track_container_t, typename trajectory_t, template <typename> class holder_t, typename metadata_t, typename container_t>
void makeTracks(
    const traccc_container_t& data, 
    Acts::TrackContainer<track_container_t, trajectory_t, holder_t>& trackContainer, 
    const detray::detector<metadata_t, container_t>& detector, 
    const Acts::TrackingGeometry& trackingGeometry, 
    const std::vector<Acts::Measurement<Acts::BoundIndices, 2>>& measurements, 
    const std::vector<Acts::SourceLink>& sourceLinks) {
    for (std::size_t i = 0; i < data.size(); i++) {
        makeTrack(data[i], trackContainer, detector, trackingGeometry, measurements, sourceLinks);
    }
    std::cout << "final size: " << trackContainer.size() << std::endl;
}

}