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

namespace Acts::TracccPlugin::TrackConversion {

template <typename algebra_t, typename metadata_t, typename container_t>
Acts::BoundTrackParameters newParams(const detray::bound_track_parameters<algebra_t>& dparams, const detray::detector<metadata_t, container_t>& det, std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry){

    Acts::GeometryContext gctx;
    Acts::ActsVector<6U> parameterVector = Utils::newVector<6U>(dparams.vector());
    typename Acts::BoundTrackParameters::CovarianceMatrix cov = Utils::newSqaureMatrix<6U>(dparams.covariance());
    Acts::ParticleHypothesis particleHypothesis = Acts::ParticleHypothesis::pion();

    auto geoID = Acts::GeometryIdentifier(det.surface(dparams.surface_link()).source);

    auto surface = trackingGeometry->findSurface(geoID);

    Acts::BoundTrackParameters params(
        std::shared_ptr<const Acts::Surface>(surface),
        parameterVector,
        cov,
        particleHypothesis
    );

    return params;
}

template <typename track_container_t, typename trajectory_t, template <typename> class holder_t, typename metadata_t, typename container_t>
void copyFittingResult(const traccc::fitting_result<traccc::transform3>& source, Acts::TrackProxy<track_container_t, trajectory_t, holder_t, false>& destination, const detray::detector<metadata_t, container_t>& det, std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry){
    const auto params = newParams(source.fit_params, det, trackingGeometry);
    //track.tipIndex() = kalmanResult.lastMeasurementIndex;
    destination.parameters() = params.parameters();
    destination.covariance() = params.covariance().value();
    destination.setReferenceSurface(params.referenceSurface().getSharedPtr());
}

template <typename transform3_t, typename trajectory_t, std::size_t M>
void copyTrackState(const traccc::track_state<transform3_t>& /*source*/, Acts::TrackStateProxy<trajectory_t, M, false>& /*destination*/) {

}

template <typename track_container_t, typename trajectory_t, template <typename> class holder_t>
void copyTrackStates(const vecmem::vector<traccc::track_state<traccc::transform3>>& source, Acts::TrackProxy<track_container_t, trajectory_t, holder_t, false>& destination){
    for (const auto& tstate : source){
        auto astate = destination.appendTrackState();
        copyTrackState(tstate, astate);
    }
}

template <typename track_container_t, typename traj_t, template <typename> class holder_t, typename metadata_t, typename container_t>
void copyTrackContainer(const traccc::track_state_container_types::host& data, Acts::TrackContainer<track_container_t, traj_t, holder_t>& trackContainer, const detray::detector<metadata_t, container_t>& det, std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {

    for (std::size_t i = 0; i < data.size(); i++) {
        auto e = data[i];
        auto fittingResult = e.header;
        auto trackStates = e.items;

        auto track = trackContainer.makeTrack();

        copyFittingResult(fittingResult, track, det, trackingGeometry);
        copyTrackStates(trackStates, track);
    }
}
}