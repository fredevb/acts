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




#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"

namespace Acts::TracccPlugin {




Acts::BoundIndices boundIndex(const traccc::bound_indices tracccBoundIndex){
    switch (tracccBoundIndex)
    {
    case traccc::bound_indices::e_bound_loc0:
        return Acts::BoundIndices::eBoundLoc0;
    case traccc::bound_indices::e_bound_loc1:
        return Acts::BoundIndices::eBoundLoc1;
    case traccc::bound_indices::e_bound_phi:
        return Acts::BoundIndices::eBoundPhi;
    case traccc::bound_indices::e_bound_theta:
        return Acts::BoundIndices::eBoundTheta;
    case traccc::bound_indices::e_bound_qoverp:
        return Acts::BoundIndices::eBoundQOverP;
    case traccc::bound_indices::e_bound_time:
        return Acts::BoundIndices::eBoundTime;
    case traccc::bound_indices::e_bound_size:
        return Acts::BoundIndices::eBoundSize;
    default:
        throw std::runtime_error("Could not convert traccc bound index");
    }
}

template <std::size_t dim, typename source_link_t>
Acts::Measurement<Acts::BoundIndices, dim> measurement(const traccc::measurement& m, const source_link_t& gsl){
    Acts::SourceLink sl{gsl};
    auto params = Utils::newVector<dim>(m.local);

    std::array<Acts::BoundIndices, dim> indices;
    for (unsigned int i = 0; i < dim; i++){
        indices[i] = boundIndex(traccc::bound_indices(m.subs.get_indices()[i]));
    }
    
    Eigen::DiagonalMatrix<Acts::ActsScalar, dim, dim> cov(Utils::newVector<dim>(m.variance));

    return Acts::Measurement<Acts::BoundIndices, dim>(std::move(sl), indices, params, cov);
}

template <std::size_t max_dim = 4UL, typename source_link_t>
Acts::BoundVariantMeasurement boundVariantMeasurement(const traccc::measurement& m, const source_link_t& gsl){
    if (max_dim == 0UL){
        std::string errorMsg = "Invalid/mismatching measurement dimension: " +
                    std::to_string(m.meas_dim);
        throw std::runtime_error(errorMsg.c_str());
    }
    if (m.meas_dim == max_dim){
        return measurement<max_dim>(m, gsl);
    }
    return boundVariantMeasurement<max_dim-1>(m, gsl);
}


template <typename algebra_t, typename metadata_t, typename container_t>
auto newParams(const detray::bound_track_parameters<algebra_t>& dparams, const detray::detector<metadata_t, container_t>& detector, const Acts::TrackingGeometry& trackingGeometry){
    Acts::ActsVector<6U> parameterVector = Utils::newVector<6U>(dparams.vector());
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

template <typename algebra_t, typename metadata_t, typename container_t, typename trajectory_t, std::size_t M>
void copyTrackState(const traccc::track_state<algebra_t>& source, Acts::TrackStateProxy<trajectory_t, M, false>& destination, const detray::detector<metadata_t, container_t>& detector, const Acts::TrackingGeometry& trackingGeometry) {
    auto geoID = Acts::GeometryIdentifier(detector.surface(source.surface_link()).source);
    auto surface = trackingGeometry.findSurface(geoID)->getSharedPtr();
    destination.setReferenceSurface(surface);

    using Covariance = typename Acts::TrackStateProxy<trajectory_t, M, false>::Covariance;

    destination.predicted() = Utils::newVector<6U>(source.predicted().vector());
    destination.predictedCovariance() = Covariance(Utils::newSqaureMatrix<6U>(source.predicted().covariance()).data());

    destination.smoothed() = Utils::newVector<6U>(source.smoothed().vector());
    destination.smoothedCovariance() = Covariance(Utils::newSqaureMatrix<6U>(source.smoothed().covariance()).data());

    destination.filtered() = Utils::newVector<6U>(source.filtered().vector());
    destination.filteredCovariance() = Covariance(Utils::newSqaureMatrix<6U>(source.filtered().covariance()).data());

    destination.jacobian() = Covariance(Utils::newSqaureMatrix<6U>(source.jacobian()).data());

    //destination.setUncalibratedSourceLink(Acts::SourceLink(source));

    auto typeFlags = destination.typeFlags();
    typeFlags.set(TrackStateFlag::ParameterFlag);
    if (surface->surfaceMaterial() != nullptr) {
        typeFlags.set(TrackStateFlag::MaterialFlag);
    }
    if (source.is_hole){
        typeFlags.set(TrackStateFlag::HoleFlag);
    }
    typeFlags.set(TrackStateFlag::MeasurementFlag);

    //destination.jacobian() = source.jacobian;
    //destination.pathLength() = source.is_hole
    //destination.chi2 = source.filtered_chi2 or source.smoothed_chi2?
    //destination.parameters()
    //destination.
    //destination.typeFlags().set(Acts::TrackStateFlag::MeasurementFlag);
}

template <typename algebra_t, typename metadata_t, typename container_t, typename track_container_t, typename trajectory_t, template <typename> class holder_t>
void copyTrackStates(const vecmem::vector<traccc::track_state<algebra_t>>& source, Acts::TrackProxy<track_container_t, trajectory_t, holder_t, false>& destination, const detray::detector<metadata_t, container_t>& detector, const Acts::TrackingGeometry& trackingGeometry){
    for (const auto& tstate : source){
        auto astate = destination.appendTrackState();
        copyTrackState(tstate, astate, detector, trackingGeometry);
    }
}

template <typename traccc_track_state_container_t, typename track_container_t, typename traj_t, template <typename> class holder_t, typename metadata_t, typename container_t>
void copyTrackContainer(const traccc_track_state_container_t& data, Acts::TrackContainer<track_container_t, traj_t, holder_t>& trackContainer, const detray::detector<metadata_t, container_t>& detector, const Acts::TrackingGeometry& trackingGeometry) {

    for (std::size_t i = 0; i < data.size(); i++) {
        auto e = data[i];
        auto fittingResult = e.header;
        auto trackStates = e.items;

        auto track = trackContainer.makeTrack();

        copyFittingResult(fittingResult, track, detector, trackingGeometry);
        copyTrackStates(trackStates, track, detector, trackingGeometry);
    }
}

}