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
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/SourceLink.hpp"

// Detray include(s)
#include "detray/tracks/bound_track_parameters.hpp"
#include "detray/core/detector.hpp"

// Traccc include(s)
#include "traccc/edm/track_state.hpp"
#include "traccc/edm/measurement.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"

// System include(s)
#include <memory>

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
    
    auto cov = Eigen::DiagonalMatrix<Acts::ActsScalar, static_cast<int>(dim)>(Utils::newVector<dim>(m.variance)).toDenseMatrix();

    return Acts::Measurement<Acts::BoundIndices, dim>(std::move(sl), indices, params, cov);
}


template <std::size_t max_dim = 4UL, typename source_link_t>
Acts::BoundVariantMeasurement boundVariantMeasurement(const traccc::measurement& m, const source_link_t& gsl){
    if constexpr (max_dim == 0UL){
        std::string errorMsg = "Invalid/mismatching measurement dimension: " +
                    std::to_string(m.meas_dim);
        throw std::runtime_error(errorMsg.c_str());
    }
    else{
        if (m.meas_dim == max_dim){
            return measurement<max_dim>(m, gsl);
        }
        return boundVariantMeasurement<max_dim-1>(m, gsl);
    }
}

}
