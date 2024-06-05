// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Plugin include(s)
#include "Acts/Plugins/Traccc/Detail/AlgebraConversion.hpp"

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
#include <variant>

namespace Acts::TracccPlugin {

/// @brief Converts a traccc bound index to an Acts bound index.
/// @param tracccBoundIndex the traccc bound index.
/// @return The Acts bound index.
inline Acts::BoundIndices boundIndex(const traccc::bound_indices tracccBoundIndex){
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

/// @brief Creates an Acts measurement from a traccc measurement.
/// @tparam dim the number of dimensions of the measurement (subspace size).
/// @param m the traccc measurement.
/// @param sl the Acts source link.
/// @return an Acts measurementwith data copied from the traccc measurement.
template <std::size_t dim>
inline Acts::Measurement<Acts::BoundIndices, dim> measurement(const traccc::measurement& m, const Acts::SourceLink sl){
    auto params = Detail::newVector<dim>(m.local);
    std::array<Acts::BoundIndices, dim> indices;
    for (unsigned int i = 0; i < dim; i++){
        indices[i] = boundIndex(traccc::bound_indices(m.subs.get_indices()[i]));
    }
    auto cov = Eigen::DiagonalMatrix<Acts::ActsScalar, static_cast<int>(dim)>(Detail::newVector<dim>(m.variance)).toDenseMatrix();
    return Acts::Measurement<Acts::BoundIndices, dim>(std::move(sl), indices, params, cov);
}

/// @brief Creates an Acts bound variant measurement from a traccc measurement.
/// Using recursion, the functions checks the dimension of the traccc measurement 
/// to determine the dimension of the Acts measurement that the bound variant measurement should hold.
/// The dimension must lie between [0; max_dim].
/// @tparam max_dim the largest possible dimension of any measurement type in the variant (default = 4)
/// @param m the traccc measurement.
/// @param sl the Acts source link.
/// @return an Acts bound variant measurement with data copied from the traccc measurement.
template <std::size_t max_dim = 4UL>
inline Acts::BoundVariantMeasurement boundVariantMeasurement(const traccc::measurement& m, const Acts::SourceLink sl){
    if constexpr (max_dim == 0UL){
        std::string errorMsg = "Invalid/mismatching measurement dimension: " +
                    std::to_string(m.meas_dim);
        throw std::runtime_error(errorMsg.c_str());
    }
    else{
        if (m.meas_dim == max_dim){
            return measurement<max_dim>(m, sl);
        }
        return boundVariantMeasurement<max_dim-1>(m, sl);
    }
}

template <std::size_t dim>
inline Acts::ActsVector<2> getLocal(const Acts::Measurement<Acts::BoundIndices, dim>& measurement){
    traccc::scalar loc0 = 0;
    traccc::scalar loc1 = 0;
    if constexpr (dim > Acts::BoundIndices::eBoundLoc0){
        loc0 = measurement.parameters()(Acts::BoundIndices::eBoundLoc0);
    }
    if constexpr (dim > Acts::BoundIndices::eBoundLoc1){
        loc1 = measurement.parameters()(Acts::BoundIndices::eBoundLoc1);
    }
    return Acts::ActsVector<2>(loc0, loc1);
}

inline Acts::ActsVector<2> getLocal(const Acts::BoundVariantMeasurement& measurement){
    return std::visit([](auto& m) { return getLocal(m); }, measurement);
}

template <std::size_t dim>
inline Acts::ActsVector<2> getVariance(const Acts::Measurement<Acts::BoundIndices, dim>& measurement){
    traccc::scalar var0 = 0;
    traccc::scalar var1 = 0;
    if constexpr (dim >= Acts::BoundIndices::eBoundLoc0){
        var0 = measurement.covariance()(Acts::BoundIndices::eBoundLoc0, Acts::BoundIndices::eBoundLoc0);
    }
    if constexpr (dim > Acts::BoundIndices::eBoundLoc1){
        var1 = measurement.covariance()(Acts::BoundIndices::eBoundLoc1, Acts::BoundIndices::eBoundLoc1);
    }
    return Acts::ActsVector<2>(var0, var1);
}

inline Acts::ActsVector<2> getVariance(const Acts::BoundVariantMeasurement& measurement){
    return std::visit([](auto& m) { return getVariance(m); }, measurement);
}

/// @brief Converts a geometry ID -> traccc cells map to traccc cells and modules.
/// @param mr The memory resource to use.
/// @param cellsMap A map from Acts geometry ID (value) to traccc cell.
/// @param geom The traccc geometry 
/// @param dconfig The traccc digitization configuration.
/// @param barcode_map A map from Acts geometry IDs (value) to detray barcodes.
/// @return A tuple containing the traccc cells (first item) and traccc modules (second item).
template <typename measurement_t>
inline auto createMeasurementsAndModules(
    vecmem::memory_resource* mr,
    std::map<std::uint64_t, measurement_t>& measurementMap,
    const traccc::geometry* geom,
    const traccc::digitization_config* dconfig,
    const std::map<std::uint64_t, detray::geometry::barcode>* barcodeMap = nullptr) {

    traccc::measurement_collection_types::host tracccMeasurements(mr);
    traccc::cell_module_collection_types::host tracccModules(mr);

    //detray::barcode::value_t;

    for (const auto& [original_geometry_id, measurements] : measurementMap) {

        auto geometry_id = barcodeMapTry(barcodeMap, original_geometry_id);

        auto mod = get_module(geometry_id, geom, dconfig, original_geometry_id);
        tracccModules.push_back(mod);
        for (auto& measurement : measurements) {
            traccc::measurement m{
                measurement.local(),
                measurement.variance(),
                geometry_id,
                measurement.measurementID(),
                tracccModules.size() - 1,
                std::numeric_limits<std::size_t>::max(),
                2u,
                {{0u, 1u}},
            };
            tracccMeasurements.push_back(std::move(m));
        }
    }
    return std::make_tuple(std::move(tracccMeasurements), std::move(tracccModules));
}

}
