// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Covfie/CovfieConversion.hpp"

namespace Acts::CovfieConversion{

/// @brief Get the value of the interpolated field at a specific position (unchecked).
/// @tparam cache_t 
/// @param magneticField 
/// @param cache 
/// @param position 
/// @return 
template <typename cache_t>
auto fieldLookup(const Acts::InterpolatedMagneticField& magneticField, [[maybe_unused]] cache_t& cache, const Acts::Vector3& position){
    return magneticField.getFieldUnchecked(position);;
}

/// @brief Get the value of the field at a specific position of a general magnetic field.
/// @tparam cache_t 
/// @param magneticField 
/// @param cache 
/// @param position 
/// @return 
template <typename cache_t>
auto fieldLookup(const Acts::MagneticFieldProvider& magneticField, cache_t& cache, const Acts::Vector3& position){
    auto lookupResult = magneticField.getField(position, cache);
    if(!lookupResult.ok()) {
        throw std::runtime_error{"Field lookup failure"};
    }
    return *lookupResult;
}


/// @brief Creates a strided covfie field that stores the values of the magnetic field in the volume given by min and max using a fixed sample spacing (determined by nBins).
/// @param magneticField The acts magnetic field.
/// @param cache The acts cache.
/// @param nBins 3D array of containing the number of bins for each axis.
/// @param min (min_x, min_y, min_z)
/// @param max (max_x, max_y, max_z)
/// @return A strided covfie field.
template <typename magnetic_field_t, typename cache_t, typename point3_1_t, typename point3_2_t, typename point3_3_t>
auto newBuilder(const magnetic_field_t& magneticField, cache_t& cache, const point3_1_t& nBins, const point3_2_t& min, const point3_3_t& max){

    using field_t = covfie::field<builder_backend_t>;

    field_t field(covfie::make_parameter_pack(
        field_t::backend_t::configuration_t{nBins[0], nBins[1], nBins[2]}
    ));

    field_t::view_t view(field);

    std::array<double, 3> sampleSpacing = {
        (max[0]-min[0])/(nBins[0]-1),
        (max[1]-min[1])/(nBins[1]-1),
        (max[2]-min[2])/(nBins[2]-1)
    };

    for (std::size_t x = 0; x < nBins[0]; x++) {
        for (std::size_t y = 0; y < nBins[1]; y++) {
            for (std::size_t z = 0; z < nBins[2]; z++) {

                auto position = Acts::Vector3{
                    x*sampleSpacing[0]+min[0],
                    y*sampleSpacing[1]+min[1],
                    z*sampleSpacing[2]+min[2]
                };

                auto result = fieldLookup(magneticField, cache, position);

                field_t::view_t::output_t &p = view.at(x, y, z);
                p[0] = static_cast<float>(result[0]);
                p[1] = static_cast<float>(result[1]);
                p[2] = static_cast<float>(result[2]);
            }
        }
    }
    return field;
}

/// @brief Generated the affine covfie configuration (scaling and rotation) given the size of the field (min and max)
/// @param nBins 3D array of containing the number of bins for each axis.
/// @param min (min_x, min_y, min_z)
/// @param max (max_x, max_y, max_z)
/// @return The affine field configuration.
template <typename backend_t, typename point3_1_t, typename point3_2_t, typename point3_3_t>
auto affineConfiguration(const point3_1_t& nBins, const point3_2_t& min, const point3_3_t& max){
    auto scaling = covfie::algebra::affine<3>::scaling(
        static_cast<float>((nBins[0] - 1) / (max[0] - min[0])),
        static_cast<float>((nBins[1] - 1) / (max[1] - min[1])),
        static_cast<float>((nBins[2] - 1) / (max[2] - min[2]))
    );

    auto translation = covfie::algebra::affine<3>::translation(
        static_cast<float>(-min[0]),
        static_cast<float>(-min[1]),
        static_cast<float>(-min[2])
    );

    return typename backend_t::configuration_t(scaling * translation);
}

/// @brief Creates a covfie field from a generic magnetic field.
/// @param magneticField The generic magnetic field.
/// @param cache The cache.
/// @param nBins 3D array of containing the number of bins for each axis.
/// @param min (min_x, min_y, min_z)
/// @param max (max_x, max_y, max_z)
/// @return An affine linear strided covfie field.
template <typename magnetic_field_t, typename cache_t, typename point3_1_t, typename point3_2_t, typename point3_3_t>
affine_linear_strided_field_t covfieFieldLinear(const magnetic_field_t& magneticField, cache_t& cache, const point3_1_t& nBins, const point3_2_t& min, const point3_3_t& max){
    auto builder = newBuilder(magneticField, cache, nBins, min, max);

    affine_linear_strided_field_t field(
        covfie::make_parameter_pack(
        affineConfiguration<affine_linear_strided_field_t::backend_t>(nBins, min, max),
        affine_linear_strided_field_t::backend_t::backend_t::configuration_t{},
        builder.backend()
    ));

    return field;
}

/// @brief Creates a covfie field from a magnetic field provider by sampling it.
/// @param magneticField The acts magnetic field provider.
/// @param cache The acts cache.
/// @param nBins 3D array of containing the number of bins for each axis.
/// @param min (min_x, min_y, min_z)
/// @param max (max_x, max_y, max_z)
/// @return An affine linear strided covfie field.
affine_linear_strided_field_t covfieField(const Acts::MagneticFieldProvider& magneticField, Acts::MagneticFieldProvider::Cache& cache, const std::vector<std::size_t>& nBins, const std::vector<double>& min, const std::vector<double>& max){
    return covfieFieldLinear(magneticField, cache, nBins, min, max);
}

/// @brief Creates a covfie field from an interpolated magnetic field.
/// @param magneticField The acts interpolated magnetic field.
/// @return An affine linear strided covfie field.
affine_linear_strided_field_t covfieField(const Acts::InterpolatedMagneticField& magneticField){
    Acts::MagneticFieldContext ctx;
    auto cache = magneticField.makeCache(ctx);
    return covfieFieldLinear(magneticField, cache, magneticField.getNBins(), magneticField.getMin(), magneticField.getMax());
}

/// @brief Creates a covfie field from a constant B field.
/// @param magneticField The acts constant magnetic field.
/// @return A constant covfie field.
constant_field_t covfieField(const Acts::ConstantBField& magneticField){
    auto B = magneticField.getField();
    constant_field_t field(covfie::make_parameter_pack(constant_field_t::backend_t::configuration_t{static_cast<float>(B[0]), static_cast<float>(B[1]), static_cast<float>(B[2])}));
    return field;
}

};
