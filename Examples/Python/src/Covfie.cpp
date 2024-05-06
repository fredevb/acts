// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Covfie/CovfieConversion.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {

namespace{

template <typename field_t>
typename field_t::view_t newView(const field_t& field){
  typename field_t::view_t view(field);
  return view;
}

template <typename field_t>
void declareCovfieField(py::module &m, const std::string &fieldName) {
  using view_t = typename field_t::view_t;
  m.def("newView", static_cast<view_t (*)(const field_t&)>(&newView));
  py::class_<field_t, std::shared_ptr<field_t>>(m, fieldName.c_str());
  py::class_<view_t, std::shared_ptr<view_t>>(m, (fieldName + std::string("View")).c_str()).def("at", &view_t::template at<float, float, float>);
}

}
void addCovfie(Context& ctx) {
  auto main = ctx.get("main");
  auto m = main.def_submodule("covfie_conversion", "Submodule for covfie conversion");

  declareCovfieField<Acts::CovfieConversion::constant_field_t>(m, "CovfieConstantField");
  declareCovfieField<Acts::CovfieConversion::affine_linear_strided_field_t>(m, "CovfieAffineLinearStridedField");

  m.def("covfieField", py::overload_cast<const Acts::InterpolatedMagneticField&>(&Acts::CovfieConversion::covfieField));
  m.def("covfieField", py::overload_cast<const Acts::ConstantBField&>(&Acts::CovfieConversion::covfieField));
  m.def("covfieField", py::overload_cast<const Acts::MagneticFieldProvider&, Acts::MagneticFieldProvider::Cache&, const std::vector<std::size_t>&, const std::vector<double>&, const std::vector<double>&>(&Acts::CovfieConversion::covfieField));
}

}  // namespace Acts::Python
