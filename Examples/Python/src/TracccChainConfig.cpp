// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Traccc/Common/TracccChainConfig.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace Acts::Python {

void addTracccChainConfig(Context& ctx) {
  auto m = ctx.get("examples");

  using ConfigType = typename ActsExamples::Traccc::Common::TracccChainConfig;

  py::class_<ConfigType, std::shared_ptr<ConfigType>>(m, "TracccChainConfig")
  .def(py::init<>());

}

}  // namespace Acts::Python
