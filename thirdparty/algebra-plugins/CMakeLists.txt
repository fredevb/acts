# This file is part of the Acts project.
#
# Copyright (C) 2021 CERN for the benefit of the Acts project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

# CMake include(s).
cmake_minimum_required(VERSION 3.14)
include(FetchContent)

# Tell the user what's happening.
message( STATUS "Building algebra-plugins as part of the Acts project" )

FetchContent_Declare( AlgebraPlugins ${ACTS_ALGEBRAPLUGINS_SOURCE} )

set( ALGEBRA_PLUGINS_BUILD_TESTING FALSE CACHE BOOL
   "Turn off the build of the Algebra Plugins unit tests" )
set( ALGEBRA_PLUGINS_INCLUDE_EIGEN TRUE CACHE BOOL
   "Turn on the build of algebra::eigen" )
set( ALGEBRA_PLUGINS_INCLUDE_VC TRUE CACHE BOOL
   "Turn on the build of algebra::vc_array" )
set( ALGEBRA_PLUGINS_INCLUDE_VECMEM TRUE CACHE BOOL
   "Turn on the build of algebra::vecmem_array" )

set( ALGEBRA_PLUGINS_SETUP_VC TRUE CACHE BOOL
   "Have Algebra Plugins set up Vc for itself" )
set( ALGEBRA_PLUGINS_USE_SYSTEM_VC FALSE CACHE BOOL
   "Have Algebra Plugins build Vc itself" )
set( ALGEBRA_PLUGINS_SETUP_EIGEN3 FALSE CACHE BOOL
   "Do not set up Eigen3 in Algebra Plugins" )
set( ALGEBRA_PLUGINS_SETUP_VECMEM FALSE CACHE BOOL
   "Do not set up VecMem in Algebra Plugins" )
set( ALGEBRA_PLUGINS_SETUP_GOOGLETEST FALSE CACHE BOOL
   "Do not set up GoogleTest in Algebra Plugins" )
set( ALGEBRA_PLUGINS_SETUP_BENCHMARK FALSE CACHE BOOL
   "Do not set up GoogleTest in Algebra Plugins" )

# Now set up its build.
FetchContent_MakeAvailable( AlgebraPlugins )