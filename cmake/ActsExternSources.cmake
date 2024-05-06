set( ACTS_ACTSVG_SOURCE
   "URL;https://github.com/acts-project/actsvg/archive/refs/tags/v0.4.40.tar.gz;URL_MD5;c8f3a7ac47db7f39cd4f215a5bbf148d" CACHE STRING "Source to take ACTSVG from")
mark_as_advanced( ACTS_ACTSVG_SOURCE )

set( ACTS_VECMEM_SOURCE
   "URL;https://github.com/acts-project/vecmem/archive/refs/tags/v1.4.0.tar.gz;URL_MD5;af5434e34ca9c084678c2c043441f174" CACHE STRING "Source to take VECMEM from")
mark_as_advanced( ACTS_VECMEM_SOURCE )

set( ACTS_ALGEBRAPLUGINS_SOURCE
   "URL;https://github.com/acts-project/algebra-plugins/archive/refs/tags/v0.22.0.tar.gz;URL_MD5;42bcaad8d19a2c773993a974816dfdf5" CACHE STRING "Source to take ALGEBRAPLUGINS from")
mark_as_advanced( ACTS_ALGEBRAPLUGINS_SOURCE )

set( ACTS_COVFIE_SOURCE
   "URL;https://github.com/acts-project/covfie/archive/refs/tags/v0.9.0.tar.gz;URL_MD5;b310712c6dd1acc8104c734086f40fc0" CACHE STRING "Source to take COVFIE from")
mark_as_advanced( ACTS_COVFIE_SOURCE )

set( ACTS_DETRAY_SOURCE
   "URL;https://github.com/acts-project/detray/archive/refs/tags/v0.63.0.tar.gz;URL_MD5;b03a737d09dc16ab9a88a34a7ac6b113" CACHE STRING "Source to take DETRAY from")
mark_as_advanced( ACTS_DETRAY_SOURCE )

set( ACTS_TRACCC_SOURCE
   "GIT_REPOSITORY;https://github.com/acts-project/traccc.git;GIT_TAG;origin/main" CACHE STRING "Source to take TRACCC from")
mark_as_advanced( ACTS_TRACCC_SOURCE )

set( ACTS_DFELIBS_SOURCE
   "URL;https://github.com/acts-project/dfelibs/archive/refs/tags/v20211029.tar.gz;URL_MD5;87fb09c5a11b98250f5e266e9cd501ea" CACHE STRING "Source to take dfelibs from")
mark_as_advanced( ACTS_DFELIBS_SOURCE )

set( ACTS_FRNN_SOURCE
   "GIT_REPOSITORY;https://github.com/lxxue/FRNN;GIT_TAG;3e370d8d9073d4e130363faf87d2370598b5fbf2" CACHE STRING "Source to take FRNN from")
mark_as_advanced( ACTS_FRNN_SOURCE )

set( ACTS_GEOMODEL_SOURCE
   "GIT_REPOSITORY;https://gitlab.cern.ch/GeoModelDev/GeoModel;GIT_TAG;4.6.0;PATCH_COMMAND;git am ${CMAKE_CURRENT_SOURCE_DIR}/0001-Add-option-to-skip-setting-up-json-completely.patch" CACHE STRING "Source to take GeoModel from")
mark_as_advanced( ACTS_GEOMODEL_SOURCE )

set( ACTS_NLOHMANNJSON_SOURCE
   "URL;https://github.com/nlohmann/json/archive/refs/tags/v3.10.5.tar.gz;URL_HASH;SHA1=8969f5ad1a422e01f040ff48dcae9c0e6ad0811d" CACHE STRING "Source to take nlohmann_json from")
mark_as_advanced( ACTS_NLOHMANN_JSON_SOURCE )

set( ACTS_THRUST_SOURCE
   "GIT_REPOSITORY;https://github.com/NVIDIA/thrust.git;GIT_TAG;v2.1.0" CACHE STRING "Source to take thrust from")
mark_as_advanced( ACTS_THRUST_SOURCE )

set( ACTS_BOOST_SOURCE
   "https://boostorg.jfrog.io/artifactory/main/release/${_acts_boost_recommended_version}/source/boost_${_acts_boost_recommended_version_}.tar.gz" CACHE STRING "Source to take boost from")
mark_as_advanced( ACTS_BOOST_SOURCE )

set( ACTS_EIGEN3_SOURCE
   "https://gitlab.com/libeigen/eigen/-/archive/${_acts_eigen3_version}/${_acts_eigen3_version}.tar.gz" CACHE STRING "Source to take eigen3 from")
mark_as_advanced( ACTS_EIGEN3_SOURCE )

set( ACTS_PYBIND11_SOURCE
   "GIT_REPOSITORY;https://github.com/pybind/pybind11.git;GIT_TAG;v2.10.1" CACHE STRING "Source to take pybind11 from")
mark_as_advanced( ACTS_PYBIND11_SOURCE )
