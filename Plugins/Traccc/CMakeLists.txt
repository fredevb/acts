add_library(
  ActsPluginTraccc SHARED
  src/TracccConversion.cpp)

target_include_directories(
  ActsPluginTraccc
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(
  ActsPluginTraccc
  PUBLIC
    ActsCore
    traccc::core
    traccc::io
    vecmem::core
    detray::utils
    detray::io)

install(
  TARGETS ActsPluginTraccc
  EXPORT ActsPluginTracccTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(
  DIRECTORY include/Acts
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
