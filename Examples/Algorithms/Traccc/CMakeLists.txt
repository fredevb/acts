add_library(
  ActsExamplesTraccc SHARED
  src/TracccChainAlgorithmHost.cpp
)

target_include_directories(
  ActsExamplesTraccc
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(
  ActsExamplesTraccc
  PUBLIC
    ActsPluginTraccc
    ActsExamplesFramework
    ActsExamplesDigitization
)

install(
  TARGETS ActsExamplesTraccc
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
