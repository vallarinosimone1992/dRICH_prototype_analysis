cmake_minimum_required(VERSION 3.5)

project(dRICH_reco)

find_package(ROOT)

add_library(dRICH_lib STATIC
  lib/fillMAPS.C
  lib/getChannel.C
  lib/MAPMTposition.C
  lib/handleFile.C
  lib/integrate.C
  lib/dRICH.C
)

target_include_directories(dRICH_lib
  PUBLIC 
    ${PROJECT_SOURCE_DIR}/lib
)


add_executable(reco
  plugin/main_test.C
)

target_link_libraries(reco
  PRIVATE
      dRICH_lib
)
