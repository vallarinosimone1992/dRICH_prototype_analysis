cmake_minimum_required(VERSION 3.5)

project(dRICH_reco)

find_package(ROOT 6.16 CONFIG REQUIRED)

include("${ROOT_USE_FILE}")

separate_arguments(ROOT_EXE_LINKER_FLAGS)

add_library(dRICH_lib
      STATIC
      lib/getChannel.cxx
      lib/MAPMTposition.cxx
      lib/handleFile.cxx
      lib/integrate.cxx
      lib/dRICH.cxx
      lib/fillMAPS.cxx
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
      PUBLIC
      ${ROOT_LIBRARIES}
      ${ROOT_EXE_LINKER_FLAGS}
      ~                                                                           
      )
