add_executable(api_check api_check.c)
target_link_libraries(api_check PRIVATE aacdec)

add_executable(compare_check compare_check.c)
target_link_libraries(compare_check PRIVATE aacdec)

add_executable(he_check he_check.c)
target_link_libraries(he_check PRIVATE aacdec)

add_executable(usac_check usac_check.c)
target_link_libraries(usac_check PRIVATE aacdec)
if(UNIX)
    target_link_libraries(usac_check PRIVATE m)
endif()

add_executable(ps_check ps_check.c)
target_link_libraries(ps_check PRIVATE aacdec)

add_executable(ld_check ld_check.c)
target_link_libraries(ld_check PRIVATE aacdec)
if(UNIX)
    target_link_libraries(ld_check PRIVATE m)
endif()

add_executable(eld_check eld_check.c)
target_link_libraries(eld_check PRIVATE aacdec)
if(UNIX)
    target_link_libraries(eld_check PRIVATE m)
endif()

add_executable(drc_check drc_check.c)
target_link_libraries(drc_check PRIVATE aacdec)
target_include_directories(drc_check PRIVATE ${CMAKE_SOURCE_DIR}/src)
if(UNIX)
    target_link_libraries(drc_check PRIVATE m)
endif()

add_executable(conceal_check conceal_check.c)
target_link_libraries(conceal_check PRIVATE aacdec)
if(UNIX)
    target_link_libraries(conceal_check PRIVATE m)
endif()

add_executable(er_check er_check.c)
target_link_libraries(er_check PRIVATE aacdec)
if(UNIX)
    target_link_libraries(er_check PRIVATE m)
endif()

add_executable(crc_check crc_check.c)
target_link_libraries(crc_check PRIVATE aacdec)
if(UNIX)
    target_link_libraries(crc_check PRIVATE m)
endif()

add_executable(ps_unit ps_unit.c)
target_include_directories(ps_unit PRIVATE ${CMAKE_SOURCE_DIR}/src)
if(UNIX)
    target_link_libraries(ps_unit PRIVATE m)
endif()

add_test(NAME api_check_stereo
         COMMAND api_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.aac)
add_test(NAME api_check_surround
         COMMAND api_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_51.aac)
add_test(NAME compare_stereo
         COMMAND compare_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.aac
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.s16)
add_test(NAME compare_surround
         COMMAND compare_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_51.aac
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_51.s16)
add_test(NAME compare_latm
         COMMAND compare_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.latm
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.s16 32)
add_test(NAME compare_adif
         COMMAND compare_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.adif
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.s16 32)
add_test(NAME compare_raw
         COMMAND compare_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.raw
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.s16 32 raw
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.asc)
add_test(NAME he_check_stereo
         COMMAND he_check ${CMAKE_CURRENT_SOURCE_DIR}/data/he_stereo.aac)
add_test(NAME usac_check_stereo
         COMMAND usac_check ${CMAKE_CURRENT_SOURCE_DIR}/data/usac_stereo.loas)
add_test(NAME ps_check_stereo
         COMMAND ps_check ${CMAKE_CURRENT_SOURCE_DIR}/data/he_v2.aac)
add_test(NAME ps_unit
         COMMAND ps_unit)
add_test(NAME ld_check
         COMMAND ld_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_ld.latm
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_ld_ref.pcm)
add_test(NAME eld_check
         COMMAND eld_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_eld.latm
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_eld_ref.pcm)
add_test(NAME drc_check
         COMMAND drc_check)
add_test(NAME conceal_check
         COMMAND conceal_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.aac
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.s16
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_stereo.asc)
add_test(NAME er_check
         COMMAND er_check)
add_test(NAME crc_check
         COMMAND crc_check ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_adts_crc4.aac
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_adts_crc_ref.s16)
add_test(NAME compare_latm_amv1
         COMMAND compare_check
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_latm_amv1.loas
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_latm_ref.s16 32)
add_test(NAME compare_latm_2sub
         COMMAND compare_check
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_latm_2sub.loas
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_latm_ref.s16 32)
add_test(NAME compare_adts_crc1
         COMMAND compare_check
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_adts_crc1.aac
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_adts_crc_ref.s16 32)
add_test(NAME compare_adts_crc4
         COMMAND compare_check
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_adts_crc4.aac
                 ${CMAKE_CURRENT_SOURCE_DIR}/data/tiny_adts_crc_ref.s16 32)
