set(CEA_EXCEL_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

add_library(cea_excel SHARED
    src/cea_excel.c
)
add_library(cea::excel ALIAS cea_excel)

target_compile_definitions(cea_excel PRIVATE
    CEA_EXCEL_BUILD
    CEA_EXCEL_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}
    CEA_EXCEL_VERSION_MINOR=${PROJECT_VERSION_MINOR}
    CEA_EXCEL_VERSION_PATCH=${PROJECT_VERSION_PATCH}
)
target_include_directories(cea_excel PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/cea/excel>
)
set_target_properties(cea_excel PROPERTIES
    OUTPUT_NAME cea_excel
    LIBRARY_OUTPUT_DIRECTORY ${CEA_EXCEL_OUTPUT_DIR}
    RUNTIME_OUTPUT_DIRECTORY ${CEA_EXCEL_OUTPUT_DIR}
    ARCHIVE_OUTPUT_DIRECTORY ${CEA_EXCEL_OUTPUT_DIR}
)

if(APPLE)
    set_target_properties(cea_excel PROPERTIES
        INSTALL_NAME_DIR "@rpath"
    )
endif()

install(TARGETS cea_excel DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES include/cea_excel.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cea/excel)

if (CEA_BUILD_TESTING)
    add_executable(cea_excel_version_test tests/cea_excel_version_test.c)
    target_include_directories(cea_excel_version_test PRIVATE include)
    target_compile_definitions(cea_excel_version_test PRIVATE
        CEA_EXCEL_TEST_LIBRARY="$<TARGET_FILE:cea_excel>"
    )

    if(UNIX AND NOT APPLE)
        target_link_libraries(cea_excel_version_test PRIVATE dl)
    endif()

    add_test(
        NAME cea_excel_version_test
        COMMAND cea_excel_version_test
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()
