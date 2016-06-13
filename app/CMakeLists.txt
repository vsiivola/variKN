FILE(GLOB srcs "*.cc")
set(libs "varikn")

foreach (file ${srcs})
        GET_FILENAME_COMPONENT(EXENAME ${file} NAME_WE)
        add_executable(${EXENAME} ${file})
        target_link_libraries(${EXENAME} ${libs})
        install(TARGETS ${EXENAME} DESTINATION bin)
endforeach(file)
