
add_subdirectory("qlaw")

include_directories(include ${CMAKE_CURRENT_SOURCE_DIR}/qlaw)

add_executable(qLawSolver main.cpp)
target_link_libraries(qLawSolver qLaw)