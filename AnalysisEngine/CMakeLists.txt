set(PROJECT_NAME AnalysisEngine)

################################################################################
# Source groups
################################################################################
set(Header_Files
    "include/Shape.h"
    "include/XYData.h"
)
source_group("Header Files" FILES ${Header_Files})

set(src
    "AnalysisEngine.cpp"
    "src/Material.cpp"
    "src/Shape.cpp"
    "src/XYData.cpp"
)
source_group("src" FILES ${src})

set(src__BoundaryConditions
    "include/Hinge.h"
    "include/Restraint.h"
)
source_group("src\\BoundaryConditions" FILES ${src__BoundaryConditions})

set(src__Elements
    "include/Element.h"
    "include/FrameMember.h"
    "include/SerendipityShell.h"
    "include/Shell.h"
    "include/TrussMember.h"
    "src/FrameMember.cpp"
    "src/SerendipityShell.cpp"
    "src/Shell.cpp"
    "src/TrussMember.cpp"
)
source_group("src\\Elements" FILES ${src__Elements})

set(src__Loadings
    "include/DistributedLoad.h"
    "include/NodalLoad.h"
)
source_group("src\\Loadings" FILES ${src__Loadings})

set(src__Material
    "include/Material.h"
)
source_group("src\\Material" FILES ${src__Material})

set(src__Nodes
    "include/Node.h"
    "src/Node.cpp"
)
source_group("src\\Nodes" FILES ${src__Nodes})

set(src__Section
    "include/Piece.h"
    "include/Section.h"
    "src/Piece.cpp"
    "src/Section.cpp"
)
source_group("src\\Section" FILES ${src__Section})

set(src__Solvers
    "include/StructureSolver.h"
    "src/StructureSolver.cpp"
)
source_group("src\\Solvers" FILES ${src__Solvers})

set(src__Structure
    "include/Structure.h"
    "src/Structure.cpp"
)
source_group("src\\Structure" FILES ${src__Structure})

set(src__Utils
    "include/GeometryHelper.h"
    "include/Matrix.h"
    "include/UtilMethods.h"
    "include/Vector.h"
    "include/XYPoint.h"
    "include/XYZPoint.h"
    "src/GeometryHelper.cpp"
    "src/UtilMethods.cpp"
    "src/Vector.cpp"
)
source_group("src\\Utils" FILES ${src__Utils})

set(ALL_FILES
    ${Header_Files}
    ${src}
    ${src__BoundaryConditions}
    ${src__Elements}
    ${src__Loadings}
    ${src__Material}
    ${src__Nodes}
    ${src__Section}
    ${src__Solvers}
    ${src__Structure}
    ${src__Utils}
)

################################################################################
# Target
################################################################################
add_executable(${PROJECT_NAME} ${ALL_FILES})

use_props(${PROJECT_NAME} "${CMAKE_CONFIGURATION_TYPES}" "${DEFAULT_CXX_PROPS}")
set(ROOT_NAMESPACE AnalysisEngine)

set_target_properties(${PROJECT_NAME} PROPERTIES
    VS_GLOBAL_KEYWORD "Win32Proj"
)
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    set_target_properties(${PROJECT_NAME} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
    )
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
    set_target_properties(${PROJECT_NAME} PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION_RELEASE "TRUE"
    )
endif()
################################################################################
# Compile definitions
################################################################################
if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Debug>:"
            "_DEBUG"
        ">"
        "$<$<CONFIG:Release>:"
            "NDEBUG"
        ">"
        "_CONSOLE;"
        "UNICODE;"
        "_UNICODE"
    )
elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
    target_compile_definitions(${PROJECT_NAME} PRIVATE
        "$<$<CONFIG:Debug>:"
            "_DEBUG"
        ">"
        "$<$<CONFIG:Release>:"
            "NDEBUG"
        ">"
        "WIN32;"
        "_CONSOLE;"
        "UNICODE;"
        "_UNICODE"
    )
endif()

################################################################################
# Compile and link options
################################################################################
if(MSVC)
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /Oi;
                /Gy
            >
            /permissive-;
            /sdl;
            /W3;
            ${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
            /bigobj;
            ${DEFAULT_CXX_EXCEPTION_HANDLING}
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
        target_compile_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Release>:
                /Oi;
                /Gy
            >
            /permissive-;
            /sdl;
            /W3;
            ${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
            ${DEFAULT_CXX_EXCEPTION_HANDLING}
        )
    endif()
    if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
        target_link_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Debug>:
                /INCREMENTAL
            >
            $<$<CONFIG:Release>:
                /OPT:REF;
                /OPT:ICF;
                /INCREMENTAL:NO
            >
            /DEBUG;
            /SUBSYSTEM:CONSOLE
        )
    elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
        target_link_options(${PROJECT_NAME} PRIVATE
            $<$<CONFIG:Debug>:
                /INCREMENTAL
            >
            $<$<CONFIG:Release>:
                /OPT:REF;
                /OPT:ICF;
                /INCREMENTAL:NO
            >
            /DEBUG;
            /SUBSYSTEM:CONSOLE
        )
    endif()
endif()

################################################################################
# Dependencies
################################################################################
# if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
#     set(ADDITIONAL_LIBRARY_DEPENDENCIES
#         "$(SolutionDir)ArmadilloSolver\OpenBLAS\lib_win64\libopenblas"
#     )
# elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
#     set(ADDITIONAL_LIBRARY_DEPENDENCIES
#         "$(SolutionDir)ArmadilloSolver\OpenBLAS\lib_win64\libopenblas"
#     )
# endif()
# target_link_libraries(${PROJECT_NAME} PRIVATE "${ADDITIONAL_LIBRARY_DEPENDENCIES}")

# if("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x64")
#     target_link_directories(${PROJECT_NAME} PRIVATE
#         "${CMAKE_SOURCE_DIR}/ArmadilloSolver/OpenBLAS/lib_win64"
#     )
# elseif("${CMAKE_VS_PLATFORM_NAME}" STREQUAL "x86")
#     target_link_directories(${PROJECT_NAME} PRIVATE
#         "${CMAKE_SOURCE_DIR}/ArmadilloSolver/OpenBLAS/lib_win64"
#     )
# endif()

