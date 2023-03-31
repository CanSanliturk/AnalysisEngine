#include <chrono>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include "Shape.h"
#include "Vector.h"
#include "UtilMethods.h"
#include "StructureSolver.h"

#define LOG(x) std::cout << x << "\n"
#define LOGW(x) std::wcout << x << "\n"
#define LOGIF(x, cond) if (cond) LOG(x)

int main()
{
    LOG(" __________________________________________________");
    LOG(" |                                                |");
    LOG(" |  3-Dimensional Finite Element Analysis Engine  |");
    LOG(" |          Created by Can Sanliturk              |");
    LOG(" |________________________________________________|");
    LOG("");

    // Start timer
    auto timenow =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    // Call test function (Later on, these guys will be moved to a unit test project)
    try
    {
        LOG(" Analysis completed without errors...");
    }
    catch (const std::runtime_error& e)
    {
        LOG(" Analysis cannot be completed...");
        LOG(" " << e.what());
    }

    // Log duration
    auto timenow2 =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    LOG(" Elapsed Time = " << timenow2 - timenow << " seconds\n");

    return 0;
}
