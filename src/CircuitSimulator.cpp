#ifdef WITH_PYTHON
#include "matplotlibcpp.h"
#endif
#include <iostream>
#include "Maths/DynamicMatrix.hpp"
#include "CircuitElements/CircuitElements.hpp"
#include "CircuitSimulator/Simulator.hpp"
#include <thread>
#include <complex>
#include <chrono>

#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#include <shobjidl.h>
#include <ShObjIdl_core.h>
#endif


/// @brief main function to launch the circuit simulator
///
/// Contains platform specific code to handle directory changes
///
/// @param argc
/// @param argv[]
///
/// @return statuscode
int
main(int argc, char * argv[]) {
#ifdef _WIN32
    // credit for this way to get exe path:
    // https://gist.github.com/karolisjan/f9b8ac3ae2d41ec0ce70f2feac6bdfaf
    char buffer[MAX_PATH];
    GetModuleFileNameA(NULL, buffer, MAX_PATH);
    std::string::size_type pos = std::string(buffer).find_last_of("\\/");
    _chdir(std::string(buffer).substr(0, pos).c_str());

    // modification to get it as a wstring
    wchar_t bufferW[MAX_PATH];
    GetModuleFileNameW(NULL, bufferW, MAX_PATH);
    std::string::size_type posW = std::wstring(bufferW).find_last_of(L"\\/");
    std::wstring exePathW = std::wstring(bufferW).substr(0, posW);
#endif


    std::string filePath = "Netlists/Diode Test.netlist";
    if (argc > 1) {
        filePath = argv[1];
        std::cout << "Using netlist: " << filePath;
    } else {
#ifdef _WIN32
        // example code from:
        // https://docs.microsoft.com/en-us/windows/win32/learnwin32/example--the-open-dialog-box
        // modified to use ascii

        HRESULT hr = CoInitializeEx(NULL, COINIT_APARTMENTTHREADED |
                                              COINIT_DISABLE_OLE1DDE);
        if (SUCCEEDED(hr)) {
            IFileOpenDialog * pFileOpen;
            IShellItem * defaultFolder;

            // Create the FileOpenDialog object.
            hr = CoCreateInstance(CLSID_FileOpenDialog, NULL, CLSCTX_ALL,
                                  IID_IFileOpenDialog,
                                  reinterpret_cast<void **>(&pFileOpen));
            std::wstring defaultNetlistPath = exePathW + L"\\Netlists";
            SHCreateItemFromParsingName(defaultNetlistPath.c_str(), NULL,
                                        IID_IShellItem,
                                        reinterpret_cast<void **>(&defaultFolder));
            pFileOpen->SetDefaultFolder(defaultFolder);
            pFileOpen->SetFolder(defaultFolder);

            if (SUCCEEDED(hr)) {
                // Show the Open dialog box.
                hr = pFileOpen->Show(NULL);

                // Get the file name from the dialog box.
                if (SUCCEEDED(hr)) {
                    IShellItem * pItem;
                    hr = pFileOpen->GetResult(&pItem);
                    if (SUCCEEDED(hr)) {
                        PWSTR pszFilePath;
                        hr = pItem->GetDisplayName(SIGDN_FILESYSPATH, &pszFilePath);

                        // Display the file name to the user.
                        if (SUCCEEDED(hr)) {
                            // MessageBoxW(NULL, pszFilePath, L"File Path", MB_OK);
                            WideCharToMultiByte(CP_UTF8, 0, pszFilePath, -1, buffer,
                                                MAX_PATH, NULL, NULL);
                            filePath = std::string(buffer);
                            CoTaskMemFree(pszFilePath);
                        }
                        pItem->Release();
                    }
                }
                pFileOpen->Release();
            }
            CoUninitialize();
        }
        std::cout << "Using netlist: " << filePath;
#else
        std::cout << "No path given. Defaulting to using netlist: " << filePath;
#endif
    }
    std::cout << std::endl;

    SimulationEnvironment<double> env(filePath);

    env.simulate();

    return 0;
}
