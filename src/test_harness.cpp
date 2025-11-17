//
// Created by Goswami, Sayan on 11/14/25.
//

#include <vector>
#include <pybind11/embed.h> // Key for embedding Python
#include <iostream>

namespace py = pybind11;

using namespace std;

int main(int argc, char *argv[]) {
    // Start the Python interpreter.
    // This is the C++ equivalent of running `python`
    // The `scoped_interpreter` handles Py_Initialize() and Py_Finalize()
    py::scoped_interpreter guard{};
    try {
        // --- Setup ---
        // Add the current directory to Python's path
        // This is so Python can find your .so and .py files
        py::module_::import("sys").attr("path").attr("append")(".");

        // We will import your python helper script to build the list
        py::module_ psp = py::module_::import("_core");

        // const char * argv[] = {"run", "--ref", "/scratch/NASExperiments/tmp/zymo", "--threads", "16", "--PML", "--minimizer-alphabet", "--classify"};
        py::list py_list;
        for (int i = 1; i < argc; i++) py_list.append(argv[i]);

        // Get the class 'Index' from the module 'example'
        py::object IndexClass = psp.attr("Index");
        py::object index_instance = IndexClass(py_list);
        index_instance.attr("validate")();


    } catch (py::error_already_set &e) {
        // Handle any Python exceptions
        std::cerr << "--- Python Error ---" << std::endl;
        std::cerr << e.what() << std::endl;
        std::cerr << "----------------------" << std::endl;
        return 1;
    }

    return 0;
}