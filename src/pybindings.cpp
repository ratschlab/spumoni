//
// Created by Sayan, Sayan on 10/22/25.
//

#include "spumoni_main.hpp"
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <numeric>
#include <parallel/algorithm>

#define GLIBCXX_PARALLEL

namespace py = pybind11;
using namespace std;

static inline char* to_c_str(std::string &s) {
    auto cs = new char[s.size()+1];
    strcpy(cs, s.c_str());
    return cs;
}

static vector<char*> kwargs_to_argv(const py::args& args, const py::kwargs& kwargs) {
    vector<string> arg_strings;

    // First argument is conventionally the program name
    arg_strings.emplace_back("program_name");

    // Handle positional arguments (flags)
    for (auto& arg : args) {
        string flag = "--" + arg.cast<string>();  // Convert to "--flag" format
        arg_strings.emplace_back(flag);
    }

    // Handle keyword arguments ("--key=value")
    for (auto& item : kwargs) {
        auto key = item.first.cast<string>();
        string value = py::str(item.second);  // Convert value to string
        if (key.size()==1) arg_strings.emplace_back("-" + key + "=" + value);
        else arg_strings.emplace_back("--" + key + "=" + value);
    }

    vector<char*> argv;  // Store pointers to C-style strings
    for (auto &s: arg_strings)
        argv.push_back(to_c_str(s));

    // Now `argc` and `argv.data()` can be used in a function expecting C-style args
    return argv;
}

struct Alignment {
    string ctg;
    int r_st = 0, r_en = 0, strand = 1;
    float pres_frac = 0.0f;

    Alignment() = default;

    Alignment(const char *header, bool fwd, int start, float pres_frac, int qry_len):
        ctg(header), r_st(start), r_en(start + qry_len), strand(fwd?1:-1), pres_frac(pres_frac) {}
};

struct Index {
    Index(const py::args& args, const py::kwargs& kwargs) {
        auto argv = kwargs_to_argv(args, kwargs);

        // Grab the run options, and validate they are not missing/don't make sense
        parse_run_options(argv.size(), argv.data(), &run_opts);
        run_opts.populate_types();
        run_opts.validate();

        // Add extension to reference file based on options
        if (run_opts.use_promotions)
            run_opts.ref_file += ".bin";
        else
            run_opts.ref_file += ".fa";

        tie(ms, max_value_thr) = setup_spumoni(&run_opts);
    }

    Alignment query(string &sequence) {
        const auto& [found, pres] = classify_read(ms, sequence, run_opts.k, run_opts.w, run_opts.use_promotions,
            run_opts.use_dna_letters, run_opts.bin_size, run_opts.result_type == PML, max_value_thr);
        if (found) return {"Found", true, 0, static_cast<float>(pres), static_cast<int>(sequence.length())};
        else return {"*", true, 0, 0, static_cast<int>(sequence.length())};
    }

    vector<Alignment> query_batch(const py::list& sequences) {
        auto nr = sequences.size();
        vector<Alignment> results(nr);

        vector<size_t> indices(nr);
        iota(indices.begin(), indices.end(), 0);

        __gnu_parallel::for_each(indices.begin(), indices.end(),
            [&](size_t i) {
                auto sequence = sequences[i].cast<string>();
                results[i] = query(sequence);
            });

        return results;
    }

private:
    SpumoniRunOptions run_opts;
    void *ms;
    size_t max_value_thr;

};

PYBIND11_MODULE(_core, m) {
    py::class_<Alignment>(m, "Alignment")
            .def(py::init<>())  // Default constructor
            .def(py::init<const char*, bool, int, float, int>(),  // Parameterized constructor
                 py::arg("header"), py::arg("fwd"), py::arg("start"),
                 py::arg("pres_frac"), py::arg("qry_len"))
            .def_readonly("ctg", &Alignment::ctg)
            .def_readonly("r_st", &Alignment::r_st)
            .def_readonly("r_en", &Alignment::r_en)
            .def_readonly("strand", &Alignment::strand)
            .def_readonly("pres_frac", &Alignment::pres_frac);

    py::class_<Index>(m, "Index")
            .def(py::init<const py::args&, const py::kwargs&>())
            .def("query", &Index::query)
            .def("query_batch", &Index::query_batch);
}