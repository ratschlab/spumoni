//
// Created by Sayan, Sayan on 10/22/25.
//

#include "spumoni_main.hpp"
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "parlay/primitives.h"
#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/io.h"

namespace py = pybind11;
using namespace std;

static inline char* to_c_str(std::string &s) {
    auto cs = new char[s.size()+1];
    strcpy(cs, s.c_str());
    return cs;
}

struct Alignment {
    string ctg;
    int r_st = 0, r_en = 0, strand = 1;
    float pres_frac = 0.0f;

    Alignment() = default;

    Alignment(const char *header, bool fwd, int start, float pres_frac, int qry_len):
        ctg(header), r_st(start), r_en(start + qry_len), strand(fwd?1:-1), pres_frac(pres_frac) {}
};

struct Request {
    int channel = 0;
    string id;
    string seq;
    Request() = default;
    Request(int channel, string &id, string &seq): channel(channel), id(id), seq(seq) {}
};

struct Response {
    int channel = 0;
    string id;
    Alignment alignment;
    Response() = default;
    Response(int channel, string &id, Alignment &alignment): channel(channel), id(id), alignment(alignment) {}
};

struct ResponseGenerator {
    parlay::sequence<Response> responses;
    explicit ResponseGenerator(parlay::sequence<Response> &responses): responses(responses) {}
    Response next() {
        if (responses.empty())
            throw py::stop_iteration();
        auto response = responses.back();
        responses.pop_back();
        return response;
    }
    ResponseGenerator& iter() {
        return *this;
    }
};

struct Index {
    void setup(vector<char*> argv) {
        printf("Arguments\n");
        for (auto& arg: argv) printf("%s ", arg);
        printf("\n");

        char *cmd = "run";
        argv.insert(argv.begin(), cmd);

        // Grab the run options, and validate they are not missing/don't make sense
        parse_run_options(argv.size(), argv.data(), &run_opts);
        run_opts.populate_types();

        // Add extension to reference file based on options
        if (run_opts.use_promotions)
            run_opts.ref_file += ".bin";
        else
            run_opts.ref_file += ".fa";

        tie(ms, max_value_thr) = setup_spumoni(&run_opts);
    }
    Index(vector<char*> argv) {
        setup(argv);
    }
    Index(const py::list& args) {
        vector<char*> argv;
        for (auto& arg: args) {
            auto s = arg.cast<string>();
            argv.push_back(to_c_str(s));
        }
        setup(argv);
    }

    void validate() {
        /* Checks the options for the run command, and makes sure it has everything it needs */
      if (run_opts.ref_file.empty()){FATAL_WARNING("A reference file must be provided.");}
      if (run_opts.result_type == NOT_CHOSEN) {FATAL_WARNING("An output type of MS or PML must be specified, only one can be used at a time.");}

      // Make sure provided files are valid
      if (!is_file(run_opts.ref_file)) {FATAL_ERROR("The following path is not valid: %s (remember to only specify output prefix)", (run_opts.ref_file).data());}

      // Make sure reference file is a valid type
      if (!run_opts.is_general_text && run_opts.ref_type == NOT_SET) {FATAL_ERROR("Reference file is an unrecognized type. It needs to be a\n"
                                                    "       FASTA file or binary file produced by spumoni build.");}

      // Verify doc array is available, if needed
      if (run_opts.use_doc && !is_file(run_opts.ref_file+".doc"))
        FATAL_WARNING("document array file (%s) is not present, so it cannot be used.", (run_opts.ref_file+".doc").data());

      switch (run_opts.result_type) {
        case MS:
            if (!is_file(run_opts.ref_file+".thrbv.ms"))
            {FATAL_WARNING("The index required for this computation is not available, please use spumoni build.");} break;
        case PML:
            if (!is_file(run_opts.ref_file+".thrbv.spumoni"))
            {FATAL_WARNING("The index required for this computation is not available, please use spumoni build.");} break;
        default:
            FATAL_WARNING("An output type with -M or -P must be specified, only one can be used at a time."); break;
      }

      // Check the values for k and w
      if (run_opts.k > 4) {FATAL_WARNING("small window size (k) cannot be larger than 4 characters.");}
      if (run_opts.w < run_opts.k) {FATAL_WARNING("large window size (w) should be larger than the small window size (k)");}

      // Check if we choose the minimizer type correctly
      if (run_opts.min_digest) {
        if (run_opts.use_promotions && run_opts.use_dna_letters) {FATAL_ERROR("Only one type of minimizer can be specified from either -m or -a.");}
        if (!run_opts.use_promotions && !run_opts.use_dna_letters) {FATAL_ERROR("A minimizer type must be specified using -m or -a.");}
      } else {
        if (run_opts.use_promotions || run_opts.use_dna_letters) {FATAL_ERROR("A minimizer type should not be specified if intending not to use minimizer digestion.");}
      }

      // Check the size of KS-test region
      if (run_opts.bin_size < 50 || run_opts.bin_size > 400) {FATAL_WARNING("the bin size used is not optimal. Re-run using a value between 50 and 400.");}
        printf("Spumoni Index arguments are valid.\n");
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
        parlay::for_each(parlay::iota(nr), [&](size_t i){
            auto sequence = sequences[i].cast<string>();
            results[i] = query(sequence);
        });
        return results;
    }

    ResponseGenerator query_stream(const py::iterator& reads) {
        parlay::sequence<Request> requests;
        for (auto &read: reads) {
            auto request = read.cast<Request>();
            requests.push_back(request);
        }
        auto responses = parlay::tabulate(requests.size(), [&](size_t i) {
            auto alignment = query(requests[i].seq);
            return Response(requests[i].channel, requests[i].id, alignment);
        });
        return ResponseGenerator(responses);
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
            .def(py::init<const py::list&>())
            .def("validate", &Index::validate)
            .def("query", &Index::query)
            .def("query_batch", &Index::query_batch)
            .def("query_stream", &Index::query_stream);

    py::class_<Request>(m, "Request")
            .def(py::init<int, string&, string&>(), py::arg("channel"), py::arg("id"), py::arg("seq"))
            .def_readwrite("channel", &Request::channel)
            .def_readwrite("id", &Request::id)
            .def_readwrite("seq", &Request::seq);

    py::class_<Response>(m, "Response")
            .def(py::init<int, string&, Alignment&>(), py::arg("channel"), py::arg("id"), py::arg("alignment"))
            .def_readwrite("channel", &Response::channel)
            .def_readwrite("id", &Response::id)
            .def_readwrite("alignment", &Response::alignment);

    py::class_<ResponseGenerator>(m, "ResponseGenerator")
            .def("__iter__", &ResponseGenerator::iter)
            .def("__next__", &ResponseGenerator::next);
}