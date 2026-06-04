# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**SPUMONI** (Streaming PseUdo MONI) is a bioinformatics tool for rapid read classification on DNA sequencing reads. It computes matching statistics (MS) and pseudo-matching lengths (PML) as discriminative features against a reference index. PMLs are ~3x faster and generally more accurate than MSs for binary classification.

The tool follows a two-phase workflow:
1. **Build phase** (`spumoni build`) — constructs a reference index using prefix-free parsing to build a run-length encoded BWT
2. **Query/classification phase** (`spumoni run`) — computes MS/PML values and classifies reads using a KS-test against an empirical null distribution

## Build Commands

```bash
# Configure and build
mkdir build && cd build
cmake ..
make install
# Executable at: build/spumoni (or build/install/)
```

Dependencies: `zlib`, `cmake`, `gcc`. All other third-party deps (sdsl-lite, r-index, Big-BWT, pfp-thresholds, bigrepair, shaped_slp, bonsai, klib) are fetched via CMake FetchContent at configure time.

## Running

```bash
# Build an index (PML + MS, with minimizer digestion)
./spumoni build -r genome.fa -M -P -m -o /path/to/index

# Run classification
./spumoni run -r /path/to/index -p reads.fa -P -c
```

## Testing

There is no automated test suite. `src/test_harness.cpp` embeds a Python interpreter to test the Python bindings (`pyspumoni._core`). Testing is done manually against benchmark datasets.

## Architecture

### Key Components

**`include/spumoni_main.hpp`** — Central hub: all shared enums (`SparseOptCodes`, `DigestOptCodes`), the `SpumoniHelperPrograms` struct listing helper binary paths, and `SpumoniRunOptions`/`SpumoniBuildOptions` option structs.

**`src/spumoni.cpp`** — CLI entry point. Parses subcommands (`build`/`run`), populates options structs, and dispatches to `build_spumoni_ms_main()` / `run_spumoni_main()` etc.

**`include/compute_ms_pml.hpp` / `src/compute_ms_pml.cpp`** — Core algorithms for MS and PML computation against the run-length encoded BWT index. This is the largest and most performance-critical file (~67KB).

**`include/refbuilder.hpp` / `src/refbuilder.cpp`** — Reference preprocessing: FASTA parsing, minimizer digestion (alphabet-promoted or DNA-letter), reverse complement generation, and null read sampling.

**`include/emp_null_database.hpp` / `src/emp_null_database.cpp`** — Builds and stores the empirical null distribution (MS/PML values from random null reads sampled from the reference). Used by the KS-test for classification.

**`include/ks_test.hpp` / `src/ks_test.cpp`** — Kolmogorov-Smirnov test for classification. Thresholds: `KS_STAT_MS_THR = 0.25`, `KS_STAT_PML_THR = 0.10`.

**`include/batch_loader.hpp` / `src/batch_loader.cpp`** — Streams FASTA/FASTQ input in batches (Kraken2-inspired); supports multi-threaded load balancing.

**`include/doc_array.hpp` / `src/doc_array.cpp`** — Optional document array mapping BWT positions back to source documents, enabling per-genome attribution in multi-genome indexes.

**`include/ms_rle_string.hpp`** — Extends `ri::rle_string` (from the r-index library) with matching-statistics–specific operations. The central data structure for index queries.

**`include/thresholds_ds.hpp`** — Stores threshold positions from prefix-free parsing, used to bound suffix array samples during MS/PML computation.

### Python Bindings

**`src/pybindings.cpp`** — pybind11 bindings exposing `Index`, `Alignment`, `Request`, `Response`, and `ResponseGenerator` to Python.

**`src/pyspumoni/aligner.py`** — High-level Python `Aligner` class wrapping the C++ `Index`. Converts Python kwargs to CLI args and provides `map_reads()` for streaming processing.

```python
from pyspumoni import Aligner
aligner = Aligner(ref='index_prefix', PML=True, threads=8)
```

Thread count is controlled via the `PARLAY_NUM_THREADS` environment variable (parlaylib).

### Important Constants (`include/spumoni_main.hpp`)

```cpp
#define THRBYTES 5          // Bytes for threshold storage
#define SSABYTES 5          // Bytes for SA sample storage
#define NULL_READ_CHUNK 150 // bp per null read sample
#define NUM_NULL_READS 800  // Number of null reads to sample
```

### Index File Extensions

| Extension | Contents |
|---|---|
| `.thrbv.ms` | Matching statistics index |
| `.thrbv.spumoni` | Pseudo-matching lengths index |
| `.doc` | Document array (multi-genome) |
| `.fdi` | FASTA document index |
| `.fa` / `.bin` | Reference (FASTA or binary minimizers) |
