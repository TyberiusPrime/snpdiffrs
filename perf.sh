#!/bin/bash

perf record --call-graph=dwarf target/release/snpdiffrs
perf script | /machine/opt/infrastructure/repos/FlameGraph/stackcollapse-perf.pl | /machine/opt/infrastructure/repos/rust_unmangle/rust-unmangle | /machine/opt/infrastructure/repos/FlameGraph/flamegraph.pl > flame.svg
