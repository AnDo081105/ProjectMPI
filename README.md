# Hybrid MPI+OpenMP DNA Sequence Alignment

## Overview
This implementation provides two DNA sequence alignment algorithms:

- Simple Sliding Window: Fast algorithm that counts exact base matches
- Smith-Waterman Local Alignment: More sophisticated algorithm for optimal local alignments with gaps and mismatches

The program uses a hybrid parallelization approach:

- MPI: Distributes query sequences across multiple processes
- OpenMP: Parallelizes alignment computations within each process

## Requirements
- MPI implementation (OpenMPT)
- Compiler support `-fopenmp`
- C++11 or later compiler
- FASTA format input files

## Usage
Implement in Makefile
- clean - Remove build artifacts"
- run - Run with default settings (simple algorithm)"
- run-sw - Run with Smith-Waterman algorithm"
- run-fasta - Run with custom FASTA files"
- test-performance - Quick performance test with different configurations"
- experiment-scaling - MPI process scaling experiment"
- experiment-threads - OpenMP thread scaling experiment"
- experiment-algorithms - Algorithm comparison experiment"
- help - Show this help message"