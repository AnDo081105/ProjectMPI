# Makefile for Hybrid MPI+OpenMP DNA Sequence Alignment

CXX = mpicxx
CXXFLAGS = -std=c++14 -O2 -Wall -fopenmp
LDFLAGS = -fopenmp
TARGET = mpi_dna_alignment
SOURCE = mpi_js.cpp

# Default target
all: $(TARGET)

# Build the DNA alignment executable with OpenMP support
$(TARGET): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE) $(LDFLAGS)

# Clean build artifacts
clean:
	rm -f $(TARGET) $(TARGET).exe

# Run with default settings using coronavirus FASTA file (4 processes, 2 threads per process)
run: $(TARGET)
	@echo "Running with coronavirus.fasta as both reference and query source"
	@echo "Using coronavirus.fasta for reference and queries"
	OMP_NUM_THREADS=2 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta simple

# Run with Smith-Waterman algorithm using coronavirus FASTA file
run-sw: $(TARGET)
	@echo "Running Smith-Waterman algorithm with coronavirus.fasta"
	@echo "Using coronavirus.fasta for reference and queries"
	OMP_NUM_THREADS=2 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw

# Run with FASTA files (configurable)
run-fasta: $(TARGET)
	@echo "Usage: make run-fasta REF=<reference.fasta> QUERIES=<queries.fasta> [ALGO=<algorithm>] [NP=<processes>] [THREADS=<threads>]"
	@echo "Example: make run-fasta REF=coronavirus.fasta QUERIES=coronavirus.fasta ALGO=sw NP=4 THREADS=2"
	@if [ -z "$(REF)" ] || [ -z "$(QUERIES)" ]; then \
		echo "Error: REF and QUERIES parameters are required"; \
		echo "Example: make run-fasta REF=coronavirus.fasta QUERIES=coronavirus.fasta"; \
		exit 1; \
	fi
	OMP_NUM_THREADS=$(or $(THREADS),2) mpirun -np $(or $(NP),4) -hostfile ./cluster ./$(TARGET) $(REF) $(QUERIES) $(or $(ALGO),simple)



# Run performance test with different hybrid configurations
test-performance: $(TARGET)
	@echo "Testing hybrid MPI+OpenMP DNA alignment performance with coronavirus.fasta..."
	@echo "Using coronavirus.fasta as both reference and query source"
	@echo ""
	@echo "2 MPI processes, 2 threads each:"
	OMP_NUM_THREADS=2 mpirun -np 2 ./$(TARGET) coronavirus.fasta coronavirus.fasta simple
	@echo ""
	@echo "4 MPI processes, 2 threads each (simple algorithm):"
	OMP_NUM_THREADS=2 mpirun -np 4 ./$(TARGET) coronavirus.fasta coronavirus.fasta simple
	@echo ""
	@echo "4 MPI processes, 2 threads each (Smith-Waterman algorithm):"
	OMP_NUM_THREADS=2 mpirun -np 4 ./$(TARGET) coronavirus.fasta coronavirus.fasta sw
	@echo ""
	@echo "8 MPI processes, 2 threads each:"
	OMP_NUM_THREADS=2 mpirun -np 8 ./$(TARGET) coronavirus.fasta coronavirus.fasta sw

# Comprehensive scaling experiment
experiment-scaling: $(TARGET)
	@echo "=== Scaling Experiment: MPI Process Scaling ===" > scaling_results.txt
	@echo "Date: $$(date)" >> scaling_results.txt
	@echo "Algorithm: Both Simple and Smith-Waterman" >> scaling_results.txt
	@echo "Dataset: coronavirus.fasta" >> scaling_results.txt
	@echo "OpenMP Threads per process: 2" >> scaling_results.txt
	@echo "" >> scaling_results.txt
	@echo "Running scaling experiment with different MPI process counts..."
	@echo "=== 1 Process (Serial baseline) ===" >> scaling_results.txt
	@echo "Simple Algorithm:" >> scaling_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 1 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta simple" 2>&1 | tee -a scaling_results.txt
	@echo "" >> scaling_results.txt
	@echo "Smith-Waterman Algorithm:" >> scaling_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 1 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a scaling_results.txt
	@echo "" >> scaling_results.txt
	@echo "=== 2 Processes ===" >> scaling_results.txt
	@echo "Simple Algorithm:" >> scaling_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 2 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta simple" 2>&1 | tee -a scaling_results.txt
	@echo "" >> scaling_results.txt
	@echo "Smith-Waterman Algorithm:" >> scaling_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 2 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a scaling_results.txt
	@echo "" >> scaling_results.txt
	@echo "=== 4 Processes ===" >> scaling_results.txt
	@echo "Simple Algorithm:" >> scaling_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta simple" 2>&1 | tee -a scaling_results.txt
	@echo "" >> scaling_results.txt
	@echo "Smith-Waterman Algorithm:" >> scaling_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a scaling_results.txt
	@echo "" >> scaling_results.txt
	@echo "=== 8 Processes ===" >> scaling_results.txt
	@echo "Simple Algorithm:" >> scaling_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 8 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta simple" 2>&1 | tee -a scaling_results.txt
	@echo "" >> scaling_results.txt
	@echo "Smith-Waterman Algorithm:" >> scaling_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 8 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a scaling_results.txt
	@echo "Scaling experiment completed! Results saved to scaling_results.txt"

# OpenMP thread scaling experiment  
experiment-threads: $(TARGET)
	@echo "=== Thread Scaling Experiment: OpenMP Thread Scaling ===" > thread_results.txt
	@echo "Date: $$(date)" >> thread_results.txt
	@echo "MPI Processes: 4" >> thread_results.txt
	@echo "Algorithm: Smith-Waterman" >> thread_results.txt
	@echo "Dataset: coronavirus.fasta" >> thread_results.txt
	@echo "" >> thread_results.txt
	@echo "Running thread scaling experiment with different OpenMP thread counts..."
	@echo "=== 1 Thread per process ===" >> thread_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=1 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a thread_results.txt
	@echo "" >> thread_results.txt
	@echo "=== 2 Threads per process ===" >> thread_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a thread_results.txt
	@echo "" >> thread_results.txt
	@echo "=== 4 Threads per process ===" >> thread_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=4 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a thread_results.txt
	@echo "" >> thread_results.txt
	@echo "=== 8 Threads per process ===" >> thread_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=8 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a thread_results.txt
	@echo "Thread scaling experiment completed! Results saved to thread_results.txt"

# Algorithm comparison experiment
experiment-algorithms: $(TARGET)
	@echo "=== Algorithm Comparison Experiment ===" > algorithm_results.txt
	@echo "Date: $$(date)" >> algorithm_results.txt
	@echo "Configuration: 4 MPI processes, 2 OpenMP threads each" >> algorithm_results.txt
	@echo "Dataset: coronavirus.fasta" >> algorithm_results.txt
	@echo "" >> algorithm_results.txt
	@echo "Running algorithm comparison experiment..."
	@echo "=== Simple Sliding Window Algorithm ===" >> algorithm_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta simple" 2>&1 | tee -a algorithm_results.txt
	@echo "" >> algorithm_results.txt
	@echo "=== Smith-Waterman Algorithm ===" >> algorithm_results.txt
	@/usr/bin/time -f "real %e\nuser %U\nsys %S" bash -c "OMP_NUM_THREADS=2 mpirun -np 4 -hostfile ./cluster ./$(TARGET) coronavirus.fasta coronavirus.fasta sw" 2>&1 | tee -a algorithm_results.txt
	@echo "Algorithm comparison completed! Results saved to algorithm_results.txt"


# Help target
help:
	@echo "=== Hybrid MPI+OpenMP DNA Sequence Alignment Makefile ==="
	@echo ""
	@echo "Available targets:"
	@echo "  all - Build DNA alignment executable"
	@echo "  $(TARGET) - Build DNA alignment executable"
	@echo "  clean - Remove build artifacts"
	@echo "  run - Run with default settings (simple algorithm)"
	@echo "  run-sw - Run with Smith-Waterman algorithm"
	@echo "  run-fasta - Run with custom FASTA files"
	@echo "  test-performance - Quick performance test with different configurations"
	@echo "  experiment-scaling - MPI process scaling experiment"
	@echo "  experiment-threads - OpenMP thread scaling experiment"
	@echo "  experiment-algorithms - Algorithm comparison experiment"
	@echo "  help - Show this help message"
	@echo ""
	@echo "Experimental Targets:"
	@echo "  experiment-scaling - Test MPI process scaling (1,2,4,8 processes)"
	@echo "  experiment-threads - Test OpenMP thread scaling (1,2,4,8 threads)"
	@echo "  experiment-algorithms - Compare simple vs Smith-Waterman algorithms"
	@echo ""
	@echo "Experiment Output:"
	@echo "  All experiments save detailed results to text files"
	@echo "  Results include timing data and performance metrics"
	@echo "  Use these for performance analysis and reporting"
	@echo ""
	@echo "FASTA File Usage:"
	@echo "  make run-fasta REF=<reference.fasta> QUERIES=<queries.fasta> [options]"
	@echo ""
	@echo "Options for run-fasta:"
	@echo "  REF=<file> - Reference genome FASTA file (required)"
	@echo "  QUERIES=<file> - Query sequences FASTA file (required)"
	@echo "  ALGO=<algorithm> - Algorithm: 'simple' or 'sw' (default: simple)"
	@echo "  NP=<processes> - Number of MPI processes (default: 4)"
	@echo "  THREADS=<threads> - OpenMP threads per process (default: 2)"
	@echo ""
	@echo "Examples:"
	@echo "  make run-fasta REF=coronavirus.fasta QUERIES=coronavirus.fasta"
	@echo "  make run-fasta REF=coronavirus.fasta QUERIES=coronavirus.fasta ALGO=sw NP=8 THREADS=4"
	@echo ""
	@echo "Direct command line usage:"
	@echo "  mpirun -np <procs> ./$(TARGET) <reference.fasta> <queries.fasta> [algorithm]"
	@echo "  Example: mpirun -np 4 ./$(TARGET) coronavirus.fasta coronavirus.fasta sw"
	@echo ""
	@echo "Environment Variables:"
	@echo "  OMP_NUM_THREADS - Override OpenMP thread count"
	@echo "  OMP_SCHEDULE - Set OpenMP scheduling (e.g., 'dynamic,1')"

.PHONY: all clean run run-sw run-fasta test-performance experiment-scaling experiment-threads experiment-algorithms experiment-comprehensive help