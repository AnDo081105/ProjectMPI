#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <random>
#include <cstring>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <sstream>

using namespace std;
// FASTA parser function
struct FASTASequence {
    string header;
    string sequence;
};

vector<FASTASequence> parseFASTA(const string& filename) {
    vector<FASTASequence> sequences;
    ifstream file(filename);
    
    if (!file.is_open()) {
        cerr << "Error: Cannot open FASTA file: " << filename << endl;
        return sequences;
    }
    
    string line, currentSeq;
    string currentHeader;

    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!currentSeq.empty()) {
                sequences.push_back({currentHeader, currentSeq});
                currentSeq.clear();
            }
            // Extract header (skip '>' character)
            currentHeader = (line.length() > 1) ? line.substr(1) : "";
        } else {
            currentSeq += line;
        }
    }
    if (!currentSeq.empty()) {
        sequences.push_back({currentHeader, currentSeq});
    }
    return sequences;
}

// DNA sequence alignment result structure
struct AlignmentResult {
    int queryIndex;
    int position;
    int score;
    string query;
    string matchedSegment;
    
    AlignmentResult() : queryIndex(-1), position(-1), score(-1) {}
    AlignmentResult(int qi, int pos, int s, const string& q, const string& match) 
        : queryIndex(qi), position(pos), score(s), query(q), matchedSegment(match) {}
};

// Structure for MPI communication
struct QueryData {
    int queryIndex;
    int queryLength;
    char querySequence[1024]; // Fixed size for MPI communication
    char header[256]; // Fixed size for header
    
    QueryData() : queryIndex(-1), queryLength(0) {
        memset(querySequence, 0, sizeof(querySequence));
        memset(header, 0, sizeof(header));
    }

    QueryData(int idx, const string& query, const string& hdr) : queryIndex(idx), queryLength(query.length()) {
        memset(querySequence, 0, sizeof(querySequence));
        strncpy(querySequence, query.c_str(), min(query.length(), sizeof(querySequence) - 1));
        memset(header, 0, sizeof(header));
        strncpy(header, hdr.c_str(), min(hdr.length(), sizeof(header) - 1));
    }
};

// Structure for result communication
struct ResultData {
    int queryIndex;
    int position;
    int score;
    int queryLength;
    int matchLength;
    char querySequence[1024];
    char matchedSegment[1024];
    char header[256];

    ResultData() : queryIndex(-1), position(-1), score(-1), queryLength(0), matchLength(0) {
        memset(querySequence, 0, sizeof(querySequence));
        memset(matchedSegment, 0, sizeof(matchedSegment));
        memset(header, 0, sizeof(header));
    }
    
    ResultData(const AlignmentResult& result) 
        : queryIndex(result.queryIndex), position(result.position), score(result.score),
          queryLength(result.query.length()), matchLength(result.matchedSegment.length()) {
        memset(querySequence, 0, sizeof(querySequence));    
        memset(matchedSegment, 0, sizeof(matchedSegment));
        memset(header, 0, sizeof(header));
        strncpy(querySequence, result.query.c_str(), min(result.query.length(), sizeof(querySequence) - 1));
        strncpy(matchedSegment, result.matchedSegment.c_str(), min(result.matchedSegment.length(), sizeof(matchedSegment) - 1));
    }
    
    ResultData(const AlignmentResult& result, const string& queryHeader) 
        : queryIndex(result.queryIndex), position(result.position), score(result.score),
          queryLength(result.query.length()), matchLength(result.matchedSegment.length()) {
        memset(querySequence, 0, sizeof(querySequence));    
        memset(matchedSegment, 0, sizeof(matchedSegment));
        memset(header, 0, sizeof(header));
        strncpy(querySequence, result.query.c_str(), min(result.query.length(), sizeof(querySequence) - 1));
        strncpy(matchedSegment, result.matchedSegment.c_str(), min(result.matchedSegment.length(), sizeof(matchedSegment) - 1));
        strncpy(header, queryHeader.c_str(), min(queryHeader.length(), sizeof(header) - 1));
    }
    
    AlignmentResult toAlignmentResult() const {
        string query(querySequence, queryLength);
        string match(matchedSegment, matchLength);
        return AlignmentResult(queryIndex, position, score, query, match);
    }
};

// Reference DNA sequence (chromosome segment)
string REFERENCE_DNA;
bool loadReferenceDNA(const string& filename) {
    auto sequences = parseFASTA(filename);
    if (sequences.empty()) {
        cerr << "Error: No sequences found in reference FASTA file." << endl;
        return false;
    }
    
    REFERENCE_DNA = sequences[0].sequence; // Use the first sequence as reference
    
    if (REFERENCE_DNA.empty()) {
        cerr << "Error: Reference sequence is empty." << endl;
        return false;
    }
    
    return true;
}

// Load query sequences from FASTA file
vector<string> loadQuerySequences(const string& filename, vector<string>& headers) {
    vector<string> queries;
    headers.clear();
    
    auto sequences = parseFASTA(filename);
    if (sequences.empty()) {
        cerr << "Error: No sequences found in query FASTA file: " << filename << endl;
        return queries;
    }
    
    for (const auto& seq : sequences) {
        // Skip empty sequences
        if (!seq.sequence.empty()) {
            queries.push_back(seq.sequence);
            headers.push_back(seq.header.empty() ? "unnamed_sequence" : seq.header);
        } else {
            cerr << "Warning: Skipping empty sequence with header: " << seq.header << endl;
        }
    }
    
    return queries;
}

// Validate DNA sequence (check for valid bases)
bool validateDNASequence(const string& sequence) {
    for (char base : sequence) {
        if (base != 'A' && base != 'T' && base != 'C' && base != 'G' && 
            base != 'a' && base != 't' && base != 'c' && base != 'g' && 
            base != 'N' && base != 'n') {
            return false;
        }
    }
    return true;
}

// Clean and validate DNA sequence (convert to uppercase, replace invalid chars with N)
string cleanDNASequence(const string& sequence) {
    string cleaned;
    for (char base : sequence) {
        char upperBase = toupper(base);
        if (upperBase == 'A' || upperBase == 'T' || upperBase == 'C' || upperBase == 'G') {
            cleaned += upperBase;
        } else {
            cleaned += 'N'; // Replace invalid characters with ambiguous base
        }
    }
    return cleaned;
}

// Calculate alignment score using simplified algorithm with OpenMP optimization
// Returns the best alignment score for the query against reference
AlignmentResult alignSequence(const string& query, int queryIndex) {
    if (REFERENCE_DNA.empty() || query.empty()) {
        return AlignmentResult(queryIndex, -1, -1, query, "");
    }
    
    if (query.length() > REFERENCE_DNA.length()) {
        return AlignmentResult(queryIndex, -1, -1, query, "");
    }
    
    int bestScore = -1;
    int bestPosition = -1;
    string bestMatch;
    
    size_t numPositions = REFERENCE_DNA.length() - query.length() + 1;
    
    // Parallelize the sliding window search
    #pragma omp parallel
    {
        int localBestScore = -1;
        int localBestPosition = -1;
        string localBestMatch;
        
        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numPositions; i++) {
            // Safety check for substring operation
            if (i + query.length() > REFERENCE_DNA.length()) {
                continue;
            }
            
            string segment = REFERENCE_DNA.substr(i, query.length());
            
            // Calculate simple alignment score (number of matching bases)
            int score = 0;
            for (size_t j = 0; j < query.length(); j++) {
                if (query[j] == segment[j]) {
                    score++;
                }
            }
            
            if (score > localBestScore) {
                localBestScore = score;
                localBestPosition = i;
                localBestMatch = segment;
            }
        }
        
        // Combine results from all threads
        #pragma omp critical(best_update)
        {
            if (localBestScore > bestScore) {
                bestScore = localBestScore;
                bestPosition = localBestPosition;
                bestMatch = localBestMatch;
            }
        }
    }
    
    return AlignmentResult(queryIndex, bestPosition, bestScore, query, bestMatch);
}

// Smith-Waterman local alignment algorithm implementation with OpenMP optimization
AlignmentResult smithWatermanAlignment(const string &query, int queryIndex)
{
    if (REFERENCE_DNA.empty() || query.empty())
    {
        return AlignmentResult(queryIndex, -1, -1, query, "");
    }

    const int MATCH_SCORE = 2;
    const int MISMATCH_SCORE = -1;
    const int GAP_PENALTY = -1;

    int queryLen = query.length();
    int refLen = REFERENCE_DNA.length();

    // Initialize scoring matrix
    vector<vector<int>> scoreMatrix(queryLen + 1, vector<int>(refLen + 1, 0));

    int maxScore = 0;
    int maxI = 0, maxJ = 0;

    // Fill the scoring matrix using Smith-Waterman algorithm
    // We can parallelize along anti-diagonals since they have no dependencies
    for (int diag = 1; diag <= queryLen + refLen; diag++)
    {
        int startI = max(1, diag - refLen);
        int endI = min(queryLen, diag - 1);

        if (startI <= endI)
        {
#pragma omp parallel for schedule(static) reduction(max : maxScore)
            for (int i = startI; i <= endI; i++)
            {
                int j = diag - i;
                if (j >= 1 && j <= refLen)
                {
                    // Safety bounds checking
                    if (i - 1 < 0 || j - 1 < 0 || i - 1 >= queryLen || j - 1 >= refLen)
                    {
                        continue;
                    }

                    // Calculate scores for three possible moves
                    int matchScore = scoreMatrix[i - 1][j - 1] +
                                     ((query[i - 1] == REFERENCE_DNA[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE);
                    int deleteScore = scoreMatrix[i - 1][j] + GAP_PENALTY;
                    int insertScore = scoreMatrix[i][j - 1] + GAP_PENALTY;

                    // Take maximum of the three scores, or 0 (local alignment)
                    scoreMatrix[i][j] = max({0, matchScore, deleteScore, insertScore});

                    // Track the maximum score and its position (thread-safe with reduction)
                    if (scoreMatrix[i][j] > maxScore)
                    {
#pragma omp critical(max_update)
                        {
                            if (scoreMatrix[i][j] > maxScore)
                            {
                                maxScore = scoreMatrix[i][j];
                                maxI = i;
                                maxJ = j;
                            }
                        }
                    }
                }
            }
        }
    }

    // Traceback to find the optimal alignment
    string alignedQuery = "";
    string alignedRef = "";
    int i = maxI, j = maxJ;

    while (i > 0 && j > 0 && scoreMatrix[i][j] > 0)
    {
        // Safety bounds checking for traceback
        if (i - 1 < 0 || j - 1 < 0 || i - 1 >= queryLen || j - 1 >= refLen)
        {
            break;
        }

        int currentScore = scoreMatrix[i][j];
        int matchScore = scoreMatrix[i - 1][j - 1] +
                         ((query[i - 1] == REFERENCE_DNA[j - 1]) ? MATCH_SCORE : MISMATCH_SCORE);
        int deleteScore = scoreMatrix[i - 1][j] + GAP_PENALTY;
        int insertScore = scoreMatrix[i][j - 1] + GAP_PENALTY;

        if (currentScore == matchScore)
        {
            // Match or mismatch
            if (i - 1 >= 0 && i - 1 < queryLen && j - 1 >= 0 && j - 1 < refLen)
            {
                alignedQuery = query[i - 1] + alignedQuery;
                alignedRef = REFERENCE_DNA[j - 1] + alignedRef;
            }
            i--;
            j--;
        }
        else if (currentScore == deleteScore)
        {
            // Deletion in reference (gap in query)
            if (j - 1 >= 0 && j - 1 < refLen)
            {
                alignedQuery = "-" + alignedQuery;
                alignedRef = REFERENCE_DNA[j - 1] + alignedRef;
            }
            i--;
        }
        else if (currentScore == insertScore)
        {
            // Insertion in reference (gap in reference)
            if (i - 1 >= 0 && i - 1 < queryLen)
            {
                alignedQuery = query[i - 1] + alignedQuery;
                alignedRef = "-" + alignedRef;
            }
            j--;
        }
        else
        {
            // Fallback case
            break;
        }
    }

    // Calculate the starting position in reference (accounting for traceback)
    int refStartPos = j;

    // Return result with Smith-Waterman score and aligned segment
    return AlignmentResult(queryIndex, refStartPos, maxScore, query, alignedRef);
}

// Enhanced alignment function that can choose between algorithms
AlignmentResult performAlignment(const string& query, int queryIndex, bool useSmithWaterman = false) {
    if (useSmithWaterman) {
        return smithWatermanAlignment(query, queryIndex);
    } else {
        return alignSequence(query, queryIndex);
    }
}

// Master process function using Scatterv/Gatherv
void masterProcess(int numProcesses, const vector<string>& queries, const vector<string>& queryHeaders, bool useSmithWaterman) {
    int numQueries = queries.size();
    
    cout << "\n=== Hybrid MPI+OpenMP DNA Sequence Alignment ===" << endl;
    cout << "Algorithm: " << (useSmithWaterman ? "Smith-Waterman Local Alignment" : "Simple Sliding Window") << endl;
    cout << "Parallelization: MPI (distributed) + OpenMP (shared memory)" << endl;
    cout << "Reference DNA length: " << REFERENCE_DNA.length() << " bases" << endl;
    cout << "Number of query sequences: " << numQueries << endl;
    cout << "MPI processes: " << numProcesses << endl;
    cout << "OpenMP threads per process: " << omp_get_max_threads() << endl;
    cout << "Total parallel workers: " << numProcesses * omp_get_max_threads() << endl;
    cout << "Distributing tasks using MPI_Scatterv..." << endl;
    
    auto startTime = chrono::high_resolution_clock::now();      
    
    // Prepare query data for scattering
    vector<QueryData> allQueryData(numQueries);
    for (int i = 0; i < numQueries; i++) {
        string header = (static_cast<size_t>(i) < queryHeaders.size()) ? queryHeaders[i] : "Query_" + to_string(i);
        allQueryData[i] = QueryData(i, queries[i], header);
    }
    
    // Calculate send counts and displacements for each process
    vector<int> sendCounts(numProcesses, 0);
    vector<int> displacements(numProcesses, 0);
    
    // Load balancing: Sort queries by computational complexity and distribute evenly
    vector<pair<long long, int>> queryComplexity;
    for (int i = 0; i < numQueries; i++) {
        long long complexity = (long long)queries[i].length() * REFERENCE_DNA.length();
        queryComplexity.push_back({complexity, i});
    }
    
    // Sort by complexity (largest first)
    sort(queryComplexity.begin(), queryComplexity.end(), greater<pair<long long, int>>());
    
    // Track computational load per process
    vector<long long> processLoad(numProcesses, 0);
    vector<vector<int>> processQueries(numProcesses);
    
    // Assign each query to the process with the least current load
    for (const auto& entry : queryComplexity) {
        long long complexity = entry.first;
        int queryIdx = entry.second;
        
        // Find process with minimum load
        int minLoadProcess = 0;
        for (int p = 1; p < numProcesses; p++) {
            if (processLoad[p] < processLoad[minLoadProcess]) {
                minLoadProcess = p;
            }
        }
        
        // Assign query to this process
        processQueries[minLoadProcess].push_back(queryIdx);
        processLoad[minLoadProcess] += complexity;
        sendCounts[minLoadProcess]++;
    }
    
    // Reorder allQueryData to match the load-balanced distribution
    vector<QueryData> balancedQueryData;
    for (int p = 0; p < numProcesses; p++) {
        for (int queryIdx : processQueries[p]) {
            string header = (static_cast<size_t>(queryIdx) < queryHeaders.size()) ? queryHeaders[queryIdx] : "Query_" + to_string(queryIdx);
            balancedQueryData.push_back(QueryData(queryIdx, queries[queryIdx], header));
        }
    }
    allQueryData = balancedQueryData;
    
    // Broadcast the send counts to all processes for load balanced distribution
    MPI_Bcast(sendCounts.data(), numProcesses, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Calculate displacements
    for (int i = 1; i < numProcesses; i++) {
        displacements[i] = displacements[i-1] + sendCounts[i-1];
    }
    
    // Convert to byte counts for MPI communication
    vector<int> sendCountsBytes(numProcesses);
    vector<int> displacementsBytes(numProcesses);
    for (int i = 0; i < numProcesses; i++) {
        sendCountsBytes[i] = sendCounts[i] * sizeof(QueryData);
        displacementsBytes[i] = displacements[i] * sizeof(QueryData);
    }
    
    // Display task distribution
    for (int i = 0; i < numProcesses; i++) {
        if (i == 0) {
            cout << "Master (rank 0) will process " << sendCounts[i] << " queries" << endl;
        } else {
            cout << "Process " << i << " will process " << sendCounts[i] << " queries" << endl;
        }
    }
    
    // Receive buffer for this process's queries
    vector<QueryData> myQueries(sendCounts[0]);
    
    // Scatter queries to all processes
    MPI_Scatterv(allQueryData.data(), sendCountsBytes.data(), displacementsBytes.data(), 
                 MPI_BYTE, myQueries.data(), sendCountsBytes[0], 
                 MPI_BYTE, 0, MPI_COMM_WORLD);
    
    // Process master's own queries using OpenMP parallelization
    vector<ResultData> myResults(myQueries.size());
    
    cout << "Master: About to process " << myQueries.size() << " queries with " 
         << omp_get_max_threads() << " OpenMP threads" << endl;
    cout.flush();
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < myQueries.size(); i++) {
        const auto& queryData = myQueries[i];
        string query(queryData.querySequence, queryData.queryLength);
        string header(queryData.header);
        
        AlignmentResult result = performAlignment(query, queryData.queryIndex, useSmithWaterman);
        myResults[i] = ResultData(result, header);
        
        // Thread-safe output
        #pragma omp critical(output)
        {
            string headerDisplay = header.empty() ? "No_Header" : 
                                 (header.length() > 20 ? header.substr(0, 20) + "..." : header);
            cout << "Master Thread " << omp_get_thread_num() << " processed query " 
                 << queryData.queryIndex << " (" << headerDisplay << ") "
                 << "(score: " << result.score << "/" << query.length() << ")" << endl;
            cout.flush();
        }
    }
    
    cout << "Master: Completed processing all queries" << endl;
    cout.flush();
    
    // Prepare for gathering results - convert to byte counts
    vector<int> recvCountsBytes(numProcesses);
    vector<int> recvDisplacementsBytes(numProcesses, 0);
    
    for (int i = 0; i < numProcesses; i++) {
        recvCountsBytes[i] = sendCounts[i] * sizeof(ResultData);
        if (i > 0) {
            recvDisplacementsBytes[i] = recvDisplacementsBytes[i-1] + recvCountsBytes[i-1];
        }
    }
    
    // Buffer to receive all results
    vector<ResultData> allResults(numQueries);
    
    // Gather results from all processes
    MPI_Gatherv(myResults.data(), myResults.size() * sizeof(ResultData), MPI_BYTE,
                allResults.data(), recvCountsBytes.data(), recvDisplacementsBytes.data(), 
                MPI_BYTE, 0, MPI_COMM_WORLD);
    
    auto endTime = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(endTime - startTime);
    
    // Filter and sort results directly using ResultData to preserve headers
    vector<ResultData> finalResults;
    for (const auto& resultData : allResults) {
        if (resultData.queryIndex >= 0) { // Valid result
            finalResults.push_back(resultData);
        }
    }
    
    sort(finalResults.begin(), finalResults.end(), 
         [](const ResultData& a, const ResultData& b) {
             return a.queryIndex < b.queryIndex;
         });
    
    // Display results with FASTA headers
    cout << "\n=== Alignment Results ===" << endl;
    cout << setw(5) << "Query" << setw(25) << "Header" << setw(10) << "Position" << setw(8) << "Score" 
         << setw(20) << "Query Sequence" << setw(20) << "Best Match" << endl;
    cout << string(88, '-') << endl;
    
    for (const auto& result : finalResults) {
        double similarity = (double)result.score / result.queryLength * 100;
        string header(result.header);
        string query(result.querySequence, result.queryLength);
        string match(result.matchedSegment, result.matchLength);
        
        // Safe string truncation
        string headerDisplay = header.empty() ? "No_Header" : 
                             (header.length() > 22 ? header.substr(0, 22) + "..." : header);
        string queryDisplay = query.empty() ? "Empty" :
                            (query.length() > 15 ? query.substr(0, 15) + "..." : query);
        string matchDisplay = match.empty() ? "Empty" :
                            (match.length() > 15 ? match.substr(0, 15) + "..." : match);
        
        cout << setw(5) << result.queryIndex 
             << setw(25) << headerDisplay
             << setw(10) << result.position 
             << setw(8) << result.score << "/" << result.queryLength
             << setw(20) << queryDisplay
             << setw(20) << matchDisplay
             << "  (" << fixed << setprecision(1) << similarity << "%)" << endl;
    }
    
    cout << "\nTotal execution time: " << duration.count() << " ms" << endl;
    cout << "Algorithm used: " << (useSmithWaterman ? "Smith-Waterman" : "Simple Sliding Window") << endl;
    cout << "Threading model: " << numProcesses << " MPI processes × " 
         << omp_get_max_threads() << " OpenMP threads = " 
         << numProcesses * omp_get_max_threads() << " total workers" << endl;
}

// Worker process function using Scatterv/Gatherv
void workerProcess(int rank, int numProcesses, int totalQueries, bool useSmithWaterman) {
    // Workers will receive the send counts from master (calculated with load balancing)
    vector<int> sendCounts(numProcesses);
    MPI_Bcast(sendCounts.data(), numProcesses, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Convert to byte counts for MPI communication
    vector<int> sendCountsBytes(numProcesses);
    for (int i = 0; i < numProcesses; i++) {
        sendCountsBytes[i] = sendCounts[i] * sizeof(QueryData);
    }
    
    // Receive buffer for this process's queries
    vector<QueryData> myQueries(sendCounts[rank]);
    
    // Participate in scatter operation
    MPI_Scatterv(nullptr, nullptr, nullptr, MPI_BYTE, 
                 myQueries.data(), sendCountsBytes[rank], 
                 MPI_BYTE, 0, MPI_COMM_WORLD);
    
    // Process assigned queries using OpenMP parallelization
    vector<ResultData> myResults(myQueries.size());
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < myQueries.size(); i++) {
        const auto& queryData = myQueries[i];
        string query(queryData.querySequence, queryData.queryLength);
        string header(queryData.header);
        AlignmentResult result = performAlignment(query, queryData.queryIndex, useSmithWaterman);
        myResults[i] = ResultData(result, header);
        
        // Thread-safe output
        #pragma omp critical(output)
        {
            string headerDisplay = header.empty() ? "No_Header" : 
                                 (header.length() > 20 ? header.substr(0, 20) + "..." : header);
            cout << "Process " << rank << " Thread " << omp_get_thread_num() 
                 << " processed query " << queryData.queryIndex 
                 << " (" << headerDisplay << ") "
                 << " (score: " << result.score << "/" << query.length() << ")" << endl;
        }
    }
    
    // Prepare receive counts for gathering (convert to byte counts)
    vector<int> recvCountsBytes(numProcesses);
    for (int i = 0; i < numProcesses; i++) {
        recvCountsBytes[i] = sendCounts[i] * sizeof(ResultData);
    }
    
    // Participate in gather operation
    MPI_Gatherv(myResults.data(), myResults.size() * sizeof(ResultData), MPI_BYTE,
                nullptr, recvCountsBytes.data(), nullptr, MPI_BYTE, 0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[]) {
    // Initialize MPI with threading support
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    
    if (provided < MPI_THREAD_FUNNELED) {
        cerr << "Warning: MPI implementation does not support threading" << endl;
    }
    
    int rank, numProcesses;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    
    // Set OpenMP thread count (can be overridden by OMP_NUM_THREADS environment variable)
    int numThreads = omp_get_max_threads();
    if (getenv("OMP_NUM_THREADS") == nullptr) {
        numThreads = min(4, (int)omp_get_num_procs()); // Default to 4 threads or number of cores
        omp_set_num_threads(numThreads);
    }
    
    if (rank == 0) {
        cout << "=== Hybrid MPI+OpenMP DNA Sequence Alignment ===" << endl;
        cout << "MPI Processes: " << numProcesses << endl;
        cout << "OpenMP Threads per process: " << numThreads << endl;
        cout << "Total parallel workers: " << numProcesses * numThreads << endl;
        cout << "MPI Threading support: " << 
                (provided >= MPI_THREAD_FUNNELED ? "Available" : "Limited") << endl << endl;
    }
    
    if (numProcesses < 2) {
        if (rank == 0) {
            cerr << "Error: Need at least 2 MPI processes for meaningful parallelization" << endl;
            cerr << "Usage: mpirun -np <procs> ./mpi_dna_alignment <reference.fasta> <queries.fasta> [algorithm]" << endl;
            cerr << "  reference.fasta: Reference genome file" << endl;
            cerr << "  queries.fasta: Query sequences file" << endl;
            cerr << "  algorithm: 'simple' or 'sw' (default: simple)" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    // New command line format:
    // ./mpi_dna_alignment <reference.fasta> <queries.fasta> [algorithm]
    
    if (argc < 3) {
        if (rank == 0) {
            cerr << "Usage: mpirun -np <procs> ./mpi_dna_alignment <reference.fasta> <queries.fasta> [algorithm]" << endl;
            cerr << "  reference.fasta: Reference genome file" << endl;
            cerr << "  queries.fasta: Query sequences file" << endl;
            cerr << "  algorithm: 'simple' or 'sw' (default: simple)" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    string refFile = argv[1];
    string queryFile = argv[2];
    string algorithm = (argc > 3) ? argv[3] : "simple";
    
    // Check file existence on master process
    if (rank == 0) {
        ifstream refTest(refFile);
        if (!refTest.good()) {
            cerr << "Error: Cannot access reference file: " << refFile << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        refTest.close();
        
        ifstream queryTest(queryFile);
        if (!queryTest.good()) {
            cerr << "Error: Cannot access query file: " << queryFile << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        queryTest.close();
    }
    
    // Parse algorithm selection
    bool useSmithWaterman = false;
    if (algorithm == "sw" || algorithm == "smith-waterman") {
        useSmithWaterman = true;
    }
    
    // Broadcast algorithm selection to all processes
    MPI_Bcast(&useSmithWaterman, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    
    // Load FASTA files and process sequences
    int numQueries = 0;
    vector<string> queries;
    vector<string> queryHeaders;
    
    if (rank == 0) {
        // Master process loads FASTA files
        cout << "\n=== Loading FASTA Files ===" << endl;
        cout << "Reference file: " << refFile << endl;
        cout << "Query file: " << queryFile << endl;
        cout << "Algorithm: " << algorithm << endl;
        
        // Load reference genome
        if (!loadReferenceDNA(refFile)) {
            cerr << "Error: Failed to load reference DNA from " << refFile << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        cout << "Reference DNA loaded: " << REFERENCE_DNA.length() << " bases" << endl;
        
        // Clean reference DNA sequence
        REFERENCE_DNA = cleanDNASequence(REFERENCE_DNA);
    }
    
    // Broadcast reference DNA length to all processes
    int refLength = REFERENCE_DNA.length();
    MPI_Bcast(&refLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Worker processes need to receive the reference DNA
    if (rank != 0) {
        REFERENCE_DNA.resize(refLength);
    }
    
    // Broadcast reference DNA to all processes
    MPI_Bcast(const_cast<char*>(REFERENCE_DNA.data()), refLength, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    // Verify reference DNA was received correctly by all processes
    if (REFERENCE_DNA.empty()) {
        cerr << "Error: Process " << rank << " - Reference DNA is empty after broadcast!" << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    if (rank == 0) {
        
        // Load query sequences
        queries = loadQuerySequences(queryFile, queryHeaders);
        if (queries.empty()) {
            cerr << "Error: Failed to load query sequences from " << queryFile << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        // Clean and validate query sequences
        for (size_t i = 0; i < queries.size(); i++) {
            queries[i] = cleanDNASequence(queries[i]);
            if (queries[i].empty()) {
                cerr << "Warning: Empty sequence found at index " << i << ", skipping." << endl;
                queries.erase(queries.begin() + i);
                queryHeaders.erase(queryHeaders.begin() + i);
                i--; // Adjust index after removal
            }
        }
        
        numQueries = queries.size();
        cout << "Query sequences loaded: " << numQueries << " sequences" << endl;
        
        if (numQueries == 0) {
            cerr << "Error: No valid query sequences found" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        
        // Display sequence statistics
        size_t totalBases = 0, minLen = SIZE_MAX, maxLen = 0;
        for (const auto& query : queries) {
            totalBases += query.length();
            minLen = min(minLen, query.length());
            maxLen = max(maxLen, query.length());
        }
        cout << "Sequence statistics:" << endl;
        cout << "  Total bases: " << totalBases << endl;
        cout << "  Average length: " << (double)totalBases / numQueries << " bases" << endl;
        cout << "  Length range: " << minLen << " - " << maxLen << " bases" << endl;
    }
    
    // Broadcast the number of queries to all processes
    MPI_Bcast(&numQueries, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        // Master process
        masterProcess(numProcesses, queries, queryHeaders, useSmithWaterman);
    } else {
        // Worker process
        workerProcess(rank, numProcesses, numQueries, useSmithWaterman);
    }
    
    MPI_Finalize();
    return 0;
}