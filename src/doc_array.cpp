/*
 * File: doc_array.cpp
 * Description: Implementation of an document array ...
 *
 * Author: Omar Ahmed
 * Start Date: January 14, 2021
 */

#include <doc_array.hpp>
#include <spumoni_main.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <math.h> 
#include <sdsl/vectors.hpp>

DocumentArray::DocumentArray(std::string ref_path, size_t num_runs): ref_file(ref_path) {
    /* Constructs a document array - maps suffix array positions to genomes they occur in */

    // Verifies the file sizes are expected, and determine num of entries
    size_t start_samples_size = grab_file_size(ref_path + ".ssa");
    size_t end_samples_size = grab_file_size(ref_path + ".esa");

    ASSERT(start_samples_size ==  end_samples_size, "The .ssa and .esa files are not equally sized as expected.");
    this->num_entries = (start_samples_size)/ (2 * SSABYTES);
    
    // Read through the samples ...
    std::vector<size_t> start_samples_orig = read_samples(ref_path + ".ssa");
    std::vector<size_t> end_samples_orig = read_samples(ref_path + ".esa");

    // Determine the ending positions for each interval
    std::vector<size_t> end_pos;
    load_seq_boundaries(); // loads the seq_lengths vector
    end_pos.resize(this->seq_lengths.size());

    size_t i = 0, sum = 0;
    std::transform(seq_lengths.begin(), seq_lengths.end(), end_pos.begin(),
                  [&](size_t size) {if (i == 0){sum = seq_lengths[0]; i++; return sum;} 
                                    else {sum += seq_lengths[i];
                                          i++; return sum;}});
    
    size_t last_pos = end_pos.back();
    end_pos.back() = last_pos + 1; // add 1 for dollar sign

    // Convert suffix samples to the positions of BWT characters
    auto convert_to_bwt_pos = [&] (size_t sample) {
        if (sample > 0) {return (sample-1);}
        else {return (end_pos.back()-1);}
    };

    std::vector<size_t> start_samples, end_samples;
    start_samples.resize(start_samples_orig.size());
    end_samples.resize(end_samples_orig.size());

    #pragma omp parallel for
    for (size_t i = 0; i < start_samples.size(); i++) {
        start_samples[i] = convert_to_bwt_pos(start_samples_orig[i]);
        end_samples[i] = convert_to_bwt_pos(end_samples_orig[i]);
    }
    ASSERT((start_samples_orig.size() == end_samples_orig.size()), "issue occurred during"
    " the document array construction.");

    // Old code
    // std::transform(start_samples_orig.begin(), start_samples_orig.end(), 
    //                start_samples.begin(), convert_to_bwt_pos);
    // std::transform(end_samples_orig.begin(), end_samples_orig.end(), 
    //                end_samples.begin(), convert_to_bwt_pos);

    // Build a bitvector that marks the end of each document with a 1
    sdsl::bit_vector doc_ends = sdsl::bit_vector(end_pos.back(), 0);
    for (auto x: end_pos)
        doc_ends[x-1] = 1;
    sdsl::rank_support_v<1> doc_ends_rank = sdsl::rank_support_v<1> (&doc_ends);

    // Perform a rank query on each suffix array position to
    // convert it to a document number
    std::vector<size_t> start_genome_ids, end_genome_ids;
    start_genome_ids.resize(start_samples.size());
    end_genome_ids.resize(end_samples.size());

    // Old code, that used binary search to find document number
    // std::transform(start_samples.begin(), start_samples.end(), 
    //                start_genome_ids.begin(), [&](size_t pos) {return binary_search_for_pos(end_pos, pos);});
    // std::transform(end_samples.begin(), end_samples.end(), 
    //                end_genome_ids.begin(), [&](size_t pos) {return binary_search_for_pos(end_pos, pos);});  

    #pragma omp parallel for
    for (size_t i = 0; i < start_samples.size(); i++) {
        start_genome_ids[i] = doc_ends_rank(start_samples[i]);
        end_genome_ids[i] = doc_ends_rank(end_samples[i]);
    }
    ASSERT((start_samples.size() == end_samples.size()), "issue occurred during"
    " the document array construction.");

    // Write the document array to int vectors
    uint32_t max_width = std::ceil(std::log2((this->seq_lengths.size() + 0.0)));
    DBG_ONLY("%s %d", "Number of bits used per document entry:", max_width);

    this->start_runs_doc = sdsl::int_vector<> (num_runs, 0, max_width);
    this->end_runs_doc = sdsl::int_vector<> (num_runs, 0, max_width);

    size_t k = 0, j = 0;
    std::for_each(start_genome_ids.begin(), start_genome_ids.end(), 
                  [&](size_t genome_id){start_runs_doc[k] = genome_id; k++;});
    std::for_each(end_genome_ids.begin(), end_genome_ids.end(), 
                  [&](size_t genome_id){end_runs_doc[j] = genome_id; j++;});
}

void DocumentArray::load_seq_boundaries() {
    /* 
     *  Takes in a FASTA document index from RefBuilder Class, and stores the length of each
     *  sequence in order to build an interval tree.
     */

    std::ifstream index_file (ref_file + ".fdi", std::ifstream::in);
    std::string line;

    while(std::getline(index_file, line)) {
        auto word_list = split(line, '\t');
        std::string seq_length = word_list[1];

        bool is_num = std::all_of(seq_length.begin(), seq_length.end(), [](char c){return std::isdigit(c);});
        ASSERT(is_num, "Issue with FASTA index, sequence length is not a number.");
        this->seq_lengths.push_back(std::atol(seq_length.data()));
    }
}

size_t DocumentArray::grab_file_size(std::string file_path) {
    /* helper method that determines the size of a file */
    std::ifstream samples (file_path, std::ifstream::binary);
    if (!samples.is_open()) {FATAL_ERROR("A path to the suffix array samples is not valid.");}

    // Get size of file and determine number of entries ...
    samples.seekg (0, samples.end);
    size_t length_bytes = samples.tellg();
    samples.seekg (0, samples.beg);
    return length_bytes;
}

void DocumentArray::print_statistics() {
    /* Prints out attributes and document array */
    std::cout << this->num_entries << std::endl;
    for (size_t i = 0; i < this->num_entries; i++) {
        std::cout << start_runs_doc[i] << " " << end_runs_doc[i] << std::endl;
    }
}

std::vector<size_t> DocumentArray::read_samples(std::string file_path) {
    /* Reads in the suffix array samples from the provided file */
    std::ifstream samples_file (file_path, std::ifstream::binary);
    std::vector<size_t> samples;

    uint64_t left = 0, right = 0;
    size_t pos = 0;

    while (samples_file.read((char*) &left, SSABYTES) && samples_file.read((char*) &right, SSABYTES)) {
        samples.push_back(right);
    }
    samples_file.close();
    return samples;
}

size_t DocumentArray::binary_search_for_pos(std::vector<size_t> end_pos, size_t sample_pos) {
    /* Performs a binary search to determine what genome a certain offset occurs in */
    size_t low = 0, mid = 0, high = end_pos.size() - 1;
    size_t true_pos = 0;

    while (low < high) {
        mid = (low + high)/2;
        if (end_pos[mid] == sample_pos) {true_pos = mid+1; break;} // The positions are non-inclusive, that explains +1
        else if (sample_pos < end_pos[mid]){high = mid;}
        else {low = mid + 1;}
    }
    if (low >= high) {true_pos = high;} // only applies update when binary search was exhausted
    if (sample_pos >= end_pos[true_pos]) {true_pos++;} // Make adjustment since these positions are not inclusive

    std::string assert_msg = "binary search during document array building has an issue.";
    if (true_pos > 0) { ASSERT((sample_pos < end_pos[true_pos] && sample_pos >= end_pos[true_pos-1]), assert_msg.data());}
    else {ASSERT((sample_pos < end_pos[true_pos]), assert_msg.data());}
    return true_pos;
}

size_t DocumentArray::serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) {
    sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written_bytes = 0;

    out.write((char *)&this->num_entries, sizeof(this->num_entries));
    written_bytes += sizeof(this->num_entries);

    written_bytes += this->start_runs_doc.serialize(out, child, "start_doc");
    written_bytes += this->end_runs_doc.serialize(out, child, "end_doc");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void DocumentArray::load(std::istream& in) {
    in.read((char *)&this->num_entries, sizeof(this->num_entries));
    start_runs_doc.load(in);
    end_runs_doc.load(in);
}
