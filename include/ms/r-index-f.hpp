/* r-index-f - Computes the simple r-index-f LF mapping table from RLE-BWT
    Copyright (C) 2020 Massimiliano Rossi

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file r-index-f.hpp
   \brief r-index-f.hpp Computes the R-Index-F table from RLE-BWT
   \author Massimiliano Rossi
   \author Nathaniel Brown
   \date 09/07/2020
*/

#ifndef _R_INDEX_F_HH
#define _R_INDEX_F_HH

#include <common.hpp>
#include <utility>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>

#include <r_index.hpp>

#include <ms_rle_string.hpp>

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd>
class r_index_f : ri::r_index<sparse_bv_type, rle_string_t>
{
public:
    // make private
    struct block_table_row 
    {
        char block_character;
        ulint block_num;
        ulint block_length;
        ulint block_offset;
    } F_block;

    typedef size_t size_type;
    vector<F_block> LF_table; 

    r_index_f() {}

    r_index_f(std::string filename) : ri::r_index<sparse_bv_type, rle_string_t>()
    {
        verbose("Building the simple r-index-f from BWT");

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        std::string bwt_fname = filename + ".bwt";

        verbose("RLE encoding BWT");

        std::string bwt_heads_fname = bwt_fname + ".heads";
        std::ifstream ifs_heads(bwt_heads_fname);
        std::string bwt_len_fname = bwt_fname + ".len";
        std::ifstream ifs_len(bwt_len_fname);
        this->bwt = rle_string_t(ifs_heads, ifs_len);
        this->r = this->bwt.number_of_runs();

        ifs_heads.seekg(0);
        ifs_len.seekg(0);
        this->build_F_(ifs_heads, ifs_len);
        // Can likely combine parsing when building F and block table
        this->build_block_table(ifs_heads, ifs_len)

        ri::ulint n = this->bwt.size();
        int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = ", this->r);
        verbose("Rate n/r = ", double(this->bwt.size()) / this->r);
        verbose("log2(r) = ", log2(double(this->r)));
        verbose("log2(n/r) = ", log2(double(this->bwt.size()) / this->r));

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Simple r-index-f construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    vector<ulint> build_F_(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->F = vector<ulint>(256, 0);
        int c;
        ulint i = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
                this->F[c] += length;
            else
            {
                this->F[TERMINATOR] += length;
                this->terminator_position = i;
            }
            i++;
        }
        for (ulint i = 255; i > 0; --i)
            this->F[i] = this->F[i - 1];
        this->F[0] = 0;
        for (ulint i = 1; i < 256; ++i)
            this->F[i] += this->F[i - 1];
        return this->F;
    }

    vector<F_block> build_LF_table(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);

        this->LF_table = vector<F_block>(this->r);
        vector<vector<int>> L_block_indices = vector<vector<int>>(256);
        
        int c;
        ulint i = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
            {
                this->LF_table[i].block_character = c;
                this->LF_table[i].block_length = length
                L_block_indices[c].push_back(i);
            }
            // check if we can simply assign terminator to C to clean up code
            else
            {
                this->LF_table[i].block_character = TERMINATOR;
                this->LF_table[i].block_length = length;
                L_block_indices[TERMINATOR].push_back(i);
            }
            i++;
        }
        
        ulint curr_L_num = 0;
        ulint L_seen = 0;
        ulint F_seen = 0;
        for(i = 0; i < L_block_indices.size(); i++) 
        {
            if (!L_block_indices[i].empty())
            {
                for(int j = 0; j < L_block_indices[i].size(); j++) 
                {
                    F_block* curr_block = &this->LF_table[L_block_indices[i][j]];

                    curr_block->block_num = curr_L_num
                    curr_block->block_offset = F_seen - L_seen;

                    F_seen += curr_block->block_length
                    
                    // Next blocks, sometimes skip
                    while (F_seen >= L_seen + this->LF_table[curr_L_num].block_length) {
                        curr_L_num++;
                        L_seen += this->LF_table[curr_L_num].block_length;
                    }
                }
            }
        }

        return this->LF_table;
    }

    void print_stats()
    {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("   terminator_position: ", sizeof(this->terminator_position));
        verbose("                     F: ", my_serialize(this->F, ns));
        verbose("                   bwt: ", this->bwt.serialize(ns));
    }

    /*
     * \param Block position (RLE blocks)
     * \param Current character offset in block
     * \return block position and offset of preceding character
     */
    std::pair<ulint, ulint> LF(ri::ulint block, ri::ulint offset)
    {
        ulint next_block = this->LF_table[block].block_num;
	    ulint next_offset = this->LF_table[block].block_offset + offset;

	    if (next_offset >= this->LF_table[next_block].block_length) {
            next_offset -= this->LF_table[next_block].block_length;
            next_block++;
        }
	
	    return std::make_pair(next_block, next_offset);
    }

    // Should serialize the table itself
    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += this->bwt.serialize(out);

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    // And load table
    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in)
    {
        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        my_load(this->F, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();
    }

#endif /* end of include guard: _R_INDEX_F_HH */