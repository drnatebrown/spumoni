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
#include <iostream> 
#include <algorithm>
#include <random>

#include <malloc_count.h>

#include <sdsl/rmq_support.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/dac_vector.hpp>

#include <r_index.hpp>

#include <ms_rle_string.hpp>
#include <thresholds_ds.hpp>

static const int BLOCK_SIZE = 1024;

template <class sparse_bv_type = ri::sparse_sd_vector,
          class rle_string_t = ms_rle_string_sd,
          class thresholds_t = thr_compressed<rle_string_t> >
class r_index_f : ri::r_index<sparse_bv_type, rle_string_t>
{
public:
    typedef size_t size_type;
    thresholds_t thresholds;

    /*
    enum
    {
        BIT_A = 0x0;
        BIT_C = 0x1;
        BIT_G = 0x2;
        BIT_T = 0x3;
        BIT_TERMINAL = 0x4;
    }
    */

    struct i_block
    {
        wt_huff<rrr_vector<63>> heads;
        std::map<char, ulint> c_map;
        std::map<char, rrr_vector<63>> c_diff;
        dac_vector<> lengths;
        dac_vector<> offsets;
    };

    vector<i_block> B_table; 

    r_index_f() {}

    r_index_f(std::string filename) : ri::r_index<sparse_bv_type, rle_string_t>()
    {
        verbose("Building the R-Index-F table from RLE-BWT");

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
        //this->build_F_(ifs_heads, ifs_len);
        build_B_table(ifs_heads, ifs_len);

        ri::ulint n = this->bwt.size();
        int log_r = bitsize(uint64_t(this->r));
        int log_n = bitsize(uint64_t(this->bwt.size()));

        verbose("Number of BWT equal-letter runs: r = ", this->r);
        verbose("Rate n/r = ", double(this->bwt.size()) / this->r);
        verbose("log2(r) = ", log2(double(this->r)));
        verbose("log2(n/r) = ", log2(double(this->bwt.size()) / this->r));

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Block-Table construction complete");
        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

        verbose("Reading thresholds from file");

        t_insert_start = std::chrono::high_resolution_clock::now();

        thresholds = thresholds_t(filename,&this->bwt);

        t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("Memory peak: ", malloc_count_peak());
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    }

    vector<i_block> build_B_table(std::ifstream &heads, std::ifstream &lengths)
    {
        heads.clear();
        heads.seekg(0);
        lengths.clear();
        lengths.seekg(0);
        
        vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(256);
        vector<char> chars = vector<char>();
        vector<ulint> lens = vector<ulint>();
        
        char c;
        ulint i = 0;
        while ((c = heads.get()) != EOF)
        {
            size_t length = 0;
            lengths.read((char *)&length, 5);
            if (c > TERMINATOR)
            {
                chars.push_back(c);
                lens.push_back(length);
                L_block_indices[c].push_back(i);
            }
            else
            {
                chars.push_back(TERMINATOR);
                lens.push_back(length);
                L_block_indices[TERMINATOR].push_back(i);
            }

            /*
            ++i;
            ++b_i;

            // End of block of intervals: compress heads and lengths
            if (i % BLOCK_SIZE > b)
            {
                ++b;
                b_i = 0;

                i_block curr;
                curr.heads = wt_huff<rrr_vector<63>>(block_chars);
                curr.lengths = dac_vector(block_lens);

                block_chars = vector<char>(BLOCK_SIZE);
                block_lens = vector<ulint>(BLOCK_SIZE);
            }
            */
        }
        
        ulint r = chars.size();

        vector<ulint> intervals = vector<ulint>(r);
        vector<ulint> offsets = vector<ulint>(r);

        ulint curr_L_num = 0;
        ulint L_seen = 0;
        ulint F_seen = 0;
        for(size_t i = 0; i < L_block_indices.size(); ++i) 
        {
            for(size_t j = 0; j < L_block_indices[i].size(); ++j) 
            {
                //F_block* curr_block = &LF_table[L_block_indices[i][j]];
                ulint pos = L_block_indices[i][j];

                intervals[pos] = curr_L_num;
                offsets[pos] = F_seen - L_seen;

                F_seen += lens[pos];
            
                while (F_seen >= L_seen + lens[curr_L_num]) 
                {
                    L_seen += lens[curr_L_num];
                    ++curr_L_num;
                }
            }
        }

        ulint B_len = (r/BLOCK_SIZE) + ((r % BLOCK_SIZE) != 0);
        B_table = vector<i_block>(B_len);

        vector<char> block_chars = vector<char>(BLOCK_SIZE);
        vector<ulint> block_lens = vector<ulint>(BLOCK_SIZE);
        vector<ulint> block_offsets = vector<ulint>(BLOCK_SIZE);
        std::map<char, ulint> block_c_map = std::map<char, ulint>();
        std::map<char, vector<bool>> bit_diff = std::map<char, vector<bool>>();

        ulint b = 0;
        ulint b_i = 0;
        for (size_t i = 0; i < r; ++i) 
        {
            cout << i << "\n";

            // End of block of intervals, update block table
            if (i/BLOCK_SIZE > b || (i+1) == r)
            {
                cout << "YEAH" << "b";
                i_block& curr = B_table[b];

                curr.heads = wt_huff<rrr_vector<63>>();
                // TODO : Write heads to file (slow), or find another way
                //std::ofstream FILE("blockchar"+b, std::ios::out | std::ofstream::binary);
                //std::copy(block_chars.begin(), block_chars.end(), std::ostreambuf_iterator<char>(FILE));
                //construct(curr.heads, "blockchar"+b, 1);

                curr.lengths = dac_vector(block_lens);
                curr.offsets = dac_vector(block_offsets);
                curr.c_map = block_c_map;

                std::map<char, rrr_vector<63>> block_c_diff;
                for (auto& kv: bit_diff) 
                {
                    block_c_diff.insert(std::pair(kv.first, rrr(kv.second)));
                }
                curr.c_diff = block_c_diff;

                block_chars = vector<char>(BLOCK_SIZE);
                block_lens = vector<ulint>(BLOCK_SIZE);
                block_offsets = vector<ulint>(BLOCK_SIZE);
                block_c_map = std::map<char, ulint>();
                bit_diff = std::map<char, vector<bool>>();

                ++b;
                b_i = 0;
            }

            char c = chars[i];
            ulint l = lens[i];
            ulint d = offsets[i];

            block_chars[b_i] = c;
            block_lens[b_i] = l;
            block_offsets[b_i] = d;

            if (!block_c_map.count(c)) {
                block_c_map.insert(std::pair<char, ulint>(c, l));
                bit_diff.insert(std::pair<char, vector<bool>>(c, vector<bool>()));
            }

            ulint diff = l - block_c_map[c];
            while (diff > 0) {
                bit_diff[c].push_back(false);
                --diff;
            }
            bit_diff[c].push_back(true);
        }

        return B_table;
    }

    rrr_vector<63> rrr(vector<bool> &b){

		if(b.size()==0) return rrr_vector<63>();

		bit_vector bv(b.size());

		for(uint64_t i=0;i<b.size();++i)
			bv[i] = b[i];

		return rrr_vector<63>(bv);
	}

    void print_stats()
    {
        sdsl::nullstream ns;

        verbose("Memory consumption (bytes).");
        verbose("   terminator_position: ", sizeof(this->terminator_position));
        //verbose("                     F: ", my_serialize(this->F, ns));
        verbose("              Block table: ", my_serialize_vector_of_structs(B_table, ns));
        verbose("                   bwt: ", this->bwt.serialize(ns));
    }

    void bwt_stats()
    {
        verbose("Number of BWT equal-letter runs: r = ", this->r);
        verbose("Length of complete BWT: n = ", this->bwt.size());
        verbose("Rate n/r = ", double(this->bwt.size()) / this->r);
        verbose("log2(r) = ", log2(double(this->r)));
        verbose("log2(n/r) = ", log2(double(this->bwt.size()) / this->r));
    }

    // Lives here for now, can move into tests if we expose the LF Table
    /*
    void invert_bwt(std::string filename) 
    {
        verbose("Inverting BWT using R-Index-F (LF table)");
        ulint num = 0;

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        //vector<char> recovered = vector<char>();
        ulint block = 0;
        ulint offset = 0;

        char c;
        while((c = LF_table[block].character) > TERMINATOR) 
        {
            //recovered.push_back(char(c));
            std::pair<ulint, ulint> block_pair = LF(block, offset);
            block = block_pair.first;
            offset = block_pair.second;
        }

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        verbose("BWT Inverted using LF Table");
        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        verbose("Average step (ns): ", std::chrono::duration<double, std::ratio<1, 1000000000>>((t_insert_end - t_insert_start)/this->bwt.size()).count());
        verbose("# of runs, r: ",this->r);
        verbose("BWT size, n: ", this->bwt.size());

        /*
        std::ofstream recovered_output(filename + ".LF_recovered");

        std::reverse(recovered.begin(), recovered.end());
        std::string recovered_string = string(recovered.begin(), recovered.end());

        recovered_output << recovered_string;
        recovered_output.close();

        verbose("Recovered text written to", filename + ".LF_recovered");
        
    }
    

    void sample_LF(size_t samples, unsigned seed)
    {
        verbose("Running random sample of LF steps for R-Index-F (LF table):");

        std::mt19937_64 gen(seed);
        std::uniform_int_distribution<ulint> dist(0, this->bwt.size());
        vector<std::pair<ulint, ulint>> pos = vector<std::pair<ulint, ulint>>(samples);
        vector<std::pair<ulint, ulint>> next_pos = vector<std::pair<ulint, ulint>>(samples);
        
        for(size_t i = 0; i < pos.size(); ++i)
        {
            pos[i] = position_to_table(dist(gen));
        }

        std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

        for(size_t i = 0; i < pos.size(); ++i)
        {
            next_pos[i] = LF(pos[i].first, pos[i].second);
        }

        std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

        /*
        for(size_t i = 0; i < next_pos.size(); ++i)
        {
            ulint pos = this->bwt.run_range(next_pos[i].first).first + next_pos[i].second;
            cerr << pos << "\n";
        }
        

        verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
        verbose("Average step (ns): ", std::chrono::duration<double, std::ratio<1, 1000000000>>((t_insert_end - t_insert_start)/samples).count());
        verbose("# of samples: ", samples);
    }

    /*
     * \param Block position (RLE blocks)
     * \param Current character offset in block
     * \return block position and offset of preceding character
    std::pair<ulint, ulint> LF(ri::ulint block, ri::ulint offset)
    {
        ulint next_block = LF_table[block].block;
	    ulint next_offset = LF_table[block].offset + offset;

	    while (next_offset >= LF_table[next_block].length) 
        {
            next_offset -= LF_table[next_block].length;
            ++next_block;
        }

	    return std::make_pair(next_block, next_offset);
    }

    // Takes a position from the BWT and returns the block position and offset in the LF table
    std::pair<ulint, ulint> position_to_table(ulint i){
        assert(i < this->bwt.size());
        ulint block = this->bwt.run_of_position(i);
        assert(block < LF_table.size());
        // MAKE FASTER, FINE FOR NOW
        ulint offset = i - this->bwt.run_pos(block);
        return std::make_pair(block, offset);
    }

    // Takes a block and offset from the LF table and returns position in the BWT
    ulint table_to_position(ulint block, ulint offset)
    {
        assert(i < LF_table.size());
        ulint pos = this->bwt.run_pos(block) + offset;
        assert(pos < this->bwt.size());
        return pos;
    }
    */

    /* serialize the structure to the ostream
     * \param out     the ostream
     */
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") // const
    {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = 0;

        out.write((char *)&this->terminator_position, sizeof(this->terminator_position));
        written_bytes += sizeof(this->terminator_position);
        //written_bytes += my_serialize(this->F, out, child, "F");
        written_bytes += my_serialize_vector_of_structs(B_table, out, child, "B_table");
        written_bytes += this->bwt.serialize(out);

        written_bytes += thresholds.serialize(out, child, "thresholds");

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    std::string get_file_extension() const
    {
        return thresholds.get_file_extension() + ".rif";
    }

    /* load the structure from the istream
     * \param in the istream
     */
    void load(std::istream &in)
    {
        in.read((char *)&this->terminator_position, sizeof(this->terminator_position));
        //my_load(this->F, in);
        my_load_vector_of_structs(B_table, in);
        this->bwt.load(in);
        this->r = this->bwt.number_of_runs();

        thresholds.load(in,&this->bwt);
    }
};

#endif /* end of include guard: _R_INDEX_F_HH */