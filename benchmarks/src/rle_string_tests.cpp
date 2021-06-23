/* rle_string_tests.cpp - Benchmark tests for R-Index using RLE String
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
   \file rle_string_tests.cpp
   \brief rle_string_tests.cpp Benchmarking using RLE String
   \author Massimiliano Rossi
   \author Nathaniel Brown
   \date 03/03/2021
*/

#include <iostream>

#define VERBOSE
#define SAMPLES 100000000
#define SEED 23

#include <common.hpp>

#include <sdsl/io.hpp>

#include <spumoni.hpp>

#include <malloc_count.h>


int main(int argc, char *const argv[])
{
    Args args;
    parseArgs(argc, argv, args);

    verbose("Loading the R-Index from RLE-BWT");
    
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();

    ms_pointers<> ms;

    std::string filename_ms = args.filename + ms.get_file_extension();

    ifstream fs_ms(filename_ms);
    ms.load(fs_ms);
    fs_ms.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

    verbose("R-Index load complete");
    verbose("Memory peak: ", malloc_count_peak());
    verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
    
    ms.bwt_stats();
    ms.print_stats();

    ms.invert_bwt(args.filename);
    ms.sample_LF(SAMPLES, SEED);

    return 0;
}