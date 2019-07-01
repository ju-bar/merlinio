// file : "merlin_hdr.h"
//
// Declares structures related to merlin file headers
//
/* -----------------------------------------------------------------------

This file is part of Merlinio.

	Merlinio is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	Merlinio is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with Merlinio.  If not, see <http://www.gnu.org/licenses/>.

----------------------------------------------------------------------- */
#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

constexpr auto MERLIN_FRAME_HDR_SIZE_MAX = 2048;


struct merlin_roi {
	int x0 = 0;
	int y0 = 0;
	int x1 = 0;
	int y1 = 0;
};

struct merlin_pix {
	int x = 0;
	int y = 0;
};

struct merlin_pos {
	double x = 0.;
	double y = 0.;
};

struct merlin_range {
	double min = 0.;
	double max = 0.;
};

struct merlin_hdr {
	int n_frames = 0; // nframes = nx*ny
	int n_columns = 0; // nx
	int n_rows = 0; // ny
	int n_files = 0; // number of files
	size_t n_fhdr_bytes = 0; // global frame header length in bytes (assuming similar headers)
	size_t n_data_bytes = 0; // global frame data length in bytes (assuming similar data)
	std::string s_timestamp;
};


struct merlin_frame_hdr {
	__int16 n_size = 0; // header size in bytes
	__int32 n_columns = 0; // nx
	__int32 n_rows = 0; // ny
	__int32 i_seq = 0; // frame acquisition sequence (zero based index)
	__int8 n_chips = 0; // number of chips
	__int8 n_bpi = 16; // number of bits per item
	__int8 n_chip_select = 0; // chip selection bits (least significant bit is first chip)
	double d_dwell = 0; // frame dwell time in seconds
	std::string s_sensor_layout; // sensor layout string
	std::string s_hid; // header id string
	std::string s_time; // frame time stamp
};

struct merlin_frame_calib {
	merlin_pos offset; // origin of the coordinate system
	merlin_pos a0; // first basis vector -> (xi',yi') = i * a0
	merlin_pos a1; // second basis vector -> (xj',yj') = j * a1
};

int imod(int i, int n);

// reads data from the merlin main header file opened as ifstream
// - fills information to the provided merlin_frame_hdr * phdr
int merlin_read_header(std::ifstream * pfin, merlin_hdr * phdr);

// gets the next parameter string from the frame header string starting
// at position ipos of the header string. ipos is expected to be
// the first character of the parameter string. The function returns
// the position of the next parameter string in the header string.
int merlin_read_frame_header_param(int ipos, std::string * pstr_hdr, std::string * prm);

// reads the frame header from the current file position
// - fills information to the provided merlin_frame_hdr * pfhdr
// - the stream position moved to after the header
int merlin_read_frame_header(std::ifstream * pfin, merlin_frame_hdr * pfhdr);

// reads data from file from given file position and returns as double
// (pre-processing can be applied in this function, switches via merlin_hdr)
// - output buf = pointer to buffer recieving pre-processed data
// - input pos = file position to read data from
// - input pfin = pointer ton an open input file stream
// - input phdr = pointer to struct merlin_hdr
// - input swapbytes = do byte swap on data transfer to double
// returns an error code > 0 in case of failures
int merlin_read_data(double *buf, std::streampos pos, std::ifstream * pfin, merlin_hdr * phdr, merlin_frame_hdr * pfhdr, bool swapbytes);


