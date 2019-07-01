// file : "merlin_hdr.cpp"
//
// Implementation of functions related to merlin file headers
//
#include "pch.h"
#include <vector>
#include "merlin_hdr.h"

int imod(int i, int n)
{
	return ((i%n + n) % n);
}


int merlin_read_header(std::ifstream * pfin, merlin_hdr * phdr)
{
	std::string str_line;
	std::vector<std::string> v_str_header;
	std::vector<std::string>::iterator i_hdr_ln;
	if ( NULL == pfin ) {
		return 1; // failure due to unknown parameter addresses
	}
	// check open file
	if (pfin->is_open()) {
		// read all header lines
		v_str_header.clear();
		while (getline(*pfin, str_line))
		{
			if (0 == str_line.find("End")) {
				v_str_header.push_back("End");
				break;
			}
			v_str_header.push_back(str_line);
		}
		pfin->close();
	}
	else {
		return 2; // failure due invalid stream state
	}
	if (phdr != NULL) {
		// parse the header lines for important information
		i_hdr_ln = v_str_header.begin();
		for (i_hdr_ln = v_str_header.begin(); i_hdr_ln != v_str_header.end(); i_hdr_ln++) {
			if (0 == i_hdr_ln->find("Time and Date Stamp (yr, mnth, day, hr, min, s):")) {
				phdr->s_timestamp = i_hdr_ln->substr(49, i_hdr_ln->size() - 49);
				continue;
			}
			if (0 == i_hdr_ln->find("Frames in Acquisition (Number):")) {
				phdr->n_frames = atoi(i_hdr_ln->substr(32, i_hdr_ln->size() - 32).c_str());
				continue;
			}
			if (0 == i_hdr_ln->find("Frames per Trigger (Number):")) {
				phdr->n_columns = atoi(i_hdr_ln->substr(29, i_hdr_ln->size() - 29).c_str());
				continue;
			}
		}
		if (phdr->n_frames > 0 && phdr->n_frames > phdr->n_columns && phdr->n_columns > 0) {
			phdr->n_rows = (phdr->n_frames - imod(phdr->n_frames, phdr->n_columns)) / phdr->n_columns;
		}
		if (0 < imod(phdr->n_frames, phdr->n_columns)) {
			phdr->n_rows++;
		}
	}
	return 0;
}

int merlin_read_frame_header_param(int ipos, std::string * pstr_hdr, std::string * prm)
{
	int lipos = (ipos>=0?ipos:0);
	int lcpos = lipos;
	int lmpos = 0;
	if (NULL == pstr_hdr) {
		return -1;
	}
	if (NULL == prm) {
		return -2;
	}
	lmpos = (int)pstr_hdr->size();
	while (pstr_hdr->at(lcpos) != ',' && lcpos < lmpos) {
		lcpos++;
	}
	if (lcpos - lipos > 1) {
		*prm = pstr_hdr->substr(lipos, lcpos - lipos);
		lcpos++;
	}
	// continue browsing the string until no more separators occur
	while (pstr_hdr->at(lcpos) == ',' && lcpos < lmpos) {
		lcpos++;
	}
	return lcpos; // return position which the possible begin of a new parameter or eos
}

int merlin_read_frame_header(std::ifstream * pfin, merlin_frame_hdr * pfhdr)
{
	std::streampos ipos = 0;
	char * cbuf = new char[MERLIN_FRAME_HDR_SIZE_MAX]; // i/o buffer
	std::string str_tmp;
	std::string str_num;
	int ihpos = 0, nhdr = 0;
	if (NULL == pfin) {
		return 1; // failure due to unknown parameter addresses
	}
	if (!pfin->is_open()) {
		return 2; // failure due to closed input steam
	}
	memset(cbuf, 0, MERLIN_FRAME_HDR_SIZE_MAX); // zero the buffer
	ipos = pfin->tellg(); // get current stream position
	if (pfin->fail()) { return 12; } // file position request failed
	// read the first part of the frame header (128 characters)
	pfin->read(cbuf, 128);
	if (pfin->fail()) { return 13; } // header test reading failed
	str_tmp = cbuf;
	// parse to 3rd entry for reading the header length
	// - 1st item = header ID (string)
	ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
	if (ihpos < 0 || ihpos >= 128) { return 3; } // pre-parsing error
	if (pfhdr) pfhdr->s_hid = str_num;
	// - 2nd item = acqusition sequence number (U32)
	ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
	if (ihpos < 0 || ihpos >= 128) { return 3; } // pre-parsing error
	if (pfhdr) pfhdr->i_seq = (__int32)(atoi(str_num.c_str())-1); // reduce by one to get the 0 based index (merlin starts with 1)
	// - 3rd item = header length in bytes (U16)
	ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
	if (ihpos < 0 || ihpos >= 128) { return 3; } // pre-parsing error
	nhdr = atoi(str_num.c_str());
	if (pfhdr) pfhdr->n_size = (__int16)nhdr;
	//
	// read the full header
	if (nhdr > 0 && nhdr <= MERLIN_FRAME_HDR_SIZE_MAX) {
		pfin->seekg(ipos);
		if (pfin->fail()) { return 14; } // file re-positioning failed
		pfin->read(cbuf, nhdr);
		if (pfin->fail()) { return 15; } // full header reading failed
		str_tmp.clear();
		str_tmp = cbuf;
	}
	else {
		return 4; // unsupported header size
	}
	//
	if (pfhdr != NULL) { // parse header content, continuing from ihpos
		// - 4th item = number of chips (U8)
		ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
		if (ihpos < 0 || ihpos >= nhdr) { return 5; } // parsing error
		pfhdr->n_chips = (__int8)atoi(str_num.c_str());
		// - 5th item = pixel dimension X (U32)
		ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
		if (ihpos < 0 || ihpos >= nhdr) { return 5; } // parsing error
		pfhdr->n_columns = (__int32)atoi(str_num.c_str());
		// - 6th item = pixel dimension Y (U32)
		ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
		if (ihpos < 0 || ihpos >= nhdr) { return 5; } // parsing error
		pfhdr->n_rows = (__int32)atoi(str_num.c_str());
		// - 7th item = pixel depth in file (string)
		ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
		if (ihpos < 0 || ihpos >= nhdr) { return 5; } // parsing error
		pfhdr->n_bpi = (__int8)atoi(str_num.substr(1,2).c_str());
		// - 8th item = sensor layout (string)
		ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
		if (ihpos < 0 || ihpos >= nhdr) { return 5; } // parsing error
		pfhdr->s_sensor_layout = str_num;
		// - 9th item = chip select (U8h)
		ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
		if (ihpos < 0 || ihpos >= nhdr) { return 5; } // parsing error
		pfhdr->n_chip_select = (__int8)std::strtoul(str_num.c_str(), 0, 16);
		// - 10th item = timestamp (string)
		ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
		if (ihpos < 0 || ihpos >= nhdr) { return 5; } // parsing error
		pfhdr->s_time = str_num;
		// - 11th item = shutter open time (double)
		ihpos = merlin_read_frame_header_param(ihpos, &str_tmp, &str_num);
		if (ihpos < 0 || ihpos >= nhdr) { return 5; } // parsing error
		pfhdr->d_dwell = (double)atof(str_num.c_str());
		// remaining items to be implemented and currently ignored,
		// essentially because we never use it and due to a messy
		// format specification.
		// - 12th item = counter (U8) for colour mode
		// - 13th item = colour mode (U8)
		// - 14th item = gain mode (U8)
		// - 15th-22nd item = 8 threashold values in keV (float)
		// - 23rd item ++ = DAC section (many numbers)
		// + extensions
		// + padding
	}
	pfin->seekg(ipos + (std::streampos)nhdr); // for some reasons the positioning after reading can fail, catch it!
	if (pfin->fail()) { return 16; } // file re-positioning failed
	return 0;
}

int merlin_read_data(double *buf, std::streampos pos, std::ifstream * pfin, merlin_hdr * phdr, merlin_frame_hdr * pfhdr, bool swapbytes)
{
	int nerr = 0;
	char * inbuf = NULL;
	char bs[4];
	unsigned __int8 * pdata8 = NULL;
	unsigned __int16 * pdata16 = NULL;
	unsigned __int32 * pdata32 = NULL;
	size_t i, j, n;
	if (NULL == buf) {
		return 1; // invalid input parameter 1
	}
	if (NULL == pfin) {
		return 3; // invalid input parameter 3
	}
	if (NULL == phdr) {
		return 4; // invalid input parameter 4
	}
	if (NULL == pfhdr) {
		return 5; // invalid input parameter 5
	}

	if (!pfin->is_open()) {
		return 13; // input file stream is not open for reading
	}
	
	pfin->seekg(pos);
	if (pfin->fail()) {
		return 23; // seeking pos ind input file stream failed
	}

	inbuf = (char*)malloc(phdr->n_data_bytes);
	if (NULL == inbuf) {
		return 100; // buffer allocation failed
	}
	
	pfin->read(inbuf, phdr->n_data_bytes);
	if (pfin->fail()) {
		nerr = 110; // reading data from file failed
		goto _exit_point;
	}

	n = pfhdr->n_columns*pfhdr->n_rows;
	// switch depending on data type
	switch (pfhdr->n_bpi) {
	case 8:
		pdata8 = (unsigned __int8*)(inbuf);
		for (i = 0; i < n; i++) {
			buf[i] = (double)pdata8[i];
		}
		break;
	case 16:
		if (swapbytes) {
			for (i = 0; i < n; i++) {
				j = 2 * i;
				bs[1] = inbuf[j];
				bs[0] = inbuf[j + 1];
				buf[i] = (double)(*((unsigned __int16*)bs));
			}
		}
		else {
			pdata16 = (unsigned __int16*)(inbuf);
			for (i = 0; i < n; i++) {
				buf[i] = (double)pdata16[i];
			}
		}
		break;
	case 32:
		if (swapbytes) {
			for (i = 0; i < n; i++) {
				j = 4 * i;
				bs[3] = inbuf[j];
				bs[2] = inbuf[j + 1];
				bs[1] = inbuf[j + 2];
				bs[0] = inbuf[j + 3];
				buf[i] = (double)(*((unsigned __int32*)bs));
			}
		}
		else {
			pdata32 = (unsigned __int32*)(inbuf);
			for (i = 0; i < n; i++) {
				buf[i] = (double)pdata32[i];
			}
		}
		
		break;
	default:
		return 200; // unsupported data type
		break;
	}

_exit_point:
	if (NULL != inbuf) {
		free(inbuf);
	}
	return nerr;
}



