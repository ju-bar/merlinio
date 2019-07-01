// file : "merlin_prm.cpp"
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// Defines parameter structures used by merlinio
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


#include "pch.h"
#include "merlin_prm.h"

merlin_params::merlin_params()
{
	btalk = true;
	binteractive = false;
	bscanframeheaders = false;
	swapbytes = false;
	gaincorrect = false;
	defects_modified = false;

	ndebug = 0;

	frame_calib.offset = { 0.,0. };
	frame_calib.a0 = { 1., 0. };
	frame_calib.a1 = { 0., 1. };

	range_annular = { 0., 0. };

	scan_rect_roi = { 0, 0, 0, 0 };

	str_file_input = "input";
	str_file_output = "output";
	str_file_ctrl = "merlinio_control";

	v_defect_corr.clear();
	img_gaincorrect = NULL;
}

merlin_params::~merlin_params()
{
	v_frm_file.clear();
	v_frm_pos.clear();
	v_str_ctrl.clear();
	v_defect_corr.clear();
	if (NULL != img_gaincorrect) { free(img_gaincorrect); }
}

int merlin_params::read_header(void)
{
	int nerr = 0;
	bool bfilefound = true;
	int nfile = 0;
	std::string str_header = str_file_input + ".hdr";
	std::string str_datfile = "";
	std::string str_line = "";
	std::ifstream fin;
	if (btalk) {
		std::cout << std::endl;
		std::cout << "Reading information from merlin header file: " << str_header << std::endl;
	}
	// open file and read all header lines
	fin.open(str_header);
	nerr = merlin_read_header(&fin, &hdr);
	fin.close();

	if (btalk) { // tell infos
		std::cout << "- timestamp: " << hdr.s_timestamp << std::endl;
		std::cout << "- # frames: " << hdr.n_frames << std::endl;
		std::cout << "- # columns: " << hdr.n_columns << std::endl;
		std::cout << "- # rows: " << hdr.n_rows << std::endl;
	}
	return 0;
}

int merlin_params::read_frame_headers(void)
{
	int nerr = 0; // error code
	bool bfilefound = true;
	bool bread = true; // header reading flag
	bool bconsistent = true;
	int ierr = 0;
	int i = 0, lfile = 0, nextfrm = 0;
	int n_frm = 0; // count number of input frames
	size_t i_pos = 0; // read position
	struct stat statbuf;
	std::string str_file = "";
	std::string str_tmp = "";
	std::streampos fpos = 0, dpos = 0, lfpos = 0;
	std::ifstream fin;
	merlin_frame_hdr fhdr;
	size_t nhsize = 0, ndsize = 0;

	hdr.n_files = 0; // reset number of files
	v_frm_file.clear(); // clear list of frame to file indices
	v_frm_pos.clear(); // clear list of frame to file positions
	hdr.n_fhdr_bytes = 0; // reset header size
	hdr.n_data_bytes = 0; // reset data size
	// get information from the first frame header -> frame size and number of bytes per item
	while (bfilefound) { // loop finding all files
		str_file = str_file_input + std::to_string(hdr.n_files + 1) + ".mib"; // file name construction
		bfilefound = (stat(str_file.c_str(), &statbuf) == 0); // check if file exists
		if (bfilefound) { // file exists
			fin.open(str_file); // open the file
			if (fin.is_open()) { // ...
				bread = true; // always read the first header of a file
				while (bread) { // loop in case of header scanning request
					
					// read the next frame header
					ierr = merlin_read_frame_header(&fin, &fhdr);
					if (0 == ierr) { // success
						if (n_frm == 0) { // this is the first frame of the whole data set
							// transfer frame data to object member ... This is the template for all headers
							hdr_frm = fhdr;
							// set size infos from frame header information to global header hdr
							hdr.n_fhdr_bytes = (size_t)hdr_frm.n_size;
							hdr.n_data_bytes = ((size_t)hdr_frm.n_columns*hdr_frm.n_rows*hdr_frm.n_bpi >> 3);
							dpos = (std::streampos)(hdr.n_fhdr_bytes + hdr.n_data_bytes); // remember the expected frame shift in the file
						}
						fpos = fin.tellg(); // position in file after the header
						if (bscanframeheaders) { // scan all headers individually ...
							// check consistency of the current header to the header template
							// ... only the numbers which matter for the header and data size
							bconsistent = (
								hdr_frm.n_size == fhdr.n_size &&
								hdr_frm.n_bpi == fhdr.n_bpi &&
								hdr_frm.n_columns == fhdr.n_columns &&
								hdr_frm.n_rows == fhdr.n_rows
								);
							bconsistent &= (fhdr.i_seq == n_frm); // include the sequence consistency
							if (bconsistent) { // we conclude that this frame belongs to the sequence
								v_frm_pos.push_back(fpos); // store data offset for this frame
								v_frm_file.push_back(hdr.n_files); // store file index for this frame
								n_frm++; // increment the frame counter
								fpos += (std::streampos)hdr.n_data_bytes; // next header position
								fin.seekg(fpos); // try stepping to the next data position in the file
								bread = fin.good(); // keep on reading as long as the file is good
							}
							else {
								bread = false;
								std::cerr << "Error: inconsistent frame header data.\n";
								nerr = 3;
								goto _cancel_file_browsing;
							}
						}
						else { // no full scan - assuming regular (same header) file content
							if (n_frm == 0) { // this is the first item
								// add fhdr right away
								v_frm_pos.push_back(fpos); // store data offset for this frame
								v_frm_file.push_back(hdr.n_files); // store file index for this frame
								n_frm++; // increment the frame counter
							}
							else { // this is a second file ...
								lfpos = v_frm_pos[n_frm-1]; // get last added frame position in previous file
								lfile = v_frm_file[n_frm - 1]; // get last added file index (previous file)
								nextfrm = n_frm; // this is the next frame index expected in the sequence
								if (fhdr.i_seq > n_frm) { // frames are missing in the sequence between the last added and the new fhdr
									for (i = nextfrm; i < fhdr.i_seq; i++) {
										// add frames
										lfpos += dpos; // shift to the next frame in the file
										v_frm_pos.push_back(lfpos); // store the position
										v_frm_file.push_back(lfile); // store the file index
										n_frm++;
									}
								}
								// add fhdr
								v_frm_pos.push_back(fpos); // store data offset for this frame
								v_frm_file.push_back(hdr.n_files); // store file index for this frame
								n_frm++; // increment the frame counter
							} // ... if (n_frm == 0) ...
							bread = false; // stop reading headers from this file
						} // ... if (bscanframeheaders) ...
					}
					else { // reading failed
						if (fin.eof()) {
							// end of file reached
							bread = false;
						}
						else {
							// another problem occurred
							std::cerr << "Error: failed to read a frame header (code " << ierr << ")." << std::endl;
							std::cerr << "Error: unknown format of file " << str_file << " " << std::endl;
							nerr = 2;
							goto _cancel_file_browsing;
						}
					}
				} // while (bread)
			}
			else { // file opening failed
				std::cerr << "Error: failed to open data file: " << str_file << std::endl;
				return 1; // (no need to close the file)
			}
		_cancel_file_browsing:
			fin.close(); // close the file
			if (nerr != 0) { // return in case of errors
				return nerr;
			}
			hdr.n_files++; // count number of files
		}
	} // while (bfilefound)

	if (hdr.n_files == 0 || n_frm == 0) {
		std::cerr << "Error: found no frame headers.\n";
		return 4;
	}

	if (hdr.n_files > 0 && n_frm < hdr.n_frames) { // expected frames are missing
		if (bscanframeheaders) { // ... errro in case of scan mode
			std::cerr << "Error: frame header scan is missing frames.\n";
			return 5;
		}
		else { // add missing frames in case of fast mode
			lfpos = v_frm_pos[n_frm - 1]; // get last added frame position in last file
			lfile = v_frm_file[n_frm - 1]; // get last added file index (previous file)
			nextfrm = n_frm; // this is the next frame index expected in the sequence
			for (i = nextfrm; i < hdr.n_frames; i++) {
				// add frames
				lfpos += dpos; // shift to the next frame in the file
				v_frm_pos.push_back(lfpos); // store the position
				v_frm_file.push_back(lfile); // store the file index
				n_frm++;
			}
		}
	}

	if (nerr == 0 && btalk && ndebug > 0) {
		std::cout << "- # files: " << hdr.n_files << std::endl;
		std::cout << "- # frame header bytes: " << hdr.n_fhdr_bytes << std::endl;
		std::cout << "- # frame data bytes: " << hdr.n_data_bytes << std::endl;
	}

	return nerr;
}


int merlin_params::read_param(int ipos, std::string * pstr, std::string * prm)
{
	std::size_t lipos = (std::size_t)(ipos >= 0 ? ipos : 0);
	std::size_t lcpos = lipos;
	std::size_t lmpos = 0;
	std::size_t lsep = 0;
	std::string str_sep = ", ";
	if (NULL == pstr) {
		return -1;
	}
	if (NULL == prm) {
		return -2;
	}
	lmpos = pstr->size();
	if (lcpos >= lmpos) {
		*prm = "";
		return (int)lcpos;
	}
	// find next separator ...
	lsep = pstr->find_first_of(str_sep, lcpos);
	if (lsep == pstr->npos) {
		lcpos = lmpos;
	}
	else {
		lcpos = lsep;
	}
	// lcpos is now the position of the next separator or eos
	if (lcpos - lipos > 0) { // extract sub-string
		*prm = pstr->substr(lipos, lcpos - lipos);
	}
	// find next non-separator character or eos
	if (lcpos < lmpos) {
		lsep = pstr->find_first_not_of(str_sep, lcpos);
		if (lsep == pstr->npos) { // eos
			lcpos = lmpos;
		}
		else { // next param
			lcpos = lsep;
		}
	}
	return (int)lcpos; // return position which the possible begin of a new parameter or eos
}



bool merlin_params::in_scan_roi(merlin_pix pos, merlin_roi roi)
{
	bool b_result = false;
	if (ndebug > 4) {
		std::cout << "merlin_params::in_scan_roi: pos=(" << pos.x << "," << pos.y << ")\n";
		std::cout << "merlin_params::in_scan_roi: roi=((" << roi.x0 << "," << roi.y0 << "),(" << roi.x1 << "," << roi.y1 << "))\n";
	}
	b_result = (pos.x >= roi.x0 && pos.x <= roi.x1 && pos.y >= roi.y0 && pos.y <= roi.y1);
	if (ndebug > 4) {
		if (b_result) {
			std::cout << "merlin_params::in_scan_roi = true\n";
		}
		else {
			std::cout << "merlin_params::in_scan_roi = false\n";
		}
	}
	return b_result;
}



int merlin_params::get_calib_pos(merlin_pix xin, merlin_pos * xout)
{
	double dx, dy;
	if (NULL == xout) {
		return 2; // missing output pointer
	}
	dx = (double)xin.x - frame_calib.offset.x;
	dy = (double)xin.y - frame_calib.offset.y;
	xout->x = dx * frame_calib.a0.x + dy * frame_calib.a1.x;
	xout->y = dx * frame_calib.a0.y + dy * frame_calib.a1.y;
	return 0;
}

int merlin_params::get_scan_pixel(int idx, int &x, int &y)
{
	if (hdr.n_columns <= 0) {
		return 1; // unreasonable # columns
	}
	x = imod(idx, hdr.n_columns);
	y = imod((idx - x) / hdr.n_columns, hdr.n_rows);
	return 0;
}

int merlin_params::get_frame_idx(int x, int y)
{
	if (hdr.n_columns <= 0) {
		return -1; // unreasonable # columns
	}
	if (x < 0 || x >= hdr.n_columns) {
		return -2; // invalid scan x position
	}
	if (y < 0 || y >= hdr.n_rows) {
		return -3; // invalid scan y position
	}
	return imod(x, hdr.n_columns) + imod(y, hdr.n_rows) * hdr.n_columns;
}

int merlin_params::get_frame_pixel(int idx, int &x, int &y)
{
	if (hdr_frm.n_columns <= 0) {
		return 1; // missing reasonable # columns
	}
	x = imod(idx, hdr_frm.n_columns);
	y = imod((idx - x) / hdr_frm.n_columns, hdr_frm.n_rows);
	return 0;
}

int merlin_params::get_frame_pixel_idx(int x, int y)
{
	if (hdr_frm.n_columns <= 0 || hdr_frm.n_rows <= 0) {
		return -1; // unreasonable # columns or rows
	}
	return imod(x, hdr_frm.n_columns) + imod(y, hdr_frm.n_rows) * hdr_frm.n_columns;
}

int merlin_params::get_frame_filepos(int idx, int &ifile, std::streampos &ipos)
{
	if (idx < 0 || idx >= hdr.n_frames) {
		return 1;
	}
	ifile = v_frm_file[idx];
	ipos = v_frm_pos[idx];
	return 0;
}


size_t merlin_params::get_scan_rect_roi_size(void)
{
	size_t dx = (scan_rect_roi.x1 > scan_rect_roi.x0 ? scan_rect_roi.x1 - scan_rect_roi.x0 : 0);
	size_t dy = (scan_rect_roi.y1 > scan_rect_roi.y0 ? scan_rect_roi.y1 - scan_rect_roi.y0 : 0);
	return dx * dy;
}


int merlin_params::set_scan_rect_roi(std::string str_roi)
{
	int i_pos = 0;
	int n_prm = 0;
	std::string str_num = "";
	if (ndebug > 3) {
		std::cout << "merlin_params::set_scan_rect_roi: str_roi=" << str_roi << std::endl;
	}
	while (i_pos >= 0 && n_prm < 4) {
		i_pos = read_param(i_pos, &str_roi, &str_num);
		if (i_pos < 0) { return 1 + n_prm; } // parsing error
		switch (n_prm) {
		case 0:
			scan_rect_roi.x0 = atoi(str_num.c_str());
			if (scan_rect_roi.x0 < 0) {
				std::cerr << "Warning: scan roi x0 out of bounds: < 0.\n";
			}
			break;
		case 1:
			scan_rect_roi.y0 = atoi(str_num.c_str());
			if (scan_rect_roi.y0 < 0) {
				std::cerr << "Warning: scan roi y0 out of bounds: < 0.\n";
			}
			break;
		case 2:
			scan_rect_roi.x1 = atoi(str_num.c_str());
			if (scan_rect_roi.x1 > hdr.n_columns) {
				std::cerr << "Warning: scan roi x1 out of bounds: > " << hdr.n_columns - 1 << ".\n";
			}
			break;
		case 3:
			scan_rect_roi.y1 = atoi(str_num.c_str());
			if (scan_rect_roi.y1 > hdr.n_rows) {
				std::cerr << "Warning: scan roi y1 out of bounds: > " << hdr.n_rows - 1 << ".\n";
			}
			break;
		}
		n_prm++;
	}
	if (ndebug > 3) {
		std::cout << "merlin_params::in_scan_roi: merlin_params::scan_rect_roi=((" << scan_rect_roi.x0 << "," << scan_rect_roi.y0 << "),(" << scan_rect_roi.x1 << "," << scan_rect_roi.y1 << "))\n";
	}
	return 0;
}

int merlin_params::set_origin(std::string str_org)
{
	int i_pos = 0;
	int n_prm = 0;
	std::string str_num = "";
	if (ndebug > 3) {
		std::cout << "merlin_params::set_origin: str_org=" << str_org << std::endl;
	}
	while (i_pos >= 0 && n_prm < 2) {
		i_pos = read_param(i_pos, &str_org, &str_num);
		if (i_pos < 0) { return 1 + n_prm; } // parsing error
		switch (n_prm) {
		case 0:
			frame_calib.offset.x = (double)atof(str_num.c_str());
			break;
		case 1:
			frame_calib.offset.y = (double)atof(str_num.c_str());
			break;
		}
		n_prm++;
	}
	if (ndebug > 3) {
		std::cout << "merlin_params::set_origin: merlin_params::frame_calib.offset=(" << frame_calib.offset.x << "," << frame_calib.offset.y << ")\n";
	}
	return 0;
}

int merlin_params::set_sampling(std::string str_samp)
{
	int i_pos = 0;
	int n_prm = 0;
	std::string str_num = "";
	if (ndebug > 3) {
		std::cout << "merlin_params::set_sampling: str_samp=" << str_samp << std::endl;
	}
	while (i_pos >= 0 && n_prm < 4) {
		i_pos = read_param(i_pos, &str_samp, &str_num);
		if (i_pos < 0) { return 1 + n_prm; } // parsing error
		switch (n_prm) {
		case 0:
			frame_calib.a0.x = (double)atof(str_num.c_str());
			break;
		case 1:
			frame_calib.a1.x = (double)atof(str_num.c_str());
			break;
		case 2:
			frame_calib.a0.y = (double)atof(str_num.c_str());
			break;
		case 3:
			frame_calib.a1.y = (double)atof(str_num.c_str());
			break;
		}
		n_prm++;
	}
	if (ndebug > 3) {
		std::cout << "merlin_params::set_sampling: merlin_params::frame_calib.a0=(" << frame_calib.a0.x << "," << frame_calib.a0.y << ")\n";
		std::cout << "merlin_params::set_sampling: merlin_params::frame_calib.a1=(" << frame_calib.a1.x << "," << frame_calib.a1.y << ")\n";
	}
	return 0;
}

int merlin_params::set_annular_range(std::string str_rng)
{
	int i_pos = 0;
	int n_prm = 0;
	std::string str_num = "";
	if (ndebug > 3) {
		std::cout << "set_annular_range::set_annular_range: str_samp=" << str_rng << std::endl;
	}
	while (i_pos >= 0 && n_prm < 2) {
		i_pos = read_param(i_pos, &str_rng, &str_num);
		if (i_pos < 0) { return 1 + n_prm; } // parsing error
		switch (n_prm) {
		case 0:
			range_annular.min = (double)atof(str_num.c_str());
			break;
		case 1:
			range_annular.max = (double)atof(str_num.c_str());
			break;
		}
		n_prm++;
	}
	if (ndebug > 3) {
		std::cout << "merlin_params::set_annular_range: merlin_params::range_annular=(" << range_annular.min << "," << range_annular.max << ")\n";
	}
	return 0;
}

bool merlin_params::is_defect_pixel(size_t idx)
{
	size_t ndef = v_defect_corr.size();
	if (ndef > 0) {
		size_t didx = 0;
		for (didx = 0; didx < ndef; didx++) {
			if (v_defect_corr[didx].idx == idx) return true; // found
		}
	}
	return false;
}

bool merlin_params::is_defect_pixel(int x, int y)
{
	size_t i = get_frame_pixel_idx(x, y);
	if (i>=0) return is_defect_pixel(i);
	return false;
}

bool merlin_params::is_defect_list_modified(void)
{
	return defects_modified;
}

int merlin_params::update_defect_correction_list(void)
{
	size_t ncorr = 0, nc = 0;
	size_t ndef = v_defect_corr.size();
	size_t i = 0, j = 0, idx = 0, idx2 = 0;
	int k = 0, l = 0, x = 0, y = 0;
	if (ndef > 0) { // update correction lists
		for (i = 0; i < ndef; i++) { // .. for all defects
			v_defect_corr[i].v_idx_corr.clear(); // delete correction table
			idx = v_defect_corr[i].idx;
			x = v_defect_corr[i].x;
			y = v_defect_corr[i].y;
			// generate new correction table
			for (k = -1; k <= 1; k++) {
				for (l = -1; l <= 1; l++) {
					idx2 = (size_t)get_frame_pixel_idx(x + l, y + k); // get the frame stream index for (i1+l,j1+k)
					if (idx2 < 0) continue; // invald index ?
					if (is_defect_pixel(idx2)) continue; // is also a defect ?
					v_defect_corr[i].v_idx_corr.push_back(idx2); // this pixel can be used for correction
					nc++;
				}
				if (ndebug > 3) {
					std::cout << "- registered " << nc << " correction pixels for defect " << i+1 << " at (" << x << "," << y << "), index = " << idx << " \n";
				}
			}
		}
	}
	defects_modified = false; // mark all defects to be fully registered with correction table
	return 0;
}

int merlin_params::set_defect_pixel(int x, int y)
{
	if (is_defect_pixel(x, y)) return 0; // this is already set as defect pixel
	defect_pixel_corr dpc;
	size_t i = get_frame_pixel_idx(x, y);
	if (i < 0) return 1; // error
	dpc.idx = i;
	get_frame_pixel((int)i, dpc.x, dpc.y); // get 2d pixel indices in the grid -> (i1,j1)
	// This is a new defect. -> Add with an empty correction list.
	dpc.v_idx_corr.clear();
	// Push the new defect to the list.
	v_defect_corr.push_back(dpc);
	// Mark the defect list as modified.
	defects_modified = true;
	return 0;
}


int merlin_params::set_defect_pixel(std::string str_pos)
{
	int x = 0, y = 0;
	int i_pos = 0;
	int n_prm = 0;
	std::string str_num = "";
	if (ndebug > 3) {
		std::cout << "merlin_params::set_defect_pixel: str_pos=" << str_pos << std::endl;
	}
	while (i_pos >= 0 && n_prm < 2) {
		i_pos = read_param(i_pos, &str_pos, &str_num);
		if (i_pos < 0) { return 1 + n_prm; } // parsing error
		switch (n_prm) {
		case 0:
			x = atoi(str_num.c_str());
			break;
		case 1:
			y = atoi(str_num.c_str());
			break;
		}
		n_prm++;
	}
	if (ndebug > 3) {
		std::cout << "merlin_params::set_defect_pixel: pos=(" << x << "," << y << ")\n";
	}
	return set_defect_pixel(x, y);
}

int merlin_params::load_defect_mask(std::string str_file)
{
	int nx = hdr_frm.n_columns;
	int ny = hdr_frm.n_rows;
	size_t npix = (size_t)ny*((size_t)nx);
	size_t nbytes = sizeof(int)*npix;
	int * img_defectmask = NULL;
	int num_defects = 0;
	size_t idx = 0, idx2 = 0;
	int i = 0, j = 0, l = 0, k = 0;
	bool bfilefound = false;
	struct stat statbuf;
	std::ifstream fin;
	defect_pixel_corr dpc;
	if (npix == 0) {
		std::cerr << "merlin_params::load_defect_mask: failed due to invalid frame size.\n";
		return 1;
	}
	bfilefound = (stat(str_file.c_str(), &statbuf) == 0); // check if file exists
	if (!bfilefound) {
		std::cerr << "merlin_params::load_defect_mask: file [" << str_file << "] not found.\n";
		return 2;
	}
	fin.open(str_file); // open the file
	if (fin.is_open()) { // ...
		img_defectmask = (int*)calloc(npix, sizeof(int));
		fin.read((char*)img_defectmask, nbytes);
		if (fin.fail() && (!fin.eof())) {
			std::cerr << "merlin_params::load_defect_mask: failed to read data.\n";
			fin.close();
			free(img_defectmask);
			return 4;
		}
		fin.close();
		//
		// add defects to the list
		//
		for (idx = 0; idx < npix; idx++) {
			if (img_defectmask[idx] != 0) { // defect pixel
				if (is_defect_pixel(idx)) continue; // this defect is already registered
				dpc.idx = idx;
				get_frame_pixel((int)idx, dpc.x, dpc.y);
				dpc.v_idx_corr.clear();
				v_defect_corr.push_back(dpc); // add the defect data to the list
				num_defects++;
			}
		}
		if (num_defects > 0) { // there are defect pixels
			if (btalk) {
				std::cout << "- " << num_defects << " pixels were added to the defect list.\n";
			}
			defects_modified = true; // mark the defect list as modiefied
		}
		else { // num_defects == 0 -> no defects
			if (btalk) {
				std::cout << "- the loaded defect mask has no effect.\n";
			}
		}
		if (img_defectmask != NULL) { free(img_defectmask); }
	}
	else {
		std::cerr << "merlin_params::load_defect_mask: failed to open file [" << str_file << "].\n";
		return 3;
	}
	return 0;
}

int merlin_params::load_defect_list(std::string str_file)
{
	int num_defects = (int)v_defect_corr.size();
	bool bfilefound = false;
	size_t nlines = 0;
	struct stat statbuf;
	std::ifstream fin;
	std::string str_line;
	defect_pixel_corr dpc;
	bfilefound = (stat(str_file.c_str(), &statbuf) == 0); // check if file exists
	if (!bfilefound) {
		std::cerr << "merlin_params::load_defect_list: file [" << str_file << "] not found.\n";
		return 2;
	}
	fin.open(str_file); // open the file
	if (fin.is_open()) { // ...
		// read lines and parse for coordinates
		while (fin.good()) {
			str_line.clear();
			getline(fin, str_line);
			nlines++;
			if (str_line.size() > 0) {
				// add defect to the list
				if (0 != set_defect_pixel(str_line)) {
					std::cerr << "Error: failed to set defect pixel from line " << nlines << " of the list file.\n";
					std::cout << "     : " << str_line << std::endl;
				}
			}
		}
		fin.close();
		//
		num_defects = (int)v_defect_corr.size() - num_defects;
		if (num_defects > 0) { // there are defect pixels added
			if (btalk) {
				std::cout << "- " << num_defects << " pixels were added to the defect list.\n";
			}
			defects_modified = true; // mark the defect list as modiefied
		}
		else { // num_defects == 0 -> no defects
			if (btalk) {
				std::cout << "- the loaded defect list has no effect.\n";
			}
		}
	}
	else {
		std::cerr << "merlin_params::load_defect_list: failed to open file [" << str_file << "].\n";
		return 3;
	}
	return 0;
}


int merlin_params::unset_defect_list(void)
{
	v_defect_corr.clear();
	defects_modified = false;
	return 0;
}

int merlin_params::unset_defect_pixel(int x, int y)
{
	size_t idx_d = get_frame_pixel_idx(x, y);
	size_t ndef = v_defect_corr.size();
	size_t i = 0;
	if (idx_d >= 0 && ndef > 0) {
		for (i = 0; i < ndef; i++) {
			if (v_defect_corr[i].idx == idx_d) { // remove this
				v_defect_corr.erase(v_defect_corr.begin()+i);
				if (v_defect_corr.size() > 0) {
					defects_modified = true;
				}
				else {
					v_defect_corr.clear();
					defects_modified = false;
				}
				break; // stop the loop
			}
		}
	}
	return 0;
}

int merlin_params::unset_defect_pixel(std::string str_pos)
{
	int x = 0, y = 0;
	int i_pos = 0;
	int n_prm = 0;
	std::string str_num = "";
	if (ndebug > 3) {
		std::cout << "merlin_params::unset_defect_pixel: str_pos=" << str_pos << std::endl;
	}
	while (i_pos >= 0 && n_prm < 2) {
		i_pos = read_param(i_pos, &str_pos, &str_num);
		if (i_pos < 0) { return 1 + n_prm; } // parsing error
		switch (n_prm) {
		case 0:
			x = atoi(str_num.c_str());
			break;
		case 1:
			y = atoi(str_num.c_str());
			break;
		}
		n_prm++;
	}
	if (ndebug > 3) {
		std::cout << "merlin_params::unset_defect_pixel: pos=(" << x << "," << y << ")\n";
	}
	return unset_defect_pixel(x, y);
}


int merlin_params::load_gain_correction(std::string str_file)
{
	size_t npix = (size_t)hdr_frm.n_rows*((size_t)hdr_frm.n_columns);
	size_t nbytes = sizeof(float)*npix;
	float * inbuf = NULL;
	bool bfilefound = false;
	struct stat statbuf;
	std::ifstream fin;
	if (npix == 0) {
		std::cerr << "merlin_params::load_gain_correction: failed due to invalid frame size.\n";
		return 1;
	}
	bfilefound = (stat(str_file.c_str(), &statbuf) == 0); // check if file exists
	if (!bfilefound) {
		std::cerr << "merlin_params::load_gain_correction: file [" << str_file << "] not found.\n";
		return 2;
	}
	if (ndebug > 0) {
		std::cout << "opening file " << str_file << std::endl;
	}
	fin.open(str_file); // open the file
	if (fin.is_open()) { // ...
		if (gaincorrect) unset_gain_correction();
		inbuf = (float*)calloc(npix, sizeof(float)); // read 32-bit float data
		if (ndebug > 0) {
			std::cout << "- loading " << nbytes << " bytes " << std::endl;
		}
		fin.read((char*)inbuf, nbytes);
		if (fin.fail() && (!fin.eof())) { // failed ?
			std::cerr << "merlin_params::load_gain_correction: failed to read data.\n";
			fin.close();
			gaincorrect = false;
			free(inbuf);
			return 4;
		}
		fin.close();
		if (ndebug > 3) {
			std::cout << "- transferring data to internal memory" << std::endl;
		}
		img_gaincorrect = (double*)calloc(npix, sizeof(double)); // type-cast factors to double for internal use
		for (size_t i = 0; i < npix; i++) {
			img_gaincorrect[i] = (double)inbuf[i];
		}
		free(inbuf);
		gaincorrect = true; // done
		if (btalk) {
			std::cout << "- gain correction factors loaded successfully.\n";
		}
	}
	else {
		std::cerr << "merlin_params::load_gain_correction: failed to open file [" << str_file << "].\n";
		return 3;
	}
	return 0;
}

int merlin_params::unset_gain_correction(void)
{
	gaincorrect = false;
	if (NULL != img_gaincorrect) { free(img_gaincorrect); }
	return 0;
}

int merlin_params::gain_correction(double * buf)
{
	size_t npix = (size_t)hdr_frm.n_rows * (size_t)hdr_frm.n_columns;
	size_t i = 0;
	if (gaincorrect) {
		if (buf == NULL || img_gaincorrect == NULL) return 1;
		if (npix == 0) return 2;
		for (i = 0; i < npix; i++) {
			buf[i] = buf[i] * img_gaincorrect[i];
		}
	}
	return 0;
}


int merlin_params::defect_correction(double * buf)
{
	size_t npix = (size_t)hdr_frm.n_rows * (size_t)hdr_frm.n_columns;
	size_t idx = 0, idx2 = 0;
	int idef = 0, nc = 0;
	int icor = 0;
	double rnc = 0.;
	double snc = 0.;
	if (v_defect_corr.size() > 0) {
		if (buf == NULL) return 1;
		if (npix == 0) return 2;
		if (ndebug > 3) {
			std::cout << "- correcting " << v_defect_corr.size() << " defect pixels\n";
		}
		for (idef = 0; idef < v_defect_corr.size(); idef++) {
			idx = v_defect_corr[idef].idx; // get defect pixel index
			nc = (int)v_defect_corr[idef].v_idx_corr.size(); // get number of correction pixels
			if (nc > 0) {
				if (ndebug > 3) {
					std::cout << "- defect " << idef << " (" << idx << ") from " << nc << " neighbor pixels\n";
				}
				rnc = 1. / ((double)nc);
				snc = 0.;
				for (icor = 0; icor < nc; icor++) { // loop over correction pixels
					idx2 = v_defect_corr[idef].v_idx_corr[icor]; // get index of correction pixel
					snc += buf[idx2]; // accumulate intensity
				}
				buf[idx] = snc * rnc; // set correction to defect pixel
			}
		}
	}
	return 0;
}

