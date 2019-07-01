// file: 'merlinio.cpp'
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// This file contains the 'main' function.
// Program execution begins and ends there.
//
// This program reads data from merlin detector output files.
// The data is dumped again to disk as binary raw data.
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
#include <algorithm>
#include <cctype>


merlin_params prm;


// -----------------------------------------------------------------------------
//
// Program parameter I/O
//
// -----------------------------------------------------------------------------

int parseoptions(int argc, char* argv[])
{
	prm.btalk = true;
	prm.bscanframeheaders = false;
	prm.ndebug = 0;
	if (argc > 1) {
		prm.str_file_input = argv[1]; // expecting input file name string as first argument (no number and extension)

		if (argc > 2) {
			std::string cmd;
			// parse other arguments
			for (int iarg = 2; iarg < argc; iarg++) {
				cmd = argv[iarg];
				if (cmd == "/silent") {
					prm.btalk = false;
					continue;
				}
				if (cmd == "/debug") { // direct debug switch
					prm.ndebug = 1; // set debug level 1
					continue;
				}
				if (cmd == "/sfh" || cmd == "/scanframeheaders") {
					prm.bscanframeheaders = true; // force careful frame header scan
					continue;
				}

				if (cmd == "-dbgl") { // leveled debug switch + level
					iarg++;
					if (iarg >= argc) {
						std::cerr << "Error: expecting a debug level number after option -dbgl.\n";
						return 1;
					}
					prm.ndebug = atoi(argv[iarg]);
					if (prm.ndebug < 0) prm.ndebug = 0;
					if (prm.ndebug > 5) prm.ndebug = 5;
					continue;
				}

				if (cmd == "-o" || cmd =="-output") { // modified output name + name
					iarg++;
					if (iarg >= argc) {
						std::cerr << "Error: expecting a file name string after option -o.\n";
						return 1;
					}
					prm.str_file_output = argv[iarg];
					continue;
				}
				if (cmd == "-c" || cmd == "-control") { // modified name of control file
					iarg++;
					if (iarg >= argc) {
						std::cerr << "Error: expecting a file name string after option -control (-c).\n";
						return 1;
					}
					prm.str_file_ctrl = argv[iarg];
					continue;
				}
			}
		}
	}
	if (prm.ndebug > 0) { // debug overrides
		prm.btalk = true;
	}
	return 0;
}


// -----------------------------------------------------------------------------
//
// DATA I/O
//
// -----------------------------------------------------------------------------


int write_data(char* buf, size_t nbytes, std::string str_file)
{
	int nerr = 0;
	std::ofstream fout;
	fout.open(str_file, std::ios::trunc | std::ios::binary);
	if (fout.is_open()) {
		fout.write(buf, nbytes);
		if (fout.fail()) {
			std::cerr << "Error: failed to write data to file " << str_file << ".\n";
			nerr = 2;
		}
		fout.close();
	}
	else {
		std::cerr << "Error: failed to open output file " << str_file << ".\n";
		nerr = 1;
	}
	return nerr;
}


// -----------------------------------------------------------------------------
//
// DATA processing functions
//
// -----------------------------------------------------------------------------



int prepare_annular_detector(size_t nlen, double * detbuf, int * dethash, size_t * nhash, merlin_frame_hdr * pfhdr)
{
	size_t i = 0, j = 0;
	double qm = 0.;
	merlin_pix p;
	merlin_pos q;
	bool bhash = false;
	if (nlen == 0) {
		return 1; // invalid parameter 1
	}
	if (NULL == detbuf) {
		return 2; // missing parameter 2
	}
	if (NULL == pfhdr) {
		return 5; // missing parameter 5
	}

	if (dethash != NULL && nhash != NULL) {
		*nhash = 0;
		bhash = true;
	}
	for (i = 0; i < nlen; i++) { // loop through array
		if (0 != prm.get_frame_pixel((int)i, p.x, p.y)) {
			return 10; // failed to get scan position
		}
		if (0 != prm.get_calib_pos(p, &q)) {
			return 20; // failed to get calibrated position
		}
		// hard mask // replace by detector transfer function later
		qm = sqrt(q.x*q.x + q.y*q.y);
		if (qm >= prm.range_annular.min && qm < prm.range_annular.max) {
			detbuf[i] = 1.;
			if (bhash) {
				j = *nhash;
				dethash[j] = (int)i;
				*nhash = j + 1;
			}
		}
		else {
			detbuf[i] = 0.;
		}
	}
	return 0;
}

int prepare_frame_coordinates(size_t nlen, double * x, double * y, merlin_frame_hdr * pfhdr)
{
	size_t i = 0, j = 0;
	double qm = 0.;
	merlin_pix p;
	merlin_pos q;
	if (nlen == 0) {
		return 1; // invalid parameter 1
	}
	if (NULL == x) {
		return 2; // missing parameter 2
	}
	if (NULL == y) {
		return 3; // missing parameter 3
	}
	if (NULL == pfhdr) {
		return 4; // missing parameter 4
	}
	for (i = 0; i < nlen; i++) { // loop through array
		if (0 != prm.get_frame_pixel((int)i, p.x, p.y)) {
			return 10; // failed to get scan position
		}
		if (0 != prm.get_calib_pos(p, &q)) {
			return 20; // failed to get calibrated position
		}
		// store positions
		x[i] = q.x;
		y[i] = q.y;
	}
	return 0;
}

int sum_annular_range(size_t nlen, double * buf, double * detbuf, int * dethash, size_t nhash, double * res)
{
	double lres = 0.0;
	size_t i = 0, j=0;
	if (NULL == buf) {
		return 2; // invalid input pointer, parameter 2
	}
	if (NULL == detbuf) {
		return 3; // invalid input pointer, parameter 3
	}
	if (NULL == res) {
		return 6; // invalid input pointer, parameter 6
	}
	*res = 0.;
	if (nlen > 0) {
		if (dethash != NULL && nhash > 0) { // hash-based summation
			for (i = 0; i < nhash; i++) {
				j = (size_t)dethash[i];
				lres += buf[j] * detbuf[j];
			}
		}
		else { // standard summation
			for (i = 0; i < nlen; i++) {
				lres += buf[i] * detbuf[i];
			}
		}
		*res = lres;
	}
	return 0;
}

int com_annular_range(size_t nlen, double * buf, double * detbuf, double * x, double * y, int * dethash, size_t nhash, double ref0, double * resx, double * resy)
{
	double lresx = 0.0;
	double lresy = 0.0;
	size_t i = 0, j = 0;
	if (NULL == buf) {
		return 2; // invalid input pointer, parameter 2
	}
	if (NULL == detbuf) {
		return 3; // invalid input pointer, parameter 3
	}
	if (NULL == x) {
		return 4; // invalid input pointer, parameter 4
	}
	if (NULL == y) {
		return 5; // invalid input pointer, parameter 5
	}
	if (NULL == resx) {
		return 9; // invalid input pointer, parameter 9
	}
	if (NULL == resy) {
		return 10; // invalid input pointer, parameter 10
	}
	*resx = 0.;
	*resy = 0.;
	if (nlen > 0 && ref0 > 0.) {
		if (dethash != NULL && nhash > 0) { // hash-based summation
			for (i = 0; i < nhash; i++) {
				j = (size_t)dethash[i];
				lresx += (x[j] * buf[j] * detbuf[j]);
				lresy += (y[j] * buf[j] * detbuf[j]);
			}
		}
		else { // standard summation
			for (i = 0; i < nlen; i++) {
				lresx += (x[i] * buf[i] * detbuf[i]);
				lresy += (y[i] * buf[i] * detbuf[i]);
			}
		}
		*resx = lresx / ref0;
		*resy = lresy / ref0;
	}
	return 0;
}

// -----------------------------------------------------------------------------
//
// SET control parameter functions
//
// -----------------------------------------------------------------------------




// -----------------------------------------------------------------------------
//
// RUN control set functions
//
// -----------------------------------------------------------------------------


int run_extract_frames()
{
	int nerr = 0;
	int i_frm = 0; // frame index
	size_t scan_pix = (size_t)prm.hdr.n_columns * prm.hdr.n_rows; // number of scan items
	char* datbuf = NULL;
	merlin_pix scan_pos;
	std::string str_file; // file name
	std::ifstream fin; // input stream
	std::ofstream fout; // output stream
	std::streampos fpos; // file position
	int ncfidx = -1;
	int fidx = -1; // file index
	int file_frm = -1; // frame index in file
	int prog_pct = 0;
	int prog_pct_old = 0;
	size_t nres = 0; // number of result items
	fout.open(prm.str_file_output, std::ios::binary | std::ios::trunc); // open output file for writing binary data
	if (!fout.is_open()) {
		std::cerr << "Error: failed to open output file " << prm.str_file_output << " for writing data.\n";
		return 1;
	}
	if (prm.hdr.n_data_bytes > 0 && prm.hdr.n_frames > 0) {
		datbuf = (char*)calloc(prm.hdr.n_data_bytes, sizeof(char));
		if (prm.btalk) {
			std::cout << "- extracting frames in current scan roi ...\n";
			std::cout << "  0 %\r";
		}
		for (i_frm = 0; i_frm < prm.hdr.n_frames; i_frm++) {
			nerr = prm.get_scan_pixel((int)i_frm, scan_pos.x, scan_pos.y);
			if (nerr != 0) { // integration failed
				std::cerr << "Error: failed to determine scan position for frame # " << i_frm << " (code " << nerr << ").\n";
				nerr = 100;
				goto _cancel_point; // stop working
			}
			if (prm.in_scan_roi(scan_pos, prm.scan_rect_roi)) {
				nerr = prm.get_frame_filepos(i_frm, fidx, fpos);
				if (nerr != 0) { // couldn't determine file index
					std::cerr << "Error: failed to determine file index and datat offset for frame # " << i_frm << std::endl;
					nerr = 101;
					goto _cancel_point; // stop working
				}
				if (fidx != ncfidx && fidx >= 0) { // not the right file is open
					if (ncfidx >= 0 && fin.is_open()) { // close the wrong file
						fin.close();
					}
					str_file = prm.str_file_input + std::to_string(fidx + 1) + ".mib"; // file name construction
					fin.open(str_file, std::ios::binary); // open the right file
					if (fin.is_open()) { // successfully opened
						ncfidx = fidx; // update current file index
					}
					else { // failure opening file
						std::cerr << "Error: failed to open input file: " << str_file << " for reading data.\n";
						nerr = 102;
						goto _cancel_point; // stop working
					}
				}
				fin.seekg(fpos);
				if (!fin.fail()) {
					fin.read(datbuf, prm.hdr.n_data_bytes); // read from input file
					if (fin.fail()) {
						std::cerr << "Error: failed reading data from input file: " << str_file << std::endl;
						nerr = 103;
						goto _cancel_point; // stop working
					}
					fout.write(datbuf, prm.hdr.n_data_bytes); // append to output file
					if (fout.fail()) {
						std::cerr << "Error: failed writing data to output file: " << prm.str_file_output << std::endl;
						nerr = 104;
						goto _cancel_point; // stop working
					}
					nres++; // increment result numbers
				}
				else {
					std::cerr << "Error: failed to position file pointer to frame # " << i_frm << ".\n";
					nerr = 106;
					goto _cancel_point; // stop working
				}
			} // if in roi
			// 
			prog_pct = (int)(100. * (double)i_frm / (double)prm.hdr.n_frames); // progress in percent
			if (prm.btalk && prog_pct > prog_pct_old) { // progress step ...
				std::cout << "  " << prog_pct << " %\r";
				prog_pct_old = prog_pct;
			}
		} // frame loop
	_cancel_point:
		if (datbuf) free(datbuf);
		if (fin.is_open()) {
			fin.close();
		}
		if (fout.is_open()) {
			fout.close();
			if (prm.btalk) {
				std::cout << "- written " << nres << " frames to file " << prm.str_file_output << std::endl;
				std::cout << "  bits per item: " << (int)prm.hdr_frm.n_bpi << std::endl;
				std::cout << "  items per frame: " << prm.hdr_frm.n_columns*prm.hdr_frm.n_rows << std::endl;
			}
		}
		if (nerr == 0) { // write the info file
			fout.open(prm.str_file_output + ".hdr", std::ios::trunc);
			if (fout.is_open()) {
				fout << "File name: " << prm.str_file_output << std::endl;
				fout << "Number of frames: " << nres << std::endl;
				fout << "Frame columns: " << prm.hdr_frm.n_columns << std::endl;
				fout << "Frame rows: " << prm.hdr_frm.n_rows << std::endl;
				fout << "Data integer bits: " << (int)prm.hdr_frm.n_bpi << std::endl;
				fout.close();
				if (prm.btalk) {
					std::cout << "- written output data info file " << prm.str_file_output + ".hdr" << std::endl;
				}
			}
			else {
				std::cerr << "Error: failed to open info file " << prm.str_file_output + ".hdr" << " for writing.\n";
				nerr = 200;
			}
		}
	}
	return nerr;
}



int run_average_frames()
{
	int nerr = 0;
	size_t i_frm = 0; // frame index
	size_t i_pix = 0; // pixel index
	size_t frm_pix = (size_t)prm.hdr_frm.n_columns * prm.hdr_frm.n_rows; // number of frame items
	size_t scan_pix = (size_t)prm.hdr.n_columns * prm.hdr.n_rows; // number of scan items
	double dtmp = 0.;
	double * datbuf = NULL; // pre-processed data buffer
	double * resbuf = NULL; // result buffer
	double * devbuf = NULL; // deviation buffer
	merlin_pix scan_pos; // scan position index
	std::string str_file; // file name
	std::ifstream fin; // input stream
	std::ofstream fout; // output stream
	std::streampos fpos; // file position
	int ncfidx = -1;
	int fidx = -1; // file index
	int file_frm = -1; // frame index in file
	int prog_pct = 0;
	int prog_pct_old = 0;
	size_t nres = 0; // number of result items
	if (frm_pix > 0 && prm.hdr.n_frames > 0) {
		datbuf = (double*)calloc(frm_pix, sizeof(double));
		resbuf = (double*)calloc(frm_pix, sizeof(double));
		devbuf = (double*)calloc(frm_pix, sizeof(double));
		if (prm.btalk) {
			std::cout << "- averaging frames in current scan roi ...\n";
			std::cout << "  0 %\r";
		}
		for (i_frm = 0; i_frm < prm.hdr.n_frames; i_frm++) { // loop over all frames
			if (0 == prm.get_scan_pixel((int)i_frm, scan_pos.x, scan_pos.y)) { // got a scan position for frame
				if (prm.in_scan_roi(scan_pos, prm.scan_rect_roi)) { // scan position is in ROI
					nerr = prm.get_frame_filepos((int)i_frm, fidx, fpos); // get the file index and position
					if (nerr != 0) {
						std::cerr << "Error: failed to determine file index and position for frame : " << i_frm << ".\n";
						nerr = 103;
						goto _cancel_point; // stop working
					}
					if (fidx != ncfidx && fidx >= 0) { // not the right file is open
						if (ncfidx >= 0 && fin.is_open()) { // close the wrong file
							fin.close();
						}
						str_file = prm.str_file_input + std::to_string(fidx + 1) + ".mib"; // file name construction
						fin.open(str_file, std::ios::binary); // open the right file
						if (fin.is_open()) { // successfully opened
							ncfidx = fidx; // update current file index
						}
						else { // failure opening file
							std::cerr << "Error: failed to open file: " << str_file << " for reading data.\n";
							nerr = 102;
							goto _cancel_point; // stop working
						}
					}
					if (fidx < 0) { // couldn't determine file index
						std::cerr << "Error: failed to determine file index for frame # " << i_frm << std::endl;
						nerr = 101;
						goto _cancel_point; // stop working
					}
					nerr = merlin_read_data(datbuf, fpos, &fin, &prm.hdr, &prm.hdr_frm, prm.swapbytes);
					if (nerr != 0) { // integration failed
						std::cerr << "Error: failed loading data of frame # " << i_frm << " (code " << nerr << ").\n";
						nerr = 106;
						goto _cancel_point; // stop working
					}
					// integrate to result buffer
					for (i_pix = 0; i_pix < frm_pix; i_pix++) {
						dtmp = datbuf[i_pix];
						resbuf[i_pix] += dtmp; // accumulate values
						devbuf[i_pix] += (dtmp*dtmp); // accumulate squares
					}
					nres++; // increment result numbers
				}
			}
			// 
			prog_pct = (int)(100. * (double)i_frm / (double)prm.hdr.n_frames); // progress in percent
			if (prm.btalk && prog_pct > prog_pct_old) { // progress step ...
				std::cout << "  " << prog_pct << " %\r";
				prog_pct_old = prog_pct;
			}
		} // frame loop
		if (nres > 0) { // normalize result buffer (otherwise we have 0 in the result)
			// check for required update of the defect correction list
			if (prm.is_defect_list_modified()) prm.update_defect_correction_list();
			//
			nerr = prm.gain_correction(resbuf); // apply gain correction if present
			if (nerr != 0) { // gain correction failed
				std::cerr << "Error: gain correction failed on accumulated data (code " << nerr << ").\n";
				nerr = 110;
				goto _cancel_point; // stop working
			}
			nerr = prm.defect_correction(resbuf); // apply defect pixel correction if present
			if (nerr != 0) { // defect correction failed
				std::cerr << "Error: defect correction failed on accumulated data (code " << nerr << ").\n";
				nerr = 111;
				goto _cancel_point; // stop working
			}
			nerr = prm.gain_correction(devbuf); // apply gain correction twice of the square accumulation
			nerr = prm.gain_correction(devbuf); // ...
			if (prm.ndebug > 0 && prm.btalk) {
				std::cout << "- rescaling output to average of " << nres << " frames.\n";
			}
			double sca = 1. / (double)nres;
			for (i_pix = 0; i_pix < frm_pix; i_pix++) {
				resbuf[i_pix] *= sca; // calculate mean
				dtmp = devbuf[i_pix] * sca;
				devbuf[i_pix] = sqrt(dtmp - resbuf[i_pix] * resbuf[i_pix]); // calculate std. deviation
			}
			nerr = prm.defect_correction(devbuf); // apply defect pixel correction on the devation image
		}
		else {
			std::cerr << "Error: averaging over zero frames.\n";
			nerr = 110;
		}
	_cancel_point:
		if (fin.is_open()) fin.close();
		if (datbuf) free(datbuf);
		if (nerr == 0) { // write the result to files
			str_file = prm.str_file_output + "_avg.dat";
			if (0 == write_data((char*)resbuf, sizeof(double)*frm_pix, str_file)) {
				if (prm.btalk) {
					std::cout << "- written average frame to file " << str_file << ".\n";
					std::cout << "  data type: floating point, 64 bit\n";
					std::cout << "  sampling: " << prm.hdr_frm.n_columns << " x " << prm.hdr_frm.n_rows << " scan points\n";
				}
			}
			else {
				nerr = 200;
			}
			str_file = prm.str_file_output + "_sdev.dat";
			if (0 == write_data((char*)devbuf, sizeof(double)*frm_pix, str_file)) {
				if (prm.btalk) {
					std::cout << "- written standard deviation frame to file " << str_file << ".\n";
					std::cout << "  data type: floating point, 64 bit\n";
					std::cout << "  sampling: " << prm.hdr_frm.n_columns << " x " << prm.hdr_frm.n_rows << " scan points\n";
				}
			}
			else {
				nerr = 210;
			}
		}
		if (resbuf) free(resbuf);
		if (devbuf) free(devbuf);
	}
	return nerr;
}



int run_integrate_annular_range()
{
	int nerr = 0;
	int i_frm = 0; // frame index
	size_t frm_pix = (size_t)prm.hdr_frm.n_columns * prm.hdr_frm.n_rows; // number of frame items
	size_t scan_pix = (size_t)prm.hdr.n_columns * prm.hdr.n_rows; // number of scan items
	int * dethash = NULL; // detector function hash
	double * datbuf = NULL; // pre-processed data buffer
	double * detbuf = NULL; // detector function buffer
	double * resbuf = NULL; // result buffer
	merlin_pix scan_pos; // scan position index
	std::string str_file; // file name
	std::ifstream fin; // input stream
	std::streampos fpos; // file position
	int ncfidx = -1;
	int fidx = -1; // file index
	int file_frm = -1; // frame index in file
	int prog_pct = 0;
	int prog_pct_old = 0;
	size_t nhash = 0; // number of hashed pixels
	size_t nres = 0; // number of result items
	if (frm_pix > 0 && prm.hdr.n_frames > 0 && prm.range_annular.max > prm.range_annular.min) {
		dethash = (int*)calloc(frm_pix, sizeof(int));
		datbuf = (double*)calloc(frm_pix, sizeof(double));
		detbuf = (double*)calloc(frm_pix, sizeof(double));
		resbuf = (double*)calloc(scan_pix, sizeof(double));
		nerr = prepare_annular_detector(frm_pix, detbuf, dethash, &nhash, &prm.hdr_frm);
		if (nerr != 0) {
			std::cerr << "Error: failed to prepare annular detector.\n";
			nerr = 10;
			goto _cancel_point;
		}
		// check for required update of the defect correction list
		if (prm.is_defect_list_modified()) prm.update_defect_correction_list();
		//
		if (prm.btalk) {
			std::cout << "- integration over annular range in current scan roi ...\n";
			std::cout << "  0 %\r";
		}
		for (i_frm = 0; i_frm < prm.hdr.n_frames; i_frm++) {
			nerr = prm.get_scan_pixel((int)i_frm, scan_pos.x, scan_pos.y);
			if (nerr != 0) { // integration failed
				std::cerr << "Error: failed to determine scan position for frame # " << i_frm << " (code " << nerr << ").\n";
				nerr = 100;
				goto _cancel_point; // stop working
			}
			if (prm.in_scan_roi(scan_pos, prm.scan_rect_roi)) {
				nerr = prm.get_frame_filepos(i_frm, fidx, fpos);
				if (nerr != 0) { // couldn't determine file index
					std::cerr << "Error: failed to determine file index and datat offset for frame # " << i_frm << std::endl;
					nerr = 101;
					goto _cancel_point; // stop working
				}
				if (fidx != ncfidx && fidx >= 0) { // not the right file is open
					if (ncfidx >= 0 && fin.is_open()) { // close the wrong file
						fin.close();
					}
					str_file = prm.str_file_input + std::to_string(fidx + 1) + ".mib"; // file name construction
					fin.open(str_file, std::ios::binary); // open the right file
					if (fin.is_open()) { // successfully opened
						ncfidx = fidx; // update current file index
					}
					else { // failure opening file
						std::cerr << "Error: failed to open file: " << str_file << " for reading data.\n";
						nerr = 102;
						goto _cancel_point; // stop working
					}
				}
				nerr = merlin_read_data(datbuf, fpos, &fin, &prm.hdr, &prm.hdr_frm, prm.swapbytes);
				if (nerr != 0) { // integration failed
					std::cerr << "Error: failed loading data of frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 106;
					goto _cancel_point; // stop working
				}
				nerr = prm.gain_correction(datbuf); // apply gain correction if present
				if (nerr != 0) { // gain correction failed
					std::cerr << "Error: gain correction failed for frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 110;
					goto _cancel_point; // stop working
				}
				nerr = prm.defect_correction(datbuf); // apply defect pixel correction if present
				if (nerr != 0) { // defect correction failed
					std::cerr << "Error: defect correction failed for frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 111;
					goto _cancel_point; // stop working
				}
				nerr = sum_annular_range(frm_pix, datbuf, detbuf, dethash, nhash, &resbuf[nres]);
				if (nerr != 0) { // integration failed
					std::cerr << "Error: detector readout failed for frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 112;
					goto _cancel_point; // stop working
				}
				nres++; // increment result numbers
			} // if in roi
			// 
			prog_pct = (int)(100. * (double)i_frm / (double)prm.hdr.n_frames); // progress in percent
			if (prm.btalk && prog_pct > prog_pct_old) { // progress step ...
				std::cout << "  " << prog_pct << " %\r";
				prog_pct_old = prog_pct;
			}
		} // frame loop
	_cancel_point:
		if (fin.is_open()) fin.close();
		if (prm.ndebug > 0) {
			if (0 == write_data((char*)detbuf, sizeof(double)*frm_pix, prm.str_file_output + ".det")) {
				std::cout << "- written detector function to file " << prm.str_file_output + ".det" << ".\n";
			}
		}
		if (dethash) free(dethash);
		if (detbuf) free(detbuf);
		if (datbuf) free(datbuf);
		if (nerr == 0 && nres == 0) {
			if (prm.btalk) {
				std::cout << "No results calculated, output skipped.\n";
			}
		}
		if (nerr == 0 && nres > 0 && resbuf) { // write the result to file
			if (0 == write_data((char*)resbuf, sizeof(double)*nres, prm.str_file_output)) {
				if (prm.btalk) {
					std::cout << "- written integrated annular range data to file " << prm.str_file_output << ".\n";
					std::cout << "  data type: floating point, 64 bit\n";
					std::cout << "  scan sampling: " << 1 + prm.scan_rect_roi.x1 - prm.scan_rect_roi.x0 << " x " << 1 + prm.scan_rect_roi.y1 - prm.scan_rect_roi.y0 << " scan points\n";
				}
			}
			else {
				nerr = 200;
			}
		}
		if (resbuf) free(resbuf);
	}
	else {
		if (frm_pix == 0) {
			std::cerr << "Error: insuffient number of frame pixels.\n";
			nerr = 1;
		}
		if (prm.hdr.n_frames <= 0) {
			std::cerr << "Error: insuffient number of frames.\n";
			nerr = 2;
		}
		if (prm.range_annular.max <= prm.range_annular.min) {
			std::cerr << "Error: invalid annular range (" << prm.range_annular.min << "," << prm.range_annular.max << ").\n";
			nerr = 3;
		}
	}
	return nerr;
}


int run_center_of_mass()
{
	int nerr = 0;
	int i_frm = 0; // frame index
	size_t frm_pix = (size_t)prm.hdr_frm.n_columns * prm.hdr_frm.n_rows; // number of frame items
	size_t scan_pix = (size_t)prm.hdr.n_columns * prm.hdr.n_rows; // number of scan items
	int * dethash = NULL; // detector function hash
	double * datbuf = NULL; // pre-processed data buffer
	double * detbuf = NULL; // detector function buffer
	double * resbuf00 = NULL; // result buffer - integral
	double * resbuf10 = NULL; // result buffer - com.x
	double * resbuf11 = NULL; // result buffer - com.y
	double * xbuf = NULL; // x-coordinates of the data frame
	double * ybuf = NULL; // y-coordinates of the data frame
	merlin_pix scan_pos; // scan position index
	std::string str_file; // file name
	std::string str_file_out; // file names for output
	std::ifstream fin; // input stream
	std::streampos fpos; // file position
	int ncfidx = -1;
	int fidx = -1; // file index
	int file_frm = -1; // frame index in file
	int prog_pct = 0;
	int prog_pct_old = 0;
	size_t nhash = 0; // number of hashed pixels
	size_t nres = 0; // number of result items
	if (frm_pix > 0 && prm.hdr.n_frames > 0 && prm.range_annular.max > prm.range_annular.min) {
		dethash = (int*)calloc(frm_pix, sizeof(int));
		datbuf = (double*)calloc(frm_pix, sizeof(double));
		detbuf = (double*)calloc(frm_pix, sizeof(double));
		xbuf = (double*)calloc(frm_pix, sizeof(double));
		ybuf = (double*)calloc(frm_pix, sizeof(double));
		resbuf00 = (double*)calloc(scan_pix, sizeof(double));
		resbuf10 = (double*)calloc(scan_pix, sizeof(double));
		resbuf11 = (double*)calloc(scan_pix, sizeof(double));
		nerr = prepare_annular_detector(frm_pix, detbuf, dethash, &nhash, &prm.hdr_frm);
		if (nerr != 0) {
			std::cerr << "Error: failed to prepare annular detector.\n";
			goto _cancel_point;
		}
		nerr = prepare_frame_coordinates(frm_pix, xbuf, ybuf, &prm.hdr_frm);
		if (nerr != 0) {
			std::cerr << "Error: failed to prepare frame coordinates.\n";
			goto _cancel_point;
		}
		// check for required update of the defect correction list
		if (prm.is_defect_list_modified()) prm.update_defect_correction_list();
		//
		if (prm.btalk) {
			std::cout << "- integration over annular range in current scan roi ...\n";
			std::cout << "  0 %\r";
		}
		for (i_frm = 0; i_frm < prm.hdr.n_frames; i_frm++) {
			nerr = prm.get_scan_pixel((int)i_frm, scan_pos.x, scan_pos.y);
			if (nerr != 0) { // integration failed
				std::cerr << "Error: failed to determine scan position for frame # " << i_frm << " (code " << nerr << ").\n";
				nerr = 100;
				goto _cancel_point; // stop working
			}
			if (prm.in_scan_roi(scan_pos, prm.scan_rect_roi)) {
				nerr = prm.get_frame_filepos(i_frm, fidx, fpos);
				if (nerr != 0) { // couldn't determine file index
					std::cerr << "Error: failed to determine file index and datat offset for frame # " << i_frm << std::endl;
					nerr = 101;
					goto _cancel_point; // stop working
				}
				if (fidx != ncfidx && fidx >= 0) { // not the right file is open
					if (ncfidx >= 0 && fin.is_open()) { // close the wrong file
						fin.close();
					}
					str_file = prm.str_file_input + std::to_string(fidx + 1) + ".mib"; // file name construction
					fin.open(str_file, std::ios::binary); // open the right file
					if (fin.is_open()) { // successfully opened
						ncfidx = fidx; // update current file index
					}
					else { // failure opening file
						std::cerr << "Error: failed to open file: " << str_file << " for reading data.\n";
						nerr = 102;
						goto _cancel_point; // stop working
					}
				}
				nerr = merlin_read_data(datbuf, fpos, &fin, &prm.hdr, &prm.hdr_frm, prm.swapbytes);
				if (nerr != 0) { // integration failed
					std::cerr << "Error: failed loading data of frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 106;
					goto _cancel_point; // stop working
				}
				nerr = prm.gain_correction(datbuf); // apply gain correction if present
				if (nerr != 0) { // gain correction failed
					std::cerr << "Error: gain correction failed for frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 110;
					goto _cancel_point; // stop working
				}
				nerr = prm.defect_correction(datbuf); // apply defect pixel correction if present
				if (nerr != 0) { // defect correction failed
					std::cerr << "Error: defect correction failed for frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 111;
					goto _cancel_point; // stop working
				}
				nerr = sum_annular_range(frm_pix, datbuf, detbuf, dethash, nhash, &resbuf00[nres]);
				if (nerr != 0) { // integration failed
					std::cerr << "Error: detector readout failed for frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 112;
					goto _cancel_point; // stop working
				}
				nerr = com_annular_range(frm_pix, datbuf, detbuf, xbuf, ybuf, dethash, nhash, resbuf00[nres], &resbuf10[nres], &resbuf11[nres]);
				if (nerr != 0) { // integration failed
					std::cerr << "Error: detector readout failed for frame # " << i_frm << " (code " << nerr << ").\n";
					nerr = 112;
					goto _cancel_point; // stop working
				}
				nres++; // increment result numbers
			} // if in roi
			// 
			prog_pct = (int)(100. * (double)i_frm / (double)prm.hdr.n_frames); // progress in percent
			if (prm.btalk && prog_pct > prog_pct_old) { // progress step ...
				std::cout << "  " << prog_pct << " %\r";
				prog_pct_old = prog_pct;
			}
		} // frame loop
	_cancel_point:
		if (fin.is_open()) fin.close();
		if (prm.ndebug > 0) {
			if (0 == write_data((char*)detbuf, sizeof(double)*frm_pix, prm.str_file_output + ".det")) {
				std::cout << "- written detector function to file " << prm.str_file_output + ".det" << ".\n";
			}
		}
		if (dethash) free(dethash);
		if (detbuf) free(detbuf);
		if (datbuf) free(datbuf);
		if (xbuf) free(xbuf);
		if (ybuf) free(ybuf);
		if (nerr == 0 && nres == 0) {
			if (prm.btalk) {
				std::cout << "No results calculated, output skipped.\n";
			}
		}
		if (nerr == 0 && nres > 0 && resbuf00 && resbuf10 && resbuf11) { // write the results to files
			// - reference integrals 0-0
			str_file_out = prm.str_file_output + "_0-0.dat";
			if (0 == write_data((char*)resbuf00, sizeof(double)*nres, str_file_out)) {
				if (prm.btalk) {
					std::cout << "- written reference integrals to file " << str_file_out << ".\n";
					std::cout << "  data type: floating point, 64 bit\n";
					std::cout << "  scan sampling: " << 1 + prm.scan_rect_roi.x1 - prm.scan_rect_roi.x0 << " x " << 1 + prm.scan_rect_roi.y1 - prm.scan_rect_roi.y0 << " scan points\n";
				}
			}
			else {
				nerr = 200;
			}
			// - center-of-mass x 1-0
			str_file_out = prm.str_file_output + "_1-0.dat";
			if (0 == write_data((char*)resbuf10, sizeof(double)*nres, str_file_out)) {
				if (prm.btalk) {
					std::cout << "- written center-of-mass x to file " << str_file_out << ".\n";
					std::cout << "  data type: floating point, 64 bit\n";
					std::cout << "  scan sampling: " << 1 + prm.scan_rect_roi.x1 - prm.scan_rect_roi.x0 << " x " << 1 + prm.scan_rect_roi.y1 - prm.scan_rect_roi.y0 << " scan points\n";
				}
			}
			else {
				nerr = 210;
			}
			// - center-of-mass x 1-1
			str_file_out = prm.str_file_output + "_1-1.dat";
			if (0 == write_data((char*)resbuf11, sizeof(double)*nres, str_file_out)) {
				if (prm.btalk) {
					std::cout << "- written center-of-mass y to file " << str_file_out << ".\n";
					std::cout << "  data type: floating point, 64 bit\n";
					std::cout << "  scan sampling: " << 1 + prm.scan_rect_roi.x1 - prm.scan_rect_roi.x0 << " x " << 1 + prm.scan_rect_roi.y1 - prm.scan_rect_roi.y0 << " scan points\n";
				}
			}
			else {
				nerr = 220;
			}
		}
		if (resbuf00) free(resbuf00);
		if (resbuf10) free(resbuf10);
		if (resbuf11) free(resbuf11);
	}
	return nerr;
}





// -----------------------------------------------------------------------------
//
// CONTROL interface
//
// -----------------------------------------------------------------------------


int input_getline(std::istream * pfin, std::string * str) {
	if (NULL == pfin) {
		return 1; // missing or invalid parameter 1
	}
	if (NULL == str) {
		return 2; // missing or invalid parameter 2
	}
	str->clear();
	if (!pfin->good()) {
		return 11; // input stream is not good
	}
	getline(*pfin, *str);
	if (!pfin->good()) {
		std::cerr << "Error reading command from input.\n";
		return 100;
	}
	return 0;
}

int file_getline(std::ifstream * pfin, std::string * str) {
	if (NULL == pfin) {
		return 1; // missing or invalid parameter 1
	}
	if (NULL == str) {
		return 2; // missing or invalid parameter 2
	}
	str->clear();
	if (!pfin->good()) {
		return 11; // input stream is not good
	}
	getline(*pfin, *str);
	if (pfin->fail() && !pfin->eof()) {
		std::cerr << "Error reading command from input file.\n";
		return 100;
	}
	return 0;
}


// reads a new command or parameter line from std::cin or the control
// file lines to (str) and increments (iline)
// (inc) is a prefix shouted by the interactive input console
int ctrl_getline(size_t & iline, std::string inc, std::string * str)
{
	bool bsuccess = false;
	int nerr = 0;
	if (NULL == str) {
		return 1;
	}
	if (prm.binteractive) { // interactive line input
		std::cout << inc << " > "; // prefix
		nerr = input_getline(&(std::cin), str); // get input from std::cin to string
		if (nerr == 0) { // success
			bsuccess = true;
			prm.v_str_ctrl.push_back(*str);
			iline++;
		}
	}
	else {
		if (iline >= 0 && iline < prm.v_str_ctrl.size()) {
			*str = prm.v_str_ctrl[iline];
			bsuccess = true;
			iline++;
			if (prm.btalk) {
				std::cout << *str << std::endl;
			}
		}
	}
	if (!bsuccess && nerr == 0) {
		nerr = 300;
	}
	return nerr;
}


int read_ctrl_file(std::string * str_ctrl = NULL)
{
	int nerr = 0;
	std::ifstream fcin;
	std::string sline;
	std::string str_file = prm.str_file_ctrl;
	if (str_ctrl) str_file = *str_ctrl; // use file name from routine parameter
	fcin.open(str_file);
	if (!fcin.is_open()) {
		std::cerr << "Error: failed to open control file " << str_file << " for reading." << std::endl;
		return 1;
	}
	prm.v_str_ctrl.clear();
	while (nerr == 0) {
		nerr = file_getline(&fcin, &sline);
		if (nerr == 0) {
			prm.v_str_ctrl.push_back(sline);
		}
		else {
			break;
		}
	}
	if (nerr != 0 && fcin.eof()) { // no reading error, just eof
		nerr = 0;
	}
	if (nerr != 0) {
		std::cerr << "Error reading command from file.\n ";
	}
	if (fcin.is_open()) fcin.close();
	return nerr;
}

int write_ctrl_file(std::string * str_ctrl = NULL)
{
	int nerr = 0;
	size_t nlines = 0, iline = 0;
	std::ofstream fcout;
	std::string sline;
	std::string str_file = prm.str_file_ctrl;
	if (str_ctrl) str_file = *str_ctrl; // use file name from routine parameter
	fcout.open(str_file, std::ios::trunc);
	if (!fcout.is_open()) {
		std::cerr << "Error: failed to open control file " << str_file << " for writing." << std::endl;
		return 1;
	}
	nlines = prm.v_str_ctrl.size();
	if (nlines > 0) {
		for (iline = 0; iline < nlines; iline++) {
			fcout << prm.v_str_ctrl[iline] << std::endl;
		}
	}
	if (fcout.is_open()) fcout.close();
	if (nerr == 0 && prm.btalk) {
		std::cout << "- command list written to file: " << str_file << std::endl;
	}
	return nerr;
}

int run_ctrl(std::string file_ctrl)
{
	int nerr = 0;
	size_t ncmd = 0, icmd = 0;
	bool bnextcommand = true;
	bool bprocessed = false;
	std::string scmd, sprm;
	std::ifstream fcin;
	std::ofstream fcout;
	
	// handle control file I/O depending on interactive mode
	if (prm.binteractive) {
		prm.v_str_ctrl.clear();
	}
	else {
		nerr = read_ctrl_file(&file_ctrl);
		if (nerr != 0) {
			std::cerr << "Error: failed to read control file " << file_ctrl << std::endl;
			return 1;
		}
		ncmd = prm.v_str_ctrl.size();
		if (ncmd == 0) {
			std::cerr << "Error: empty control file " << file_ctrl << std::endl;
			bnextcommand = false; // no commands found in file
			return 2;
		}
	}

	while (bnextcommand) {

		bprocessed = false;

		// get command
		nerr = ctrl_getline(icmd, "", &scmd);
		if (nerr != 0) {
			std::cerr << "Error reading command.\n ";
			break;
		}

		// transform command to lower case for easier identification
		std::transform(scmd.begin(), scmd.end(), scmd.begin(), ::tolower);

		// identify command and react
		//
		// - set commands
		if (scmd == "set_scan_rect_roi") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) nerr = prm.set_scan_rect_roi(sprm);
			bprocessed = true;
		}

		if (scmd == "set_origin") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) nerr = prm.set_origin(sprm);
			bprocessed = true;
		}

		if (scmd == "set_sampling") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) nerr = prm.set_sampling(sprm);
			bprocessed = true;
		}

		if (scmd == "set_annular_range") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) nerr = prm.set_annular_range(sprm);
			bprocessed = true;
		}

		if (scmd == "set_output_file") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) prm.str_file_output = sprm;
			bprocessed = true;
		}

		if (scmd == "set_defect_mask") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) {
				nerr = prm.load_defect_mask(sprm);
			}
			bprocessed = true;
		}

		if (scmd == "set_defect_list") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) {
				nerr = prm.load_defect_list(sprm);
			}
			bprocessed = true;
		}

		if (scmd == "set_defect_pixel") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) {
				nerr = prm.set_defect_pixel(sprm);
			}
			bprocessed = true;
		}

		if (scmd == "unset_defect_pixel") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) {
				nerr = prm.unset_defect_pixel(sprm);
			}
			bprocessed = true;
		}

		if (scmd == "unset_defect_list") {
			nerr = prm.unset_defect_list();
			bprocessed = true;
		}

		if (scmd == "set_gain_correction") {
			nerr = ctrl_getline(icmd, scmd, &sprm); // get parameter input
			if (nerr == 0) {
				nerr = prm.load_gain_correction(sprm);
			}
			bprocessed = true;
		}

		if (scmd == "unset_gain_correction") {
			nerr = prm.unset_gain_correction();
			bprocessed = true;
		}

		//
		// - processing commands
		if (scmd == "extract_frames") {
			nerr = run_extract_frames();
			bprocessed = true;
		}

		if (scmd == "average_frames") {
			nerr = run_average_frames();
			bprocessed = true;
		}

		if (scmd == "integrate_annular_range") {
			nerr = run_integrate_annular_range();
			bprocessed = true;
		}

		if (scmd == "center_of_mass") {
			nerr = run_center_of_mass();
			bprocessed = true;
		}

		if (scmd == "exit" || scmd == "quit") {
			if (prm.btalk) {
				std::cout << "Exiting program.\n";
			}
			bnextcommand = false;
			bprocessed = true;
		}

		if (!bprocessed && scmd.size() > 0) {
			std::cerr << "Unknown or invalid command: " << scmd << std::endl;
		}
		if (bprocessed && nerr != 0) {
			std::cerr << "Error while processing command: " << scmd << " (code: " << nerr << ")\n";
		}

	}

	if (prm.binteractive) {
		std::cout << "Do you want to store the commands to a new control file? <1> Yes. <2> No. ";
		std::cin >> scmd;
		if (scmd == "1") {
			nerr = write_ctrl_file(&file_ctrl);
			if (nerr != 0) {
				std::cerr << "Error: failed to write command strings to control file " << file_ctrl << std::endl;
				return 3;
			}
		}
	}

	return 0;
}


int control_interface(void)
{
	bool bfilefound = false;
	prm.binteractive = true;
	struct stat statbuf;

	if (prm.btalk) {
		std::cout << std::endl;
	}

	bfilefound = (stat(prm.str_file_ctrl.c_str(), &statbuf) == 0); // check if file exists
	if (bfilefound) {
		prm.binteractive = false;
		if (prm.btalk) {
			std::cout << "Running control file: " << prm.str_file_ctrl << std::endl;
		}
	}
	else {
		prm.binteractive = true;
		if (prm.btalk) {
			std::cout << "Interactive control input" << std::endl;
		}
	}
	return run_ctrl(prm.str_file_ctrl);
}



// -----------------------------------------------------------------------------
//
// MAIN
//
// -----------------------------------------------------------------------------


int main(int argc, char* argv[])
{
	// preset error code
	int nerr = 0;
	unsigned __int16 testendian;
	char bl[2];
	

	nerr = parseoptions(argc, argv);
	if (0 < nerr) {
		std::cerr << "Error while parsing call options (code " << nerr << ").\n";
		return 1;
	}
		
	if (prm.btalk) {
		std::cout << "Running program MERLINIO\n";
		std::cout << "  "<< MERLINIO_VER <<"."<< MERLINIO_VER_SUB << "." << MERLINIO_VER_SUBSUB << " ("<< MERLINIO_VER_BUILD << ")\n";
		std::cout << "  by J. Barthel, Copyright (c) 2019\n";
		std::cout << "  ju.barthel@fz-juelich.de\n";
		std::cout << "  Forschungszentrum Juelich GmbH, Juelich, Germany\n";
		std::cout << "\n";

		
		if (prm.ndebug > 0) {
			std::cout << "- running in debug mode (level " << prm.ndebug << ")\n";
			std::cout << "- control file: " << prm.str_file_ctrl << std::endl;
			if (prm.bscanframeheaders) std::cout << "- scanning frame headers.\n";
		}
		std::cout << "- input files: " << prm.str_file_input << std::endl;
		std::cout << "- output files: " << prm.str_file_output << std::endl;
		bl[0] = 1; bl[1] = 0;
		testendian = *((unsigned __int16*)bl);
		if (testendian == 1) {
			std::cout << "- I'm little endian. Assuming that Merlin delivers big endian: swapping bytes.\n";
			prm.swapbytes = true;
		}
	}
	
	nerr = prm.read_header();
	if (0 < nerr) {
		std::cerr << "Error while reading the header file (code " << nerr << ").\n";
		return 2;
	}

	nerr = prm.read_frame_headers();
	if (0 < nerr) {
		std::cerr << "Error while reading the data files (code " << nerr << ").\n";
		return 3;
	}

	// preset scan roi again, now that we know the frame size
	prm.scan_rect_roi.x0 = 0;
	prm.scan_rect_roi.y0 = 0;
	prm.scan_rect_roi.x1 = prm.hdr.n_columns - 1;
	prm.scan_rect_roi.y1 = prm.hdr.n_rows - 1;

	nerr = control_interface();
	if (0 < nerr) {
		std::cerr << "Error in the control interface (code " << nerr << ").\n";
		return 4;
	}

    
	if (prm.btalk) {
		std::cout << "\n";
		std::cout << "Done.\n";
	}

	return nerr;
}


