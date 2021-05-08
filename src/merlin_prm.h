// file : "merlin_prm.h"
// author: J. Barthel, Forschungszentrum Juelich GmbH, Juelich, Germany
//         ju.barthel@fz-juelich.de
//
// Declares parameter structures used by merlinio
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
#include "merlin_hdr.h"

constexpr auto MERLINIO_VER = 1;
constexpr auto MERLINIO_VER_SUB = 1;
constexpr auto MERLINIO_VER_SUBSUB = 0;
constexpr auto MERLINIO_VER_BUILD = 3105081359;

struct defect_pixel_corr {
	size_t idx;
	int x;
	int y;
	std::vector<size_t> v_idx_corr;
};

class merlin_params
{
public:
	// constructor
	merlin_params();
	// destructor
	~merlin_params();
	// member parameters
	bool btalk; // flag talkative mode
	bool binteractive; // flag interactive control mode
	bool bscanframeheaders; // flag causing a careful frame header scan
	bool swapbytes; // flag for swapping bytes when converting to floats
	int ndebug; // debug level
	
	merlin_hdr hdr;
	merlin_frame_hdr hdr_frm;
	std::vector<int> v_frm_file;
	std::vector<std::streampos> v_frm_pos;
	
	merlin_frame_calib frame_calib;
	merlin_range range_annular;
	merlin_pos offset_annular;
	merlin_roi scan_rect_roi;
	std::string str_file_input;
	std::string str_file_output;
	std::string str_file_ctrl;

	std::vector<std::string> v_str_ctrl; // list of commands read from control file

protected:
	bool gaincorrect; // indicates that a gain correction is present and used
	bool defects_modified; // indicates that the defect list has been modified and that the correction lists need updates
	std::vector<defect_pixel_corr> v_defect_corr; // list of registered defect pixels with correction data
	double * img_gaincorrect; // gain correction factors (size determine by hdr_frm)
	
	// member functions
public:
	// separates a single parameter string (*prm)
	// from a list of parameters in a string (*pstr)
	// beginning with position (ipos) of (*pstr)
	int read_param(int ipos, std::string * pstr, std::string * prm);
	
	// reads information from a merlin header file ".hdr"
	int read_header(void);
	// reads information from merlin frame headers in merlin ".mib" files.
	int read_frame_headers(void);

	
	// applies the given frame calibration
	//   dx = xin.x - pcalib->offset.x
	//   dy = xin.y - pcalib->offset.y
	//   xout->x = dx * pcalib->a0.x + dy * pcalib->a1.x
	//   xout->y = dx * pcalib->a0.y + dy * pcalib->a1.y
	// returns 0 in case of success
	int get_calib_pos(merlin_pix xin, merlin_pos * xout);

	// determines the scan x,y position from the frame index
	// - input idx = frame index
	// - output x, y = scan position
	// - return value = error code (0: success)
	int get_scan_pixel(int idx, int &x, int &y);

	// determines the frame index from scan position x,y
	// - input x, y = scan position
	// - return value = frame index (<0 = error code)
	int get_frame_idx(int x, int y);

	// determines the frame pixel x,y position from the frame pixel index
	// - input idx = frame pixel index
	// - output x, y = frame pixel position
	// - return value = error code (0: success)
	int get_frame_pixel(int idx, int &x, int &y);

	// determines the frame pixel index from the frame pixel x,y position 
	// - input x, y = frame pixel position
	// - return value = frame pixel index (<0: error)
	int get_frame_pixel_idx(int x, int y);

	// determines the file index (ifile) and the file position (fpos)
	// for a given global frame index (idx)
	// - input idx = global frame index (0 = first frame)
	// - output ifile = file index (0 = first file)
	// - output ipos  = position of frame in file (0 = begin of file)
	// - return value = error code (0: success)
	int get_frame_filepos(int idx, int &ifile, std::streampos &ipos);

	// returns the number of pixels in the current rectangular scan roi
	size_t get_scan_rect_roi_size(void);


	bool in_scan_roi(merlin_pix pos, merlin_roi roi);


	int set_scan_rect_roi(std::string str_roi);

	int set_origin(std::string str_org);

	int	set_sampling(std::string str_samp);

	int set_annular_range(std::string str_rng);

	int set_annular_offset(std::string str_pos);

	bool is_defect_pixel(size_t idx);
	
	bool is_defect_pixel(int x, int y);

	bool is_defect_list_modified(void);

	// writes new defect correction tables
	int update_defect_correction_list(void);

	int set_defect_pixel(int x, int y);

	// adds the pixel at x,y frame coordinates to the defect list
	// the coordinates x,y are given as parameter string (str_pos)
	int set_defect_pixel(std::string str_pos);

	// loads a defect mask from file
	int load_defect_mask(std::string str_file);

	// loads a defect list from file
	int load_defect_list(std::string str_file);

	// unsets and frees all memory related to defect pixel correction
	int unset_defect_list(void);

	// removes the pixel at given x,y frame coordinates from the defect list
	int unset_defect_pixel(int x, int y);

	// removes the pixel at x,y frame coordinates from the defect list
	// the coordinates x,y are given as parameter string (str_pos)
	int unset_defect_pixel(std::string str_pos);

	// loads a gain correction image from file
	int load_gain_correction(std::string str_file);

	// unsets and frees all memory related to gain correction
	int unset_gain_correction(void);

	// applies the gain correction to frame data
	int gain_correction(double * buf);

	// applies the defect pixel correction to frame data
	int defect_correction(double * buf);
};