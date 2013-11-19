/**
 * @file  imageio.cxx
 * @brief Image I/O functions.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <iostream>            // cerr, endl
#include <fstream>             // ifstream
#include <cstdlib>             // abort()
#include <limits>              // numeric_limits<T>::min(), numeric_limits<T>::max()

#define _NIFTI1_IO_C_          // to get definition of LSB_FIRST
#include <nifti1_io.h>         // NIfTI I/O

#include <dramms/basis.h>      // BASIS_THROW, assert(), makedirs(), dirname(), hasext(), RELEASE

#include "general.h"           // round()
#include "utilities.h"         // CastImage(), ScaleImage(), ResizeImage()

#include "imageio.h"


// acceptable in .cxx file
using namespace std;
using namespace basis;


namespace dramms {


// NOTE: Function names in all lowercase and underscore (_) as word separator
//       denote functions that are only defined in this module and which are
//       not exposed to users of this image I/O module in the header file.
//
//       Exception: convert_mhd2nhdr() - to be compatible to nifti_convert_*().


// ===========================================================================
// image - either NIfTI-1, ANALYZE 7.5, or Meta
// ===========================================================================

// ---------------------------------------------------------------------------
template <typename T>
static bool read_image_data(znzFile fp, nifti_image* const nim, T*** p, bool merge)
{
    const int N  = (nim->nu > 0) ? nim->nu : 1; // components per voxel
    const int nz = (nim->nz > 0) ? nim->nz : 1; // number of slices, can be zero for 2D images
    if (merge) {
        const size_t nbytes = nim->nx * nim->nbyper; // bytes per row
        T* row = new T[nbytes / sizeof(T)];
        for (int n = 0; n < N; n++) {
            for (int k = 0; k < nz; k++) {
                for (int j = 0; j < nim->ny; j++) {
                    if (nifti_read_buffer(fp, row, nbytes, nim) != nbytes) {
                        return false;
                    }
                    for (int i = 0; i < nim->nx; i++) {
                        p[k][j][i * N + n] = row[i];
                    }
                }
            }
        }
        delete [] row;
    } else {
        const size_t nbytes = N * nim->nx * nim->nbyper; // bytes per row
        for (int k = 0; k < nz; k++) {
            for (int j = 0; j < nim->ny; j++) {
                if (nifti_read_buffer(fp, p[k][j], nbytes, nim) != nbytes) {
                    return false;
                }
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
template <typename T>
static bool write_image_data(znzFile fp, const nifti_image* nim, T*** const p, const Image* m, bool split)
{
    const int N  = (nim->nu > 0) ? nim->nu : 1; // components per voxel
    const int nz = (nim->nz > 0) ? nim->nz : 1; // number of slices, can be zero for 2D images
    T v;
    if (split) {
        for (int n = 0; n < N; n++) {
            for (int z = 0; z < nz; z++) {
                for (int y = 0; y < nim->ny; y++) {
                    for (int x = 0; x < nim->nx; x++) {
                        if (!m || m->img.uc[z][y][x]) {
                            v = p[z][y][x + N * n];
                        } else {
                            v = 0;
                        }
                        if (nifti_write_buffer(fp, &v, sizeof(T)) != sizeof(T)) {
                            return false;
                        }
                    }
                }
            }
        }
    } else {
        const size_t nbytes = N * nim->nx * sizeof(T); // bytes per row
        for (int z = 0; z < nz; z++) {
            for (int y = 0; y < nim->ny; y++) {
                if (m) {
                    for (int x = 0; x < nim->nx; x++) {
                        for (int n = 0; n < N; n++) {
                            if (m->img.uc[z][y][x]) {
                                v = p[z][y][x + N * n];
                            } else {
                                v = 0;
                            }
                            if (nifti_write_buffer(fp, &v, sizeof(T)) != sizeof(T)) {
                                return false;
                            }
                        }
                    }
                } else {
                    if (nifti_write_buffer(fp, p[z][y], nbytes) != nbytes) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
bool ReadRawImage(const char* filename, Image* image)
{
    assert(image != NULL);
    // open file
    znzFile fp = znzopen(filename, "rb", image->compress);
    if (znz_isnull(fp)) return false;
    // convert NIfTI-1 header to nifti_image structure
    nifti_1_header nhdr = image->hdr;
    if (image->imgfmt == Image::FORMAT_DRAMMS) nhdr = PermuteXY(nhdr);
    nifti_image *nim = nifti_convert_nhdr2nim(nhdr, NULL);
    if (!nim) return false;
    // read image data
    switch (image->hdr.datatype) {
        case DT_UNSIGNED_CHAR:
            read_image_data(fp, nim, image->img.uc, image->imgfmt == Image::FORMAT_FSL);
            break;
        case DT_SIGNED_SHORT:
            read_image_data(fp, nim, image->img.ss, image->imgfmt == Image::FORMAT_FSL);
            break;
        case DT_UINT16:
            read_image_data(fp, nim, image->img.us, image->imgfmt == Image::FORMAT_FSL);
            break;
        case DT_SIGNED_INT:
            read_image_data(fp, nim, image->img.si, image->imgfmt == Image::FORMAT_FSL);
            break;
        case DT_FLOAT:
            read_image_data(fp, nim, image->img.fl, image->imgfmt == Image::FORMAT_FSL);
            break;
        default:
            DRAMMS_ERROR("Cannot read raw image data of type: " << image->hdr.datatype);
            return false;
    }
    nifti_image_free(nim);
    znzclose(fp);
    return true;
}

// ---------------------------------------------------------------------------
bool WriteRawImage(const char* filename, const Image* image, const Image* mask)
{
    assert(image != NULL);
    if (mask && mask->hdr.datatype != DT_UNSIGNED_CHAR) {
        BASIS_THROW(invalid_argument, "Mask image must be of datatype unsigned char!");
    }
    if (mask && !mask->HasSameSizeAs(image)) {
        BASIS_THROW(invalid_argument, "Mask image has not the same size as the image to write!");
    }
    // -----------------------------------------------------------------------
    // If image data covers only parts of the image domain, resize the image,
    // i.e., "resample" it such that the image data is defined on the grid
    // specified by the image header. See also ResizeImage().
    if (image->region.ox != 0 ||
            image->region.oy != 0 ||
            image->region.oz != 0 ||
            image->region.nx != image->hdr.dim[1] ||
            image->region.ny != image->hdr.dim[2] ||
            image->region.nz != image->hdr.dim[3]) {
        Image* tmp = ResizeImage(image);
        if (tmp == NULL) {
            DRAMMS_ERROR("Failed to resize image to its original size before writing it!");
            return false;
        }
        Image* tmp_mask = NULL;
        if (mask) {
            tmp_mask = ResizeImage(mask);
            if (tmp_mask == NULL) {
                DRAMMS_ERROR("Failed to resize mask to its original size before writing the image!");
                delete tmp;
                return false;
            }
        }
        bool ok = WriteRawImage(filename, tmp, tmp_mask);
        delete tmp;
        if (tmp_mask) delete tmp_mask;
        return ok;
    }
    // -----------------------------------------------------------------------
    // write image data
    nifti_1_header hdr = image->hdr;
    if (image->imgfmt == Image::FORMAT_DRAMMS) hdr = PermuteXY(hdr);
    nifti_image *nim = nifti_convert_nhdr2nim(hdr, NULL);
    if (!nim) return false;
    // create output directory if not existent
    if (!os::makedirs(os::path::dirname(filename))) {
        nifti_image_free(nim);
        return false;
    }
    // open file
    string iname    = filename;
    string ext      = os::path::splitext(iname)[1];
    bool   compress = image->compress;
    if (ext == ".raw" || ext == ".img") compress = false;
    if (ext == ".gz")                   compress = true;
    if (compress && ext != ".gz") iname += ".gz";
    znzFile fp = znzopen(iname.c_str(), const_cast<char*>("wb"), image->compress);
    if (znz_isnull(fp)) {
        nifti_image_free(nim);
        return false;
    }
    // write image data
    switch (image->hdr.datatype) {
        case DT_UNSIGNED_CHAR:
            write_image_data(fp, nim, image->img.uc, mask, image->imgfmt == Image::FORMAT_FSL);
            break;
        case DT_SIGNED_SHORT:
            write_image_data(fp, nim, image->img.ss, mask, image->imgfmt == Image::FORMAT_FSL);
            break;
        case DT_UINT16:
            write_image_data(fp, nim, image->img.us, mask, image->imgfmt == Image::FORMAT_FSL);
            break;
        case DT_SIGNED_INT:
            write_image_data(fp, nim, image->img.si, mask, image->imgfmt == Image::FORMAT_FSL);
            break;
        case DT_FLOAT:
            write_image_data(fp, nim, image->img.fl, mask, image->imgfmt == Image::FORMAT_FSL);
            break;
        default:
            znzclose(fp);
            nifti_image_free(nim);
            return false;
    }
    nifti_image_free(nim);
	znzclose(fp);
    return true;
}

// ---------------------------------------------------------------------------
bool CanReadImage(const char* filename)
{
    if      (IsMetaImage (filename) != -1) return CanReadMetaImage (filename);
    else if (IsNiftiImage(filename) != -1) return CanReadNiftiImage(filename);
    else                                   return false;
}

// ---------------------------------------------------------------------------
Image* ReadImage(const char* filename, Image::Format fmt)
{
    return ReadImage(filename, DT_UNKNOWN, false, fmt);
}

// ---------------------------------------------------------------------------
Image* ReadImage(const char*   filename,
                 short         datatype,
                 bool          scale,
                 Image::Format fmt)
{
    if      (IsMetaImage (filename) != -1) return ReadMetaImage (filename, datatype, scale, fmt);
    else if (IsNiftiImage(filename) != -1) return ReadNiftiImage(filename, datatype, scale, fmt);
    else                                   return NULL;
}

// ---------------------------------------------------------------------------
Sequence ReadSequence(const char* filename, Image::Format fmt)
{
    return ReadSequence(filename, DT_UNKNOWN, false, fmt);
}

// ---------------------------------------------------------------------------
Sequence ReadSequence(const char*   filename,
                      short         datatype,
                      bool          scale,
                      Image::Format fmt)
{
    if (IsNiftiImage(filename) != -1) return ReadNiftiSequence(filename, datatype, scale, fmt);
    else                              return Sequence();
}

// ---------------------------------------------------------------------------
Image* ReadDeformationField(const char* filename, Image::Format fmt)
{
    if      (IsMetaImage (filename) != -1) return ReadMetaDeformationField (filename, fmt);
    else if (IsNiftiImage(filename) != -1) return ReadNiftiDeformationField(filename, fmt);
    else                                   return NULL;
}

// ---------------------------------------------------------------------------
bool WriteImage(const char* filename, const Image* image, Image::FileFormat fmt)
{
    return WriteImage(filename, image, NULL, fmt);
}

// ---------------------------------------------------------------------------
bool WriteImage(const char* filename, const Image* image, const Image* mask, Image::FileFormat fmt)
{
    if (filename == NULL || filename[0] == '\0') {
        BASIS_THROW(invalid_argument, "No filename specified for output image!");
    }
    if (fmt == Image::FILE_FORMAT_META) {
        return WriteMetaImage(filename, image, mask);
    } else if (fmt == Image::FILE_FORMAT_UNKNOWN) {
        set<string> mhdexts;
        mhdexts.insert(".mhd");
        mhdexts.insert(".def");
        set<string> rawexts;
        mhdexts.insert(".raw");
        mhdexts.insert(".raw.gz");
        if (os::path::hasext(filename, &mhdexts)) {
            return WriteMetaImage(filename, image, mask);
        } else if (os::path::hasext(filename, &rawexts)) {
            return WriteRawImage(filename, image, mask);
        } else {
            set<string> niftiexts;
            niftiexts.insert(".hdr");
            niftiexts.insert(".hdr.gz");
            niftiexts.insert(".img");
            niftiexts.insert(".img.gz");
            niftiexts.insert(".nii");
            niftiexts.insert(".nii.gz");
            niftiexts.insert(".nia");
            niftiexts.insert(".nia.gz");
            if (!os::path::hasext(filename, &niftiexts) &&
                    image->filefmt == Image::FILE_FORMAT_META) {
                return WriteMetaImage(filename, image, mask);
            }
        }
    }
    return WriteNiftiImage(filename, image, mask, FileFormatToNiftiType(fmt));
}

// ---------------------------------------------------------------------------
bool WriteSequence(const char* filename, const ConstSequence images, Image::FileFormat fmt)
{
    return WriteSequence(filename, images, NULL, fmt);
}

// ---------------------------------------------------------------------------
bool WriteSequence(const char* filename, const ConstSequence images, const Image* mask, Image::FileFormat fmt)
{
    if (filename == NULL || filename[0] == '\0') {
        BASIS_THROW(invalid_argument, "No filename specified for output image!");
    }
    if (images.empty()) {
        BASIS_THROW(invalid_argument, "No images to write!");
    }
    if (fmt == Image::FILE_FORMAT_META) {
        return false; // Meta image sequence output not supported
    } else if (fmt == Image::FILE_FORMAT_UNKNOWN) {
        set<string> mhdexts;
        mhdexts.insert(".mhd");
        mhdexts.insert(".def");
        set<string> rawexts;
        mhdexts.insert(".raw");
        mhdexts.insert(".raw.gz");
        if (os::path::hasext(filename, &mhdexts)) {
            return false; // Meta image sequence output not supported
        } else if (os::path::hasext(filename, &rawexts)) {
            return false; // Raw image sequence output not supported
        } else {
            set<string> niftiexts;
            niftiexts.insert(".hdr");
            niftiexts.insert(".hdr.gz");
            niftiexts.insert(".img");
            niftiexts.insert(".img.gz");
            niftiexts.insert(".nii");
            niftiexts.insert(".nii.gz");
            niftiexts.insert(".nia");
            niftiexts.insert(".nia.gz");
            if (!os::path::hasext(filename, &niftiexts) &&
                    images[0]->filefmt == Image::FILE_FORMAT_META) {
                return false; // Meta image sequence output not supported
            }
        }
    }
    return WriteNiftiSequence(filename, images, mask, FileFormatToNiftiType(fmt));
}

// ===========================================================================
// transformation - either affine or deformable
// ===========================================================================

// ---------------------------------------------------------------------------
bool ReadTransform(const char*       filename,
                   Image*&           deformation_field,
                   Image::Transform& transform,
                   Image::Format     fmt)
{
    // reset affine transformation matrix
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            transform.m[r][c] = 0;
        }
    }
    // try to read an affine 4x4 transformation matrix from the file
    if (ReadAffineTransform(filename, transform, fmt) == 1) {
        deformation_field = NULL;
    } else {
        // reset affine transformation matrix
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                transform.m[r][c] = 0;
            }
        }
        // try to read deformation field from the file
        deformation_field = ReadDeformationField(filename, fmt);
    }
    // if both read operations failed, the input file must be invalid
    return transform.m[3][3] != 0 || deformation_field != NULL;
}

// ===========================================================================
// NIfTI-1 or ANALYZE 7.5 image
// ===========================================================================

// ---------------------------------------------------------------------------
static bool can_read_image(const nifti_1_header& hdr)
{
    // major image dimensions, i.e., either 2D or 3D image
    if (hdr.dim[1] <= 1 || hdr.dim[2] <= 1 || hdr.dim[3] < 0) return false;
    // image sequences (4D image) are mainly processed sequentially
    if (hdr.dim[4] < 0) return false;
    // unused dimensions may not be set
    if (hdr.dim[6] > 1 || hdr.dim[7] > 1) return false;
    // datatype
    if (hdr.datatype != DT_UNSIGNED_CHAR &&
        hdr.datatype != DT_SIGNED_SHORT  &&
        hdr.datatype != DT_UINT16        &&
		hdr.datatype != DT_SIGNED_INT    &&
        hdr.datatype != DT_FLOAT) return false;
    // everything alright, can read the image
    return true;
}

// ---------------------------------------------------------------------------
int IsNiftiImage(const char* filename, nifti_1_header* hdr)
{
    nifti_set_debug_level(0); // discard debug message of nifti_is_complete_filename()
    int retval = is_nifti_file(filename);
    if (hdr && retval != -1) {
        nifti_1_header* nhdr = nifti_read_header(filename, NULL, 1);
        if (nhdr != NULL) {
            memcpy(hdr, nhdr, sizeof(nifti_1_header));
            free(nhdr);
        } else {
            retval = -1;
        }
    }
    nifti_set_debug_level(1); // reset debug level to default
                              // unfortunately, the nifticlib does not provide
                              // a way to retrieve the current debug level,
                              // so we cannot restore it to whatever it has
                              // been set on before calling this function
    return retval;
}

// ---------------------------------------------------------------------------
int IsNiftiDeformationField(const char* filename, nifti_1_header* hdr)
{
    nifti_set_debug_level(0); // discard debug message of nifti_is_complete_filename()
    int retval = is_nifti_file(filename);
    if (retval != -1) {
        nifti_1_header* nhdr = nifti_read_header(filename, NULL, 1);
        if (nhdr != NULL) {
            if ((nhdr->dim[0] != 3 && nhdr->dim[0] != 5) ||                       // 3D image
                    nhdr->dim[1] <= 1 || nhdr->dim[2] <= 1 || nhdr->dim[3] < 0 || // valid size
                    nhdr->dim[4] < 0 || nhdr->dim[4] > 1 ||                       // no time series
                    nhdr->dim[5] != 3) {                                          // 3 components
                retval = -1;
            }
            if (retval != -1 && hdr) memcpy(hdr, nhdr, sizeof(nifti_1_header));
            free(nhdr);
        } else {
            retval = -1;
        }
    }
    nifti_set_debug_level(1); // reset debug level to default
                              // unfortunately, the nifticlib does not provide
                              // a way to retrieve the current debug level,
                              // so we cannot restore it to whatever it has
                              // been set on before calling this function
    return retval;
}

// ---------------------------------------------------------------------------
static bool can_read_nifti_image(const char* filename, nifti_image** nim)
{
    nifti_image* _nim = nifti_image_read(filename, 0);
    if (!_nim) return false;
    nifti_1_header hdr = nifti_convert_nim2nhdr(_nim);
    if (can_read_image(hdr)) {
        if (nim) *nim = _nim;
        else nifti_image_free(_nim);
        return true;
    } else {
        nifti_image_free(_nim);
        return false;
    }
}

// ---------------------------------------------------------------------------
bool CanReadNiftiImage(const char* filename)
{
    return can_read_nifti_image(filename, NULL);
}

// ---------------------------------------------------------------------------
Image* ReadNiftiImage(const char* filename, Image::Format fmt)
{
    return ReadNiftiImage(filename, DT_UNKNOWN, false, fmt);
}

// ---------------------------------------------------------------------------
Image* ReadNiftiImage(const char*   filename,
                      short         datatype,
                      bool          scale,
                      Image::Format fmt)
{
  Sequence images = ReadNiftiSequence(filename, datatype, scale, fmt);
  if (images.size() == 1) return images[0];
  images.Clear();
  return NULL;
}

// ---------------------------------------------------------------------------
Sequence ReadNiftiSequence(const char* filename, Image::Format fmt)
{
    return ReadNiftiSequence(filename, DT_UNKNOWN, false, fmt);
}

// ---------------------------------------------------------------------------
Sequence ReadNiftiSequence(const char*   filename,
                           short         datatype,
                           bool          scale,
                           Image::Format fmt)
{
    Sequence     images;
    bool         ok  = true;
    nifti_image* nim = NULL;
    // check if input image is a ANALYZE or NIfTI image
    if (IsNiftiImage(filename) < 0) {
        DRAMMS_ERROR("Either failed to open image header of file " << filename << endl
                << "or this file is no image stored in the ANALYZE 7.5 or NIfTI-1 file format!");
        return images;
    }
    // read header and verify that we can read the image
    if (!can_read_nifti_image(filename, &nim)) {
        DRAMMS_ERROR("Cannot read image " << filename << "!" << endl
                << "Note that the image must be a 2D or 3D image of either datatype" << endl
                << "DT_UNSIGNED_CHAR, DT_SIGNED_SHORT, DT_UINT16, DT_SIGNED_INT, or DT_FLOAT.");
        return images;
    }
    // allocate image memory
    if (ok) {
        nifti_1_header hdr = nifti_convert_nim2nhdr(nim);
	    if (fmt == Image::FORMAT_DRAMMS) hdr = PermuteXY(hdr);
	    hdr.dim[0] = 5;
        hdr.dim[4] = 1;
        for (int i = 1; i <= 5; i++) {
            if (hdr.dim[i] <= 0) hdr.dim[i] = 1;
        }
        if (hdr.dim[5] > 1 && hdr.intent_code == NIFTI_INTENT_NONE) {
            // otherwise, the itk::NiftiImageIO reader will think it is a scalar image
            hdr.intent_code = NIFTI_INTENT_VECTOR;
        }
        int N = (nim->nt > 0 ? nim->nt : 1);
        for (int i = 0; i < N; i++) {
            Image *image = new Image(hdr.dim[1], hdr.dim[2], hdr.dim[3], hdr.datatype, hdr.dim[5], fmt);
            if (image) {
                image->hdr      = hdr;
                image->compress = nifti_is_gzfile(nim->iname);
                image->filefmt  = static_cast<Image::FileFormat>(nim->nifti_type + 1);
            } else {
                ok = false;
                images.Clear();
                break;
            }
            images.push_back(image);
        }
    }
    // open image data file
    znzFile fp = NULL;
    if (ok) {
        fp = nifti_image_open(filename, const_cast<char*>("rb"), &nim);
        ok = !znz_isnull(fp);
    }
    // determine offset of image data
    int offset = 0;
    if (ok) {
        if (nim->iname_offset < 0) {
            if (nifti_is_gzfile(nim->iname)) {
                // invalid offset for compressed data file of image
                ok = false;
            } else {
                int size = nifti_get_filesize(nim->iname);
                if (size <= 0) {
                    // empty data file
                    ok = false;
                } else {
                    int volsize = static_cast<int>(nifti_get_volsize(nim));
                    offset = ((size > volsize) ? (size - volsize) : 0);
                }
            }
        } else {
            offset = nim->iname_offset;
        }
    }
    // seek to begin of image data
    ok = ok && znzseek(fp, static_cast<long>(offset), SEEK_SET) != -1;
    // read and if necessary cast image sequence
    for (size_t i = 0; i < images.size(); i++) {
        // read image data
        if (ok) {
            switch (nim->datatype) {
                case DT_UNSIGNED_CHAR:
                    read_image_data(fp, nim, images[i]->img.uc, fmt == Image::FORMAT_FSL);
                    break;
                case DT_SIGNED_SHORT:
                    read_image_data(fp, nim, images[i]->img.ss, fmt == Image::FORMAT_FSL);
                    break;
                case DT_UINT16:
                    read_image_data(fp, nim, images[i]->img.us, fmt == Image::FORMAT_FSL);
                    break;
                case DT_SIGNED_INT:
                    read_image_data(fp, nim, images[i]->img.si, fmt == Image::FORMAT_FSL);
                    break;
                case DT_FLOAT:
                    read_image_data(fp, nim, images[i]->img.fl, fmt == Image::FORMAT_FSL);
                    break;
                default:
                    // should have been checked by CanReadNiftiImage() already
                    ASSERT(false, "Cannot read image of datatype: " << nim->datatype);
                    abort();
            }
        }
        // cast image if necessary and/or scale intensities
        if (ok) {
            if (datatype != DT_UNKNOWN && nim->datatype != datatype) {
                Image* tmp = CastImage(images[i], datatype, scale);
                delete images[i];
                images[i] = tmp;
            } else if (scale) {
                ScaleImage(images[i], images[i]);
            }
        }
    }
    // clean up
    if (!znz_isnull(fp)) znzclose(fp);
    if (nim) nifti_image_free(nim);
    if (!ok) images.Clear();
    return images;
}

// ---------------------------------------------------------------------------
Image* ReadNiftiDeformationField(const char* filename, Image::Format fmt)
{
    Image* def = NULL;
	Sequence input_images = ReadSequence(filename);
	int nt = input_images.size();
	if (nt>1) fmt=Image::FORMAT_FSL;
	
	if (IsNiftiDeformationField(filename)) {
        if (fmt == Image::FORMAT_FSL) { 
			// re-organize fsl deformation (4D sequence image, dim[4]=3) into 3D deformation vector file (dim[5]=3)
			if ( nt<2 || nt>3 )
				BASIS_THROW(invalid_argument, "The input deformation seems to be FSL deformation field, but the number of components should be either 2 or 3. Please double check.");
			
			cout << "Caution: Handling FSL deformation field for image warping or deformation operation. Please be aware that the results may or may not have orientation issues in the header." << endl;
			
			Image* input_imagex = input_images[1];
			Image* input_imagey = input_images[0];
			Image* input_imagez = NULL;
			int   nx = input_imagex->hdr.dim[1];
			int   ny = input_imagex->hdr.dim[2];
			int   nz = 1;
			float sx = input_imagex->hdr.pixdim[1];
			float sy = input_imagex->hdr.pixdim[2];
			float sz = 0; 
			if (sx<=0.0) sx=1.0;
			if (sy<=0.0) sy=1.0;
			if (nt==3) {
				input_imagez = input_images[2];
				nz = input_imagex->hdr.dim[3];
				if (nz<1) nz=1;
				sz = input_imagex->hdr.pixdim[3];
				if (sz<=0.0) sz=1.0;
				}
			//cout << "nx,ny,nz=" << ny << ", " << nx <<  ", " << nz << "; sx,sy,sz=" << sy << ", " << sx << ", " << sz << endl;
			def = new Image(nx, ny, nz, sx, sy, sz, DT_FLOAT, nt, Image::FORMAT_DRAMMS);
			
			vector<float> v(nt);
			for (int k=0; k<nz; k++)
			  for (int i=0; i<nx; i++)
				for (int j=0; j<ny; j++) {
					v[0] = input_imagex->img.fl[k][i][j]/sx;
					v[1] = input_imagey->img.fl[k][i][j]/sy;
					if (nt==3) 
						v[2] = input_imagez->img.fl[k][i][j]/sz;
					def->set(i, j, k, v);
					}
			
			if (input_imagex->hdr.qform_code == NIFTI_XFORM_UNKNOWN) 
				input_imagex->hdr.qform_code = NIFTI_XFORM_SCANNER_ANAT;
			def->CopyRegion(input_imagex);
			def->CopyTransform(input_imagex);
			def->CopyDataScaling(input_imagex);
			def->CopyMetaData(input_imagex);
			// int x_orient, y_orient, z_orient;
			// GetImageOrientation(input_imagex->hdr, x_orient, y_orient, z_orient);
			// def->SetOrientation(x_orient, y_orient, z_orient);
			}
		else
			def = ReadNiftiImage(filename, DT_FLOAT, false, fmt);
    }
	return def;
}

// ---------------------------------------------------------------------------
bool WriteNiftiImage(const char* filename, const Image* image, Image::NiftiType type)
{
    return WriteNiftiImage(filename, image, NULL, type);
}

// ---------------------------------------------------------------------------
bool WriteNiftiImage(const char* filename, const Image* image, const Image* mask, Image::NiftiType type)
{
  ConstSequence images;
  images.push_back(image);
  return WriteNiftiSequence(filename, images, mask, type);
}

// ---------------------------------------------------------------------------
bool WriteNiftiSequence(const char* filename, const ConstSequence images, Image::NiftiType type)
{
    return WriteNiftiSequence(filename, images, NULL, type);
}

// ---------------------------------------------------------------------------
bool WriteNiftiSequence(const char* filename, const ConstSequence images, const Image* mask, Image::NiftiType type)
{
    if (filename == NULL || filename[0] == '\0') {
        BASIS_THROW(invalid_argument, "No filename specified for output image!");
    }
    if (mask && mask->hdr.datatype != DT_UNSIGNED_CHAR) {
        BASIS_THROW(invalid_argument, "Mask image must be of datatype unsigned char!");
    }
    // -----------------------------------------------------------------------
    // If image data covers only parts of the image domain, resize the image,
    // i.e., "resample" it such that the image data is defined on the grid
    // specified by the image header. See also ResizeImage().
    ConstSequence img;        // (resized) images
    const Image  *msk = NULL; // (resized) mask
    for (size_t i = 0; i < images.size(); i++) {
        const Image *image = images[i];
        if (image->region.ox != 0 ||
                image->region.oy != 0 ||
                image->region.oz != 0 ||
                image->region.nx != image->hdr.dim[1] ||
                image->region.ny != image->hdr.dim[2] ||
                image->region.nz != image->hdr.dim[3]) {
            image = ResizeImage(image);
            if (image == NULL) {
                for (size_t j = 0; j < i; j++) {
                    if (img[j] != images[j]) delete img[j];
                }
                DRAMMS_ERROR("Failed to resize image to its original size before writing it!");
                return false;
            }
        }
        if (i > 0 && !image->HasSameSizeAs(img[0])) {
            for (size_t j = 0; j < i; j++) {
                if (img[j] != images[j]) delete img[j];
            }
            BASIS_THROW(invalid_argument, "Not all images of the sequence have the same size!");
        }
        img.push_back(image);
    }
    if (mask) {
        if (mask->region.ox != 0 ||
                mask->region.oy != 0 ||
                mask->region.oz != 0 ||
                mask->region.nx != mask->hdr.dim[1] ||
                mask->region.ny != mask->hdr.dim[2] ||
                mask->region.nz != mask->hdr.dim[3]) {
            msk = ResizeImage(mask);
            if (msk == NULL) {
                for (size_t i = 0; i < images.size(); i++) {
                    if (img[i] != images[i]) delete img[i];
                }
                DRAMMS_ERROR("Failed to resize mask to its original size before writing the image(s)!");
                return false;
            }
        } else {
            msk = mask;
        }
        if (mask && !mask->HasSameSizeAs(img[0])) {
            for (size_t i = 0; i < images.size(); i++) {
                if (img[i] != images[i]) delete img[i];
            }
            if (msk != mask) delete msk;
            BASIS_THROW(invalid_argument, "Mask image has not the same size as the image(s) to write!");
        }
    }
    // -----------------------------------------------------------------------
    // first frame of sequence used as reference for header attributes
    const Image * const image = img[0];
    // -----------------------------------------------------------------------
    // determine the correct/desired output nifti_type and set explicit
    // filename extension if required to force compression
    string fname = filename;
    set<string> analyzeexts; // note: can be also NIfTI file extension!
    analyzeexts.insert(".hdr");
    analyzeexts.insert(".hdr.gz");
    analyzeexts.insert(".img");
    analyzeexts.insert(".img.gz");
    set<string> niftiexts;
    niftiexts.insert(".nii");
    niftiexts.insert(".nii.gz");
    set<string> niaexts;
    niaexts.insert(".nia");
    niaexts.insert(".nia.gz");
    Image::NiftiType nifti_type = (type != Image::NIFTI_TYPE_UNKNOWN)
                                       ? type
                                       : FileFormatToNiftiType(image->filefmt);
    if (os::path::hasext(filename, &niftiexts)) {
        nifti_type = Image::NIFTI_TYPE_NIFTI_1;
    } else if (os::path::hasext(filename, &niaexts)) {
        nifti_type = Image::NIFTI_TYPE_NIFTI_A;
    } else if (os::path::hasext(filename, &analyzeexts)) {
        if (nifti_type != Image::NIFTI_TYPE_ANALYZE) {
            nifti_type = Image::NIFTI_TYPE_NIFTI_2;
        }
    } else {
        if (nifti_type == Image::NIFTI_TYPE_ANALYZE || strncmp(image->hdr.magic, "ni1", 4) == 0) {
            fname += ".hdr";
            if (image->compress) fname += ".gz";
            if (nifti_type == Image::NIFTI_TYPE_UNKNOWN) nifti_type = Image::NIFTI_TYPE_NIFTI_2;
        } else {
            fname += ".nii";
            if (image->compress) fname += ".gz";
            if (nifti_type == Image::NIFTI_TYPE_UNKNOWN) nifti_type = Image::NIFTI_TYPE_NIFTI_1;
        }
    }
    // set orientation information for NIfTI
    nifti_1_header nhdr = image->hdr;
    if      (nifti_type == Image::NIFTI_TYPE_NIFTI_1) strncpy(nhdr.magic, "n+1",    4);
    else if (nifti_type == Image::NIFTI_TYPE_NIFTI_2) strncpy(nhdr.magic, "ni1",    4);
    else                                              strncpy(nhdr.magic, "\0\0\0", 4);
    if (nhdr.qform_code == NIFTI_XFORM_UNKNOWN) {
        Image::Transform qto_xyz = GetQFormTransform(image->hdr, image->imgfmt);
        nifti_mat44_to_quatern(qto_xyz, &nhdr.quatern_b, &nhdr.quatern_c, &nhdr.quatern_d,
                                        &nhdr.qoffset_x, &nhdr.qoffset_y, &nhdr.qoffset_z,
                                        NULL, NULL, NULL, &nhdr.pixdim[0]);
        nhdr.qform_code = NIFTI_XFORM_SCANNER_ANAT;
    }
    // swap data to radiological order in case of ANALYZE
    if (nifti_type == Image::NIFTI_TYPE_ANALYZE && GetLeftRightOrder(image->hdr, image->imgfmt) == Image::LR_NEUROLOGICAL) {
        for (size_t i = 0; i < images.size(); i++) {
            // revert direction of y axis to turn RAI into RPI...
            Image *oriented_image = (img[i] == images[i]) ? CopyImage(image) : const_cast<Image *>(img[i]);
            if (oriented_image == NULL) {
                for (size_t j = 0; j < images.size(); j++) {
                    if (img[j] != images[j]) delete img[j];
                }
                if (msk != mask) delete msk;
                DRAMMS_ERROR("Failed to allocate temporary memory to swap image data to radiological order!");
                return false;
            }
            oriented_image->SetFormat(Image::FORMAT_ITK);
            int x_orient, y_orient, z_orient;
            GetImageOrientation(oriented_image->hdr, x_orient, y_orient, z_orient, oriented_image->imgfmt);
            if (oriented_image->imgfmt == Image::FORMAT_DRAMMS) {
                x_orient += (x_orient % 2) == 0 ? -1 : +1; // see NIFTI_A2P and NIFTI_P2A codes
            } else {
                y_orient += (y_orient % 2) == 0 ? -1 : +1;
            }
            OrientImage(oriented_image, x_orient, y_orient, z_orient);
            img[i] = oriented_image;
        }
    }
    // adjust number of dimensions, e.g.,
    // 3 in case of 3D images instead of 5 with 4th and 5th dimension of size 1 only
    nhdr.dim[0] = 7;
    nhdr.dim[4] = static_cast<int>(images.size());
    while (nhdr.dim[0] > 0 && nhdr.dim[nhdr.dim[0]] <= 1) {
        nhdr.dim   [nhdr.dim[0]] = 0;
        //nhdr.pixdim[nhdr.dim[0]] = 0.0f;
        nhdr.dim[0]--;
    }
    for (int i = 1; i < nhdr.dim[0]+2; i++) {
        if (nhdr.dim[i]    <= 0)    nhdr.dim   [i] = 1;
        if (nhdr.pixdim[i] == 0.0f) nhdr.pixdim[i] = 1.0f;
	if (nhdr.pixdim[i] < 0.0f)  nhdr.pixdim[i] = -nhdr.pixdim[i];
    }
    // set intent to NIFTI_INTENT_VECTOR for deformation fields, otherwise the ITK
    // itk::NiftiImageIO class will read the image as scalar image
    if (nhdr.dim[5] > 1 && nhdr.intent_code == NIFTI_INTENT_NONE) {
        nhdr.intent_code = NIFTI_INTENT_VECTOR;
    } else if (nhdr.dim[5] <= 1 && nhdr.intent_code == NIFTI_INTENT_VECTOR) {
        nhdr.intent_code = NIFTI_INTENT_NONE;
    }
    // set default description if none is set
    if (nhdr.descrip[0] == '\0') {
        string descrip = string("DRAMMS ") + RELEASE;
        strncpy(nhdr.descrip, descrip.c_str(), min(static_cast<unsigned int>(descrip.size()), 78U));
        nhdr.descrip[79] = '\0';
    }
    // convert NIfTI-1 header to nifti_image structure
    if (image->imgfmt == Image::FORMAT_DRAMMS) nhdr = PermuteXY(nhdr);
    nifti_image *nim = nifti_convert_nhdr2nim(nhdr, fname.c_str());
    if (!nim) {
        for (size_t i = 0; i < images.size(); i++) {
            if (img[i] != images[i]) delete img[i];
        }
        if (msk != mask) delete msk;
        return false;
    }
    nim->nifti_type = static_cast<int>(nifti_type) - 1;
    // make output directory if it does not exist yet
    if (!os::makedirs(os::path::dirname(filename))) {
        for (size_t i = 0; i < images.size(); i++) {
            if (img[i] != images[i]) delete img[i];
        }
        if (msk != mask) delete msk;
        return false;
    }
    // write image header
    znzFile fp = nifti_image_write_hdr_img(nim, 2, "wb");
    if (znz_isnull(fp)) {
        nifti_image_free(nim);
        for (size_t i = 0; i < images.size(); i++) {
            if (img[i] != images[i]) delete img[i];
        }
        if (msk != mask) delete msk;
        return false;
    }
    // write image data
    for (size_t i = 0; i < images.size(); i++) {
        switch (image->hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                write_image_data(fp, nim, img[i]->img.uc, mask, image->imgfmt == Image::FORMAT_FSL);
                break;
            case DT_SIGNED_SHORT:
                write_image_data(fp, nim, img[i]->img.ss, mask, image->imgfmt == Image::FORMAT_FSL);
                break;
            case DT_UINT16:
                write_image_data(fp, nim, img[i]->img.us, mask, image->imgfmt == Image::FORMAT_FSL);
                break;
            case DT_SIGNED_INT:
                write_image_data(fp, nim, img[i]->img.si, mask, image->imgfmt == Image::FORMAT_FSL);
                break;
            case DT_FLOAT:
                write_image_data(fp, nim, img[i]->img.fl, mask, image->imgfmt == Image::FORMAT_FSL);
                break;
            default:
                znzclose(fp);
                nifti_image_free(nim);
                for (size_t i = 0; i < images.size(); i++) {
                    if (img[i] != images[i]) delete img[i];
                }
                if (msk != mask) delete msk;
                return false;
        }
    }
	znzclose(fp);
	nifti_image_free(nim);
    for (size_t i = 0; i < images.size(); i++) {
        if (img[i] != images[i]) delete img[i];
    }
    if (msk != mask) delete msk;
    return true;
}

// ===========================================================================
// Meta image
// ===========================================================================

/**
 * @brief Enumeration of supported Meta image element types.
 *
 * Note that these enumeration values are set to the corresponding NIfTI-1
 * DT_* datatype values to simplify the conversion from one to the other.
 */
enum MetaType {
    MET_NONE   = DT_UNKNOWN,
    MET_UCHAR  = DT_UNSIGNED_CHAR,
    MET_SHORT  = DT_INT16,
    MET_USHORT = DT_UINT16,
    MET_INT    = DT_INT32,
    MET_UINT   = DT_UINT32,
    MET_FLOAT  = DT_FLOAT,
    MET_OTHER  = DT_ALL ///< Unsupported Meta image type.
}; // enum MetaType

/**
 * @brief Contains information read from Meta image header file.
 */
struct MetaHeader
{
    std::string datafile;    ///< Name of image data file.
    int         header_size; ///< Number of bytes to skip at top of data file.
    int         ndims;       ///< Number of image dimensions.
    int         nchannels;   ///< Number of channels.
    int*        size;        ///< Size of each image dimension.
    float*      spacing;     ///< Spacing between elements, i.e., voxels.
    float*      elemsize;    ///< Size of each element, i.e., voxel.
    float*      offset;      ///< Image offset of size ndims;
    float**     transform;   ///< Transformation matrix of size ndims x ndims.
    MetaType    type;        ///< Type of element data.
    bool        compressed;  ///< Whether image data is compressed.
    bool        binary;      ///< Whether the data is stored in binary format.
    bool        msb;         ///< Whether the binary data is stored in MSB order.
    short       bitpix;      ///< Number of bits per element.

    /**
     * @brief Default constructor.
     */
    MetaHeader();

    /**
     * @brief Copy constructor.
     *
     * @param [in] other Other object.
     */
    MetaHeader(const MetaHeader& other);

    /**
     * @brief Destructor.
     */
    ~MetaHeader();

    /**
     * @brief Set number of dimensions.
     *
     * This method (re-)allocates all fields which depend on the number of dimensions.
     */
    void SetNDims(int n);

    /**
     * @brief Reset header to initial/invalid state.
     */
    void Clear();

    /**
     * @brief Assignment operator.
     *
     * @param [in] rhs Right-hand side.
     *
     * @returns Reference to this instance after the assignment.
     */
    MetaHeader& operator=(const MetaHeader& rhs);

}; // struct MetaHeader

// ---------------------------------------------------------------------------
MetaHeader::MetaHeader()
:
    size     (NULL),
    spacing  (NULL),
    elemsize (NULL),
    offset   (NULL),
    transform(NULL)
{
    Clear();
}

// ---------------------------------------------------------------------------
MetaHeader::MetaHeader(const MetaHeader& other)
:
    size     (NULL),
    spacing  (NULL),
    elemsize (NULL),
    offset   (NULL),
    transform(NULL)
{
    *this = other;
}

// ---------------------------------------------------------------------------
MetaHeader::~MetaHeader()
{
    Clear();
}

// ---------------------------------------------------------------------------
void MetaHeader::SetNDims(int n)
{
    if (size) {
        delete [] size;
        size = NULL;
    }
    if (spacing) {
        delete [] spacing;
        spacing = NULL;
    }
    if (elemsize) {
        delete [] elemsize;
        elemsize = NULL;
    }
    if (offset) {
        delete [] offset;
        offset = NULL;
    }
    if (transform) {
        for (int i = 0; i < ndims; i++) delete [] transform[i];
        delete [] transform;
        transform = NULL;
    }
    ndims = n > 0 ? n : 0;
    if (ndims > 0) {
        size      = new int  [ndims];
        spacing   = new float[ndims];
        elemsize  = new float[ndims];
        offset    = new float[ndims];
        transform = new float*[ndims];
        for (int i = 0; i < ndims; i++) {
            size     [i] = 0;
            spacing  [i] = 0;
            elemsize [i] = 0;
            offset   [i] = 0;
            transform[i] = new float[ndims];
            for (int c = 0; c < ndims; c++) {
                transform[i][c] = 0;
            }
        }
    }
}

// ---------------------------------------------------------------------------
MetaHeader& MetaHeader::operator=(const MetaHeader& rhs)
{
    header_size  = rhs.header_size;
    type         = rhs.type;
    compressed   = rhs.compressed;
    binary       = rhs.binary;
    msb          = rhs.msb;
    bitpix       = rhs.bitpix;
    nchannels    = rhs.nchannels;
    datafile     = rhs.datafile;
    SetNDims(rhs.ndims);
    for (int i = 0; i < ndims; i++) {
        size    [i] = rhs.size    [i];
        spacing [i] = rhs.spacing [i];
        elemsize[i] = rhs.elemsize[i];
        offset  [i] = rhs.offset  [i];
        for (int c = 0; c < ndims; c++) {
            transform[i][c] = rhs.transform[i][c];
        }
    }
    return *this;
}

// ---------------------------------------------------------------------------
void MetaHeader::Clear()
{
    header_size  = 0;
    type         = MET_NONE;
    compressed   = false;
    binary       = true;
    msb          = false;
    bitpix       = 0;
    nchannels    = 0;
    datafile.clear();
    SetNDims(0);
}

// ---------------------------------------------------------------------------
static MetaType string_to_meta_type(const string& str)
{
    if      (str == "MET_NONE")                            return MET_NONE;
    else if (str == "MET_UCHAR")                           return MET_UCHAR;
    else if (str == "MET_SHORT")                           return MET_SHORT;
    else if (str == "MET_USHORT")                          return MET_USHORT;
    else if (str == "MET_INT")                             return MET_INT;
    else if (str == "MET_UINT")                            return MET_UINT;
    else if (str == "MET_FLOAT")                           return MET_FLOAT;
    else if (str.size() > 4 && str.substr(0, 4) == "MET_") return MET_OTHER;
    else                                                   return MET_NONE;
}

// ---------------------------------------------------------------------------
static const char* meta_type_to_string(MetaType type)
{
    switch (type) {
        case MET_NONE:   return "MET_NONE";
        case MET_UCHAR:  return "MET_UCHAR";
        case MET_SHORT:  return "MET_SHORT";
        case MET_USHORT: return "MET_USHORT";
        case MET_INT:    return "MET_INT";
        case MET_UINT:   return "MET_UINT";
        case MET_FLOAT:  return "MET_FLOAT";
        default:         return "MET_OTHER";
    };
}

// ---------------------------------------------------------------------------
static short meta_to_nifti_type(MetaType type)
{
    return static_cast<short>(type);
}

// ---------------------------------------------------------------------------
static MetaType nifti_to_meta_type(short datatype)
{
    return static_cast<MetaType>(datatype);
}

// ---------------------------------------------------------------------------
static short meta_type_bitpix(const MetaType& type)
{
    switch (type)
    {
        case MET_UCHAR:  return 8;
        case MET_SHORT:  return 16;
        case MET_USHORT: return 16;
        case MET_INT:    return 32;
        case MET_UINT:   return 32;
        case MET_FLOAT:  return 32;
        default:         return 0;
    };
}

// ---------------------------------------------------------------------------
static short meta_type_swapsize(MetaType type)
{
    switch (type)
    {
        case MET_UCHAR:  return 1;
        case MET_SHORT:  return 2;
        case MET_USHORT: return 2;
        case MET_INT:    return 4;
        case MET_UINT:   return 4;
        case MET_FLOAT:  return 4;
        default:         return 0;
    };
}

// ---------------------------------------------------------------------------
static bool read_meta_header(const char* filename, MetaHeader& mhd)
{
    bool   ok = true;
    string key;
    string sval;
    float  fval;
    int    ival;
    char   cval;

    // reset header
    mhd.Clear();
    // open header file
    string mhd_file = filename;
    string mhd_ext  = os::path::splitext(mhd_file)[1];
    if      (mhd_ext == ".def") mhd_file.replace(mhd_file.size() - 3, 3, "mhd");
    else if (mhd_ext != ".mhd") mhd_file += ".mhd";
    ifstream ifs(mhd_file.c_str());
    if (!ifs.is_open()) return false;
    // parse header entries, line-by-line
    while (ifs.good()) {
        // get key of next header entry
        ifs >> key;
        if (!ifs.good()) break;
        // require presence of equal sign (=) as key/value separator
        while (ifs.good()) { if ((cval = ifs.get()) != ' ') break; }
        if (!ifs.good() || cval != '=') { ok = false; break; }
        // ObjectType
        if (key == "ObjectType") {
            ifs >> sval;
            if (sval != "Image") { ok = false; break; }
        // NDims
        } else if (key == "NDims") {
            if (mhd.ndims != 0) { ok = false; break; }
            ifs >> ival;
            if (!ifs.good() || ival <= 0) { ok = false; break; }
            mhd.SetNDims(ival);
        // BinaryData
        } else if (key == "BinaryData") {
            ifs >> sval;
            mhd.binary = (sval == "True" || sval == "true" || sval == "TRUE");
        // BinaryDataByteOrderMSB
        } else if (key == "BinaryDataByteOrderMSB") {
            ifs >> sval;
            mhd.msb = (sval == "True" || sval == "true" || sval == "TRUE");
        // CompressedData
        } else if (key == "CompressedData") {
            ifs >> sval;
            mhd.compressed = (sval == "True" || sval == "true" || sval == "TRUE");
        // AnatomicalOrientation
        } else if (key == "AnatomicalOrientation") {
            // TODO What to do about this three letter orientation code ?
            //      How does it relate to the TransformMatrix ?
        // TransformMatrix
        } else if (key == "TransformMatrix" || key == "Orientation" | key == "Rotation") {
            if (mhd.transform == NULL) { ok = false; break; }
            for (int r = 0; r < mhd.ndims; r++) {
                for (int c = 0; c < mhd.ndims; c++) {
                    ifs >> fval;
                    mhd.transform[r][c] = fval;
                }
            }
        // Offset
        } else if (key == "Offset" || key == "Position" || key == "Origin") {
            if (mhd.offset == NULL) { ok = false; break; }
            for (int i = 0; i < mhd.ndims; i++) {
                ifs >> fval;
                mhd.offset[i] = fval;
            }
        // ElementSpacing
        } else if (key == "ElementSpacing") {
            if (mhd.spacing == NULL) { ok = false; break; }
            for (int i = 0; i < mhd.ndims; i++) {
                ifs >> fval;
                mhd.spacing[i] = fval;
            }
        // ElementSize
        } else if (key == "ElementSize") {
            if (mhd.elemsize == NULL) { ok = false; break; }
            for (int i = 0; i < mhd.ndims; i++) {
                ifs >> fval;
                mhd.elemsize[i] = fval;
            }
        // DimSize
        } else if (key == "DimSize") {
            if (mhd.size == NULL) { ok = false; break; }
            for (int i = 0; i < mhd.ndims; i++) {
                ifs >> ival;
                mhd.size[i] = ival;
            }
        // ElementNumberOfChannels
        } else if (key == "ElementNumberOfChannels") {
            ifs >> mhd.nchannels;
        // ElementType
        } else if (key == "ElementType") {
            ifs >> sval;
            mhd.type   = string_to_meta_type(sval);
            mhd.bitpix = meta_type_bitpix(mhd.type);
        // ElementDataFile
        } else if (key == "ElementDataFile") {
            ifs >> mhd.datafile;
            if (!mhd.datafile.empty()) {
                mhd.datafile = os::path::join(
                                    os::path::abspath(os::path::dirname(mhd_file)),
                                    mhd.datafile);
            }
        }
        // ignore all other keys...
        // require that each entry is on separate line
        while (ifs.good()) { if (ifs.get() == '\n') break; }
    }
    ok = ok && ifs.eof();
    ok = ok && mhd.ndims > 0;
    if (ok) {
        if (mhd.spacing[0] <= 0) {
            if (mhd.elemsize[0] > 0) {
                for (int i = 0; i < mhd.ndims; i++) mhd.spacing[i] = mhd.elemsize[i];
            } else {
                for (int i = 0; i < mhd.ndims; i++) mhd.spacing[i] = 1.0f;
            }
        }
    }
    if (!ok) mhd.Clear();
    return ok;
}

// ---------------------------------------------------------------------------
static nifti_1_header* convert_mhd2nhdr(const MetaHeader& mhd)
{
    // initialize invalid header
    nifti_1_header* hdr = NULL;
    // bail out if image is not supported by this function
    if (mhd.ndims < 2 || mhd.ndims > 3) return NULL;
    // create new NIfTI-1 header
    {
        int dim[8] = {0, 0, 0, 0, 0, 0, 0, 0};
        if (mhd.nchannels > 1) {
            dim[0] = 5;
            dim[1] = mhd.size[0];
            dim[2] = mhd.size[1];
            dim[3] = mhd.ndims > 2 ? mhd.size[2] : 1;
            dim[4] = 1;
            dim[5] = mhd.nchannels;
        } else {
            dim[0] = mhd.ndims;
            dim[1] = mhd.size[0];
            dim[2] = mhd.size[1];
            dim[3] = mhd.ndims > 2 ? mhd.size[2] : 1;
        }
        short datatype = meta_to_nifti_type(mhd.type);
        if (datatype == DT_UNKNOWN) return NULL;
        hdr = nifti_make_new_header(dim, datatype);
        if (hdr == NULL) return NULL;
    }
    // Q-form matrix
    bool transform_set = false;
    for (int r = 0; r < mhd.ndims; r++) {
        for (int c = 0; c < mhd.ndims; c++) {
            if (mhd.transform[r][c] != 0) transform_set = true;
        }
    }
    if (transform_set) {
        mat44 qmat;
        for (int r = 0; r < mhd.ndims; r++) {
            for (int c = 0; c < mhd.ndims; c++) {
                qmat.m[r][c] = mhd.transform[r][c];
                if (r == c) qmat.m[r][c] *= mhd.spacing[r];
            }
        }
        for (int r = mhd.ndims; r < 4; r++) {
            for (int c = mhd.ndims; c < 4; c++) {
                qmat.m[r][c] = 0.0f;
            }
        }
        for (int r = 0; r < mhd.ndims; r++) {
            qmat.m[r][3] = mhd.offset[r];
        }
        // convert Q-form matrix to quaternion and offset representation
        hdr->qform_code = NIFTI_XFORM_SCANNER_ANAT;
        nifti_mat44_to_quatern(qmat, &hdr->quatern_b, &hdr->quatern_c, &hdr->quatern_d,
                                     &hdr->qoffset_x, &hdr->qoffset_y, &hdr->qoffset_z,
                                     &hdr->pixdim[1], &hdr->pixdim[2], &hdr->pixdim[3],
                                     &hdr->pixdim[0]);
    }
    return hdr;
}

// ---------------------------------------------------------------------------
static MetaHeader convert_nhdr2mhd(const nifti_1_header& hdr)
{
    MetaHeader mhd;
    mhd.SetNDims(hdr.dim[3] > 1 ? 3 : 2);
    mhd.nchannels = hdr.dim[5] > 1 ? hdr.dim[5] : 0;
    for (int i = 0; i < mhd.ndims; i++) {
        mhd.size   [i] = hdr.dim   [i + 1];
        mhd.spacing[i] = hdr.pixdim[i + 1];
    }
    mhd.offset[0] = hdr.qoffset_x;
    mhd.offset[1] = hdr.qoffset_y;
    if (mhd.ndims > 2) mhd.offset[2] = hdr.qoffset_z;
    if (hdr.qform_code != NIFTI_XFORM_UNKNOWN) {
        mat44 qmat = nifti_quatern_to_mat44(hdr.quatern_b,
                                            hdr.quatern_c,
                                            hdr.quatern_d,
                                            hdr.qoffset_x,
                                            hdr.qoffset_y,
                                            hdr.qoffset_z,
                                            hdr.pixdim[1],
                                            hdr.pixdim[2],
                                            hdr.pixdim[3],
                                            hdr.pixdim[0] >= 0.0 ? 1.0 : -1.0);
        for (int r = 0; r < mhd.ndims; r++) {
            for (int c = 0; c < mhd.ndims; c++) {
                mhd.transform[r][c] = qmat.m[r][c];
                if (r == c) mhd.transform[r][c] /= mhd.spacing[r];
            }
        }
    }
    mhd.type   = nifti_to_meta_type(hdr.datatype);
    mhd.msb    = (nifti_short_order() == MSB_FIRST);
    mhd.bitpix = hdr.bitpix;
    return mhd;
}

// ---------------------------------------------------------------------------
static size_t meta_read_buffer(znzFile fp, void* dataptr, size_t ntot, short swapsize, const MetaHeader& mhd)
{
    if (dataptr == NULL) return -1;
    size_t ii = znzread(dataptr, 1, ntot, fp);
    if (ii < ntot) return -1 ;
    if  (swapsize > 1 && (nifti_short_order() == LSB_FIRST) && mhd.msb) {
        nifti_swap_Nbytes(ntot / swapsize, swapsize, dataptr);
    }
    return ii;
}

// ---------------------------------------------------------------------------
template <typename T>
static bool read_image_data(znzFile fp, const MetaHeader& mhd, T*** p, bool merge)
{
    const int   N      = (mhd.nchannels > 0) ? mhd.nchannels : 1; // components per voxel
    const int   size_x = mhd.size[0];
    const int   size_y = mhd.size[1];
    const int   size_z = mhd.ndims > 2 ? mhd.size[2] : 1;
    const short swapsize = meta_type_swapsize(mhd.type);
    if (merge) {
        const size_t nbytes = size_x * mhd.bitpix / 8; // bytes per row
        T* row = new T[nbytes / sizeof(T)];
        for (int n = 0; n < N; n++) {
            for (int k = 0; k < size_z; k++) {
                for (int j = 0; j < size_y; j++) {
                    if (meta_read_buffer(fp, row, nbytes, swapsize, mhd) != nbytes) {
                        return false;
                    }
                    for (int i = 0; i < size_x; i++) {
                        p[k][j][i * N + n] = row[i];
                    }
                }
            }
        }
        delete [] row;
    } else {
        const size_t nbytes = N * size_x * mhd.bitpix / 8; // bytes per row
        for (int k = 0; k < size_z; k++) {
            for (int j = 0; j < size_y; j++) {
                if (meta_read_buffer(fp, p[k][j], nbytes, swapsize, mhd) != nbytes) {
                    return false;
                }
            }
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
static int is_meta_image(const char* filename, MetaHeader* mhd = NULL)
{
    MetaHeader _mhd;
    if (!read_meta_header(filename, _mhd)) return -1;
    if (mhd) *mhd = _mhd;
    return 1;
}

// ---------------------------------------------------------------------------
int IsMetaImage(const char* filename, nifti_1_header* hdr)
{
    if (hdr) {
        MetaHeader mhd;
        int retval = is_meta_image(filename, &mhd);
        if (retval != -1 && hdr) {
            nifti_1_header* _hdr = convert_mhd2nhdr(mhd);
            if (_hdr == NULL) return -1;
            *hdr = *_hdr;
            free(_hdr);
        }
        return retval;
    } else {
        return is_meta_image(filename);
    }
}

// ---------------------------------------------------------------------------
int IsMetaDeformationField(const char* filename, nifti_1_header* hdr)
{
    MetaHeader mhd;
    int retval = is_meta_image(filename, &mhd);
    if (retval != -1) {
        if (mhd.ndims != 3 || mhd.nchannels != 3 ||                        // 3D vector field
                mhd.size[0] <= 1 || mhd.size[1] <= 1 || mhd.size[2] < 0) { // valid size
            retval = -1;
        }
    }
    if (retval != -1 && hdr) {
        nifti_1_header* _hdr = convert_mhd2nhdr(mhd);
        if (_hdr == NULL) return -1;
        *hdr = *_hdr;
        free(_hdr);
    }
    return retval;
}

// ---------------------------------------------------------------------------
bool CanReadMetaImage(const char* filename)
{
    MetaHeader mhd;
    if (!read_meta_header(filename, mhd) || !mhd.binary || mhd.type == MET_OTHER) return false;
    nifti_1_header* hdr = convert_mhd2nhdr(mhd);
    if (hdr == NULL) return false;
    bool ok = can_read_image(*hdr);
    free(hdr);
    return ok;
}

// ---------------------------------------------------------------------------
Image* ReadMetaImage(const char* filename, Image::Format fmt)
{
    return ReadMetaImage(filename, DT_UNKNOWN, false, fmt);
}

// ---------------------------------------------------------------------------
Image* ReadMetaImage(const char*   filename,
                     short         datatype,
                     bool          scale,
                     Image::Format fmt)
{
    bool       ok  = true;
    MetaHeader mhd;
    // check if input image is a Meta image
    if (is_meta_image(filename, &mhd) < 0) {
        DRAMMS_ERROR("Either failed to open image header of file " << filename << endl
                << "or this file is no image stored in the Meta image file format!");
        return NULL;
    }
    // verify that we can read the image
    nifti_1_header* hdr = convert_mhd2nhdr(mhd);
    if (hdr == NULL || !mhd.binary || !can_read_image(*hdr)) {
        DRAMMS_ERROR("Cannot read image " << filename << "!" << endl
                << "Note that the image must be a 2D or 3D image of either datatype" << endl
                << "MET_UCHAR, MET_SHORT, MET_USHORT, MET_INT, or MET_FLOAT.");
        return NULL;
    }
    // allocate image memory
    Image* image = NULL;
    if (ok) {
        if (fmt == Image::FORMAT_DRAMMS) *hdr = PermuteXY(*hdr);
        hdr->dim[0] = 5;
        for (int i = 1; i <= 5; i++) {
            if (hdr->dim[i] <= 0) hdr->dim[i] = 1;
        }
        image = new Image(hdr->dim[1], hdr->dim[2], hdr->dim[3], hdr->datatype, hdr->dim[5], fmt);
        if (image) {
            image->hdr      = *hdr;
            image->filefmt  = Image::FILE_FORMAT_META;
            image->imgfmt   = fmt;
            image->compress = mhd.compressed;
        } else {
            ok = false;
        }
    }
    // open image data file
    znzFile fp = NULL;
    if (ok) {
        fp = znzopen(mhd.datafile.c_str(), const_cast<char*>("rb"), mhd.compressed);
        if (znz_isnull(fp)) {
            DRAMMS_ERROR("Cannot open image data file " << mhd.datafile << " for reading!");
            ok = false;
        }
    }
    // read image data
    if (ok) {
        switch (hdr->datatype) {
            case DT_UNSIGNED_CHAR:
                read_image_data(fp, mhd, image->img.uc, fmt == Image::FORMAT_FSL);
                break;
            case DT_SIGNED_SHORT:
                read_image_data(fp, mhd, image->img.ss, fmt == Image::FORMAT_FSL);
                break;
            case DT_UINT16:
                read_image_data(fp, mhd, image->img.us, fmt == Image::FORMAT_FSL);
                break;
            case DT_SIGNED_INT:
                read_image_data(fp, mhd, image->img.si, fmt == Image::FORMAT_FSL);
                break;
            case DT_FLOAT:
                read_image_data(fp, mhd, image->img.fl, fmt == Image::FORMAT_FSL);
                break;
            default:
                // should have been checked by CanReadMetaImage() already
                ASSERT(false, "Cannot read image of datatype: " << hdr->datatype);
                abort();
        }
    }
    // cast image if necessary and/or scale intensities
    if (ok) {
        if (datatype != DT_UNKNOWN && hdr->datatype != datatype) {
            Image* tmp = CastImage(image, datatype, scale);
            delete image;
            image = tmp;
        } else if (scale) {
            ScaleImage(image, image);
        }
    }
    // clean up
    free(hdr);
    if (!znz_isnull(fp)) znzclose(fp);
    if (!ok) {
        delete image;
        image = NULL;
    }
    return image;
}

// ---------------------------------------------------------------------------
Image* ReadMetaDeformationField(const char* filename, Image::Format fmt)
{
    Image* def = NULL;
    if (IsMetaDeformationField(filename)) {
        def = ReadMetaImage(filename, DT_FLOAT, false, fmt);
    }
    return def;
}

// ---------------------------------------------------------------------------
inline static void print_mhd_entry(ostream& os, const char* key, float* v, int n)
{
    os << key << " =";
    for (int i = 0; i < n; i++) os << ' ' << v[i];
    os << '\n';
}

// ---------------------------------------------------------------------------
inline static void print_mhd_entry(ostream& os, const char* key, int* v, int n)
{
    os << key << " =";
    for (int i = 0; i < n; i++) os << ' ' << v[i];
    os << '\n';
}

// ---------------------------------------------------------------------------
bool WriteMetaImage(const char* filename, const Image* image, const Image* mask)
{
    assert(image != NULL);
    if (filename == NULL || filename[0] == '\0') {
        BASIS_THROW(invalid_argument, "No filename specified for output image!");
    }
    // known extensions to remove from input filename
    set<string> exts;
    exts.insert(".hdr");
    exts.insert(".hdr.gz");
    exts.insert(".img");
    exts.insert(".img.gz");
    exts.insert(".nii");
    exts.insert(".nii.gz");
    exts.insert(".nia");
    exts.insert(".nia.gz");
    exts.insert(".mhd");
    exts.insert(".raw");
    // filename of header and data file
    string bname = os::path::splitext(os::path::basename(filename), &exts)[0];
    string fname = bname + ".mhd";
    string iname = bname + ".raw";
    // make output directory if it does not exist
    if (!os::makedirs(os::path::dirname(filename))) return false;
    // convert NIfTI-1 header to Meta image header
    nifti_1_header hdr = image->hdr;
    if (image->imgfmt == Image::FORMAT_DRAMMS) hdr = PermuteXY(hdr);
    MetaHeader mhd = convert_nhdr2mhd(hdr);
    mhd.compressed = false; // image viewers such as ITK-SNAP do not support
                            // .mhd/.raw.gz ... so do not compress the data
    // write Meta image header
    ofstream ofs(fname.c_str());
    if (!ofs.is_open()) return false;
    ofs << "ObjectType = Image\n";
    ofs << "NDims = " << mhd.ndims << '\n';
    print_mhd_entry(ofs, "DimSize", mhd.size, mhd.ndims);
    ofs << "BinaryData = True\n";
    ofs << "BinaryDataByteOrderMSB = " << (mhd.msb ? "True" : "False") << '\n';
    ofs << "CompressedData = " << (mhd.compressed ? "True" : "False") << '\n';
    if (image->hdr.qform_code != NIFTI_XFORM_UNKNOWN) {
        /*

          TODO: The anatomical orientation written by this code does not match
                the one of the input Meta image. The first two letters are
                flipped. As long as it is not clear what the right orientation
                is to write, just leave out this header entry.

        mat44 qmat;
        for (int r = 0; r < mhd.ndims; r++) {
            for (int c = 0; c < mhd.ndims; c++) {
                qmat.m[r][c] = mhd.transform[r][c];
            }
            for (int c = mhd.ndims; c < 4; c++) {
                qmat.m[r][c] = (r == c) ? 1 : 0;
            }
        }
        for (int r = mhd.ndims; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                qmat.m[r][c] = (r == c) ? 1 : 0;
            }
        }
        int orient[3];
        nifti_mat44_to_orientation(qmat, &orient[0], &orient[1], &orient[2]);
        ofs << "AnatomicalOrientation = " << OrientationToString(orient[0], orient[1], orient[2]).c_str() << '\n';
        */
        ofs << "TransformMatrix =";
        for (int r = 0; r < mhd.ndims; r++) {
            for (int c = 0; c < mhd.ndims; c++) {
                ofs << ' ' << mhd.transform[r][c];
            }
        }
        ofs << '\n';
        ofs << "CenterOfRotation = 0 0 0\n";
    }
    print_mhd_entry(ofs, "Offset", mhd.offset, mhd.ndims);
    if (mhd.spacing [0] > 0) print_mhd_entry(ofs, "ElementSpacing", mhd.spacing,  mhd.ndims);
    if (mhd.elemsize[0] > 0) print_mhd_entry(ofs, "ElementSize",    mhd.elemsize, mhd.ndims);
    if (mhd.nchannels   > 1) ofs << "ElementNumberOfChannels = " << mhd.nchannels << '\n';
    ofs << "ElementType = " << meta_type_to_string(mhd.type) << '\n';
    ofs << "ElementDataFile = " << os::path::basename(iname) << '\n';
    ofs.close();
    // write raw image data
    Image tmp(image->img,
              image->region.nx, image->region.ny, image->region.nz,
              image->hdr.datatype, image->GetNumberOfComponents(),
              image->imgfmt,
              false); // the tmp image does not own the data
    tmp.compress = false; // therefore, we use the tmp image, to set compress
                          // to false as apparently some image viewers such as
                          // ITK-SNAP cannot deal with compressed Meta images.
    return WriteRawImage(iname.c_str(), &tmp, mask);
}

// ===========================================================================
// affine transformation
// ===========================================================================

// ---------------------------------------------------------------------------
int ReadAffineTransform(const char* filename, Image::Transform& transform, Image::Format fmt)
{
    // read 4x4 floating point values from text file
    Image::Transform m;
    FILE* fp = fopen(filename, "r");
    if (fp == NULL) return -1;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            if (fscanf(fp, "%f", &m.m[r][c]) != 1) {
                fclose(fp);
                return 0;
            }
        }
    }
    // check if there are no more non-whitespace characters in the file
    // otherwise, this was probably not really a 4x4 matrix stored in a text file
    int c = EOF;
    while ((c = fgetc(fp)) != EOF) {
        if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
            fclose(fp);
            return 0;
        }
    }
    if (ferror(fp)) {
        fclose(fp);
        return 0;
    }
    fclose(fp);
    // change output transform only after reading of matrix was successful
    if (fmt == Image::FORMAT_DRAMMS) {
        transform = PermuteXY(m);
    } else {
        transform = m;
    }
    return 1;
}

// ---------------------------------------------------------------------------
bool WriteAffineTransform(const char* filename, Image::Transform transform, Image::Format fmt)
{
    if (fmt == Image::FORMAT_DRAMMS) transform = PermuteXY(transform);
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) return false;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            if (c > 0) {
                if (fprintf(fp, " ") != 1) {
                    fclose(fp);
                    return false;
                }
            }
            if (fprintf(fp, "%f", transform.m[r][c]) < 1) {
                fclose(fp);
                return false;
            }
        }
        if (fprintf(fp, "\n") != 1) {
            fclose(fp);
            return false;
        }
    }
    fclose(fp);
    return true;
}


} // namespace dramms
