/**
 * @file  image.cxx
 * @brief Definition of image type and declaration of auxiliary image functions.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <cstdlib>     // abort()
#include <string.h>    // memset(), strncpy()
#include <nifti1_io.h> // nifti_make_new_header()

#include <basis/except.h> // BASIS_THROW, std::invalid_argument

#include "image.h"


// acceptable in .cxx file
using namespace std;


namespace dramms {


// NOTE: Function names in all lowercase and underscore (_) as word separator
//       denote functions that are only defined in this module and which are
//       not exposed to users of this image I/O module in the header file.


// ===========================================================================
// private helpers
// ===========================================================================

// ---------------------------------------------------------------------------
template <typename T>
static T*** allocate_memory(int x_size, int y_size, int z_size, int n, Image::Format fmt)
{
    if (fmt == Image::FORMAT_DRAMMS) {
        int tmp = x_size;
        x_size = y_size;
        y_size = tmp;
    }
    T*** p = new T**[z_size];
    for (int k = 0; k < z_size; k++) {
        p[k] = new T*[y_size];
        for (int j = 0; j < y_size; j++) {
            p[k][j] = new T[n * x_size];
            memset (p[k][j], 0, n * x_size * sizeof(T));
        }
    }
    return p;
}

// ---------------------------------------------------------------------------
template <typename T>
static void free_memory(T*** p, int x_size, int y_size, int z_size, Image::Format fmt)
{
    if (fmt == Image::FORMAT_DRAMMS) y_size = x_size;
    if (p) {
        for (int k = 0; k < z_size; k++) {
            for (int j = 0; j < y_size; j++) {
                delete [] p[k][j];
            }
            delete [] p[k];
        }
        delete [] p;
    }
}

// ===========================================================================
// construction / destruction
// ===========================================================================

// ---------------------------------------------------------------------------
Image::Image()
{
    memset(&hdr, 0, sizeof(hdr));
    filefmt  = FILE_FORMAT_NIFTI_1;
    imgfmt   = FORMAT_DRAMMS;
    ownsimg  = true;
    compress = true;
}

// ---------------------------------------------------------------------------
Image::Image(int x_dim, int y_dim, int z_dim, short datatype, int n, Format fmt)
{
    Init(Data(), x_dim, y_dim, z_dim, 1.0f, 1.0f, 1.0f, datatype, n, fmt, true);
}

// ---------------------------------------------------------------------------
Image::Image(int   x_dim,    int   y_dim,    int   z_dim,
             float x_pixdim, float y_pixdim, float z_pixdim,
             short datatype, int n,
             Image::Format fmt)
{
    Init(Data(), x_dim, y_dim, z_dim, x_pixdim, y_pixdim, z_pixdim, datatype, n, fmt, true);
}

// ---------------------------------------------------------------------------
Image::Image(Data   img,
             int    x_dim,
             int    y_dim,
             int    z_dim,
             short  datatype,
             int    n,
             Format fmt,
             bool   ownsimg)
{
    Init(img, x_dim, y_dim, z_dim, 1.0f, 1.0f, 1.0f, datatype, n, fmt, ownsimg);
}

// ---------------------------------------------------------------------------
Image::Image(Data img,
             int   x_dim,    int   y_dim,    int   z_dim,
             float x_pixdim, float y_pixdim, float z_pixdim,
             short datatype, int n,
             Image::Format fmt,
             bool ownsimg)
{
    Init(img, x_dim, y_dim, z_dim, x_pixdim, y_pixdim, z_pixdim, datatype, n, fmt, ownsimg);
}

// ---------------------------------------------------------------------------
Image::~Image()
{
    ReleaseData();
}

// ===========================================================================
// format of image representation
// ===========================================================================

// ---------------------------------------------------------------------------
void Image::SetFormat(const Format fmt)
{
    // do nothing if image format does not change
    if (imgfmt == fmt) return;
    // If either current or desired format is FORMAT_DRAMMS, we have to
    // permute the x and y axes in the header and each vector if this
    // image is a vector field (i.e., deformation field). Moreover, the
    // units of displacements in vector fields have to be converted.
    if (imgfmt == FORMAT_DRAMMS || fmt == FORMAT_DRAMMS) {
        // permute x and y components of vector field
        // and convert units of displacements
        const int N = hdr.dim[5]; // number of components per voxel
        if (N > 1) {
            float         tmp;
            vector<float> v(N);       // vector components
            vector<float> s(N, 1.0f); // scaling factors
            if (imgfmt == FORMAT_DRAMMS) {
                for (int n = 0; n < N; n++) {
                    s[n] = hdr.pixdim[n + 1] > 0.0 ? hdr.pixdim[n + 1] : 1.0f;
                }
            } else if (fmt == FORMAT_DRAMMS) {
                for (int n = 0; n < N; n++) {
                    s[n] = hdr.pixdim[n + 1] > 0.0 ? (1.0 / hdr.pixdim[n + 1]) : 1.0f;
                }
            }
            for (int k = 0; k < region.nz; k++) {
                for (int j = 0; j < region.ny; j++) {
                    for (int i = 0; i < region.nx; i++) {
                        // get vector components
                        get(i, j, k, v);
                        // convert units *before* permuting x and y
                        for (int n = 0; n < N; n++) v[n] *= s[n];
                        // swap x and y components
                        tmp  = v[0];
                        v[0] = v[1];
                        v[1] = tmp;
                        // set modified vector components
                        set(i, j, k, v);
                    }
                }
            }
        }
        // permute x and y in image header
        hdr = PermuteXY(hdr);
        // permute x and y of image region
        int nx    = region.nx;
        region.nx = region.ny;
        region.ny = nx;
    }
    // update image format
    imgfmt = fmt;
}

// ===========================================================================
// sub-image
// ===========================================================================

// ---------------------------------------------------------------------------
Image* Image::GetSlice(int k, bool copy)
{
    assert(0 <= k && k < region.nz);

    Image* slice = NULL;

    const int N = hdr.dim[5];

    if (copy) {
        slice = new Image(region.nx, region.ny, region.nz, hdr.datatype, N, imgfmt);
        if (slice) {
            vector<float> v;
            for (int j = 0; j < region.ny; j++) {
                for (int i = 0; i < region.nx; i++) {
                    get(i, j, k, v);
                    slice->set(i, j, v);
                }
            }
        }
    } else {
        Image::Data slice_img;
        switch (hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                slice_img.uc = &img.uc[k];
                break;
            case DT_SIGNED_SHORT:
                slice_img.ss = &img.ss[k];
                break;
            case DT_UINT16:
                slice_img.us = &img.us[k];
                break;
            case DT_SIGNED_INT:
                slice_img.si = &img.si[k];
                break;
            case DT_FLOAT:
                slice_img.fl = &img.fl[k];
                break;
        };
        slice = new Image(slice_img, region.nx, region.ny, 1, hdr.datatype, N, imgfmt, copy);
    }

    if (slice) {
        slice->CopyTransform  (this);
        slice->CopyDataScaling(this);
        slice->CopyMetaData   (this);
        // TODO: Adjust origin!
    }

    return slice;
}

// ---------------------------------------------------------------------------
const Image* Image::GetSlice(int k, bool copy) const
{
    // call non-const version of this method
    return const_cast<Image*>(this)->GetSlice(k, copy);
}

// ===========================================================================
// queries
// ===========================================================================

// ---------------------------------------------------------------------------
int Image::GetNumberOfComponents() const
{
    return hdr.dim[5] > 0 ? hdr.dim[5] : 1;
}

// ---------------------------------------------------------------------------
bool Image::HasSameSizeAs(const Image* rhs) const
{
    assert(rhs != NULL);
    if (imgfmt != rhs->imgfmt && (imgfmt == FORMAT_DRAMMS || rhs->imgfmt == FORMAT_DRAMMS)) {
        return region.nx == rhs->region.ny &&
               region.ny == rhs->region.nx &&
               region.nz == rhs->region.nz;
    } else {
        return region.nx == rhs->region.nx &&
               region.ny == rhs->region.ny &&
               region.nz == rhs->region.nz;
    }
}

// ---------------------------------------------------------------------------
bool Image::HasSmallerSizeAs(const Image* rhs) const
{
    assert(rhs != NULL);
    if (imgfmt != rhs->imgfmt && (imgfmt == FORMAT_DRAMMS || rhs->imgfmt == FORMAT_DRAMMS)) {
        return region.nx < rhs->region.ny &&
               region.ny < rhs->region.nx &&
               region.nz < rhs->region.nz;
    } else {
        return region.nx < rhs->region.nx &&
               region.ny < rhs->region.ny &&
               region.nz < rhs->region.nz;
    }
}

// ---------------------------------------------------------------------------
bool Image::HasSmallerOrSameSizeAs(const Image* rhs) const
{
    return HasSmallerSizeAs(rhs) || HasSameSizeAs(rhs);
}

// ===========================================================================
// coordinate transformation
// ===========================================================================

// ---------------------------------------------------------------------------
void Image::SetOrientation(int x_orient, int y_orient, int z_orient)
{
    string orient_code = OrientationToString(x_orient, y_orient, z_orient);
    if (orient_code.empty()) BASIS_THROW(invalid_argument, "Invalid orientation code!");
    // axial orientations
    if (orient_code == "LPI") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "LPS") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "LAS") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = 1.0;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "LAI") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 1.0;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "RPS") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 1.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "RPI") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 1.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "RAI") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 1.0;
    } else if (orient_code == "RAS") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 1.0;
    // sagittal orientations
    } else if (orient_code == "PIL") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = 0.5;
        hdr.quatern_c = 0.5;
        hdr.quatern_d = 0.5;
    } else if (orient_code == "PSL") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = -0.5;
        hdr.quatern_c = -0.5;
        hdr.quatern_d = 0.5;
    } else if (orient_code == "ASL") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = -0.5;
        hdr.quatern_c = 0.5;
        hdr.quatern_d = -0.5;
    } else if (orient_code == "AIL") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.5;
        hdr.quatern_c = -0.5;
        hdr.quatern_d = -0.5;
    } else if (orient_code == "PSR") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = -0.5;
        hdr.quatern_c = -0.5;
        hdr.quatern_d = 0.5;
    } else if (orient_code == "PIR") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.5;
        hdr.quatern_c = 0.5;
        hdr.quatern_d = 0.5;
    } else if (orient_code == "ASR") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = 0.5;
        hdr.quatern_c = -0.5;
        hdr.quatern_d = -0.5;
    } else if (orient_code == "AIR") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = -0.5;
        hdr.quatern_c = 0.5;
        hdr.quatern_d = -0.5;
    // coronal orientations
    } else if (orient_code == "LIP") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.70710676908493042;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "LSP") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = -0.70710676908493042;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "LSA") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = -0.70710676908493042;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "LIA") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.70710676908493042;
        hdr.quatern_c = 0.0;
        hdr.quatern_d = 0.0;
    } else if (orient_code == "RSP") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 0.70710676908493042;
        hdr.quatern_d = -0.70710676908493042;
    } else if (orient_code == "RIP") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 0.70710676908493042;
        hdr.quatern_d = -0.70710676908493042;
    } else if (orient_code == "RIA") {
        hdr.pixdim[0] = -1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 0.70710676908493042;
        hdr.quatern_d = 0.70710676908493042;
    } else if (orient_code == "RSA") {
        hdr.pixdim[0] = 1;
        hdr.quatern_b = 0.0;
        hdr.quatern_c = 0.70710676908493042;
        hdr.quatern_d = -0.70710676908493042;
    } else {
        BASIS_THROW(invalid_argument, "Unsupported orientation: " << orient_code)
    }
    hdr.qform_code = NIFTI_XFORM_ALIGNED_ANAT;
    // set sform to qform
    Transform sto_xyz = GetQFormTransform(hdr, imgfmt);
    memcpy(hdr.srow_x, sto_xyz.m[0], 4 * sizeof(float));
    memcpy(hdr.srow_y, sto_xyz.m[1], 4 * sizeof(float));
    memcpy(hdr.srow_z, sto_xyz.m[2], 4 * sizeof(float));
    hdr.sform_code = NIFTI_XFORM_SCANNER_ANAT;
    // update transformation matrices
    UpdateTransforms();
}

// ---------------------------------------------------------------------------
Image::Transform GetQFormTransform(Image::Header hdr, Image::Format fmt)
{
    Image::Transform qto_xyz;
    if (hdr.qform_code != NIFTI_XFORM_UNKNOWN) {
        qto_xyz = nifti_quatern_to_mat44(hdr.quatern_b,
                                         hdr.quatern_c,
                                         hdr.quatern_d,
                                         hdr.qoffset_x,
                                         hdr.qoffset_y,
                                         hdr.qoffset_z,
                                         hdr.pixdim[1],
                                         hdr.pixdim[2],
                                         hdr.pixdim[3],
                                         hdr.pixdim[0] < 0.0 ? -1.0 : 1.0);
    } else {
        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 4; c++) {
                if (r == c) {
                    if (r == 3) qto_xyz.m[r][c] = 1.0;
                    else        qto_xyz.m[r][c] = hdr.pixdim[r + 1];
                } else          qto_xyz.m[r][c] = 0.0;
            }
        }
        if (NIFTI_VERSION(hdr) == 0) {
            // ANALYZE 7.5 - radiological order RPI
            if (fmt == Image::FORMAT_DRAMMS) qto_xyz.m[1][1] *= -1;
            else                             qto_xyz.m[0][0] *= -1;
        } else {
            // NIfTI-1     - neurological order RAI
            qto_xyz.m[0][0] *= -1;
            qto_xyz.m[1][1] *= -1;
        }
    }
    return qto_xyz;
}

// ---------------------------------------------------------------------------
Image::Transform GetSFormTransform(Image::Header hdr, Image::Format fmt)
{
    Image::Transform sto_xyz;
    if (hdr.sform_code != NIFTI_XFORM_UNKNOWN) {
        memcpy(sto_xyz.m[0], hdr.srow_x, 4 * sizeof(float));
        memcpy(sto_xyz.m[1], hdr.srow_y, 4 * sizeof(float));
        memcpy(sto_xyz.m[2], hdr.srow_z, 4 * sizeof(float));
        sto_xyz.m[3][0] = sto_xyz.m[3][1] = sto_xyz.m[3][2] = 0; sto_xyz.m[3][3] = 1;
    } else {
        sto_xyz = GetQFormTransform(hdr, fmt);
    }
    return sto_xyz;
}

// ---------------------------------------------------------------------------
Image::LROrder GetLeftRightOrder(Image::Header hdr, Image::Format fmt)
{
    mat33 m;
    // qform
    Image::Transform qto_xyz = GetQFormTransform(hdr, Image::FORMAT_ITK);
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            m.m[r][c] = qto_xyz.m[r][c];
        }
    }
    float qdet = nifti_mat33_determ(m);
    // sform
    Image::Transform sto_xyz = GetSFormTransform(hdr, Image::FORMAT_ITK);
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            m.m[r][c] = sto_xyz.m[r][c];
        }
    }
    float sdet = nifti_mat33_determ(m);
    return ((sdet * qdet) < 0 || fabs(sdet * qdet) < 1.e-12)
                ? Image::LR_UNKNOWN // inconsistency, invalid transformation matrices
                : ((sdet < 0) ? Image::LR_RADIOLOGICAL : Image::LR_NEUROLOGICAL);
}

// ---------------------------------------------------------------------------
std::string LROrderToString(Image::LROrder order)
{
    if      (order == Image::LR_RADIOLOGICAL) return "Radiological";
    else if (order == Image::LR_NEUROLOGICAL) return "Neurological";
    else                                      return "Unknown";
}

// ---------------------------------------------------------------------------
void Image::UpdateTransforms()
{
    qto_xyz = GetQFormTransform(hdr, imgfmt);
    qto_ijk = nifti_mat44_inverse(qto_xyz);
}

// ===========================================================================
// copy image attributes
// ===========================================================================

// ---------------------------------------------------------------------------
void Image::CopyRegion(const Image* src)
{
    assert(src != NULL);
    assert(src->region.nx == region.nx);
    assert(src->region.ny == region.ny);
    assert(src->region.nz == region.nz);
    hdr.dim[1] = src->hdr.dim[1];
    hdr.dim[2] = src->hdr.dim[2];
    hdr.dim[3] = src->hdr.dim[3];
    region     = src->region;
}

// ---------------------------------------------------------------------------
void Image::CopyTransform(const Image* src, bool qform, bool sform)
{
    // attention: Do not just copy the src header fields as the default
    //            q/sform transformation also depends on the magic field,
    //            which, however, is not copied by this function.
    //            It is saver to actually set a q/sform transformation
    //            and therefore also a valid q/sform_code.
    assert(src != NULL);
    Image::Header shdr;
    if (src->imgfmt != imgfmt &&
            (src->imgfmt == FORMAT_DRAMMS || imgfmt == FORMAT_DRAMMS)) {
        shdr = PermuteXY(src->hdr);
    } else {
        shdr = src->hdr;
    }
    if (qform) {
        memcpy(hdr.pixdim, src->hdr.pixdim, 8 * sizeof(float));
        hdr.xyzt_units = src->hdr.xyzt_units;
        nifti_mat44_to_quatern(GetQFormTransform(src->hdr, src->imgfmt),
                               &hdr.quatern_b, &hdr.quatern_c, &hdr.quatern_d,
                               &hdr.qoffset_x, &hdr.qoffset_y, &hdr.qoffset_z,
                               NULL, NULL, NULL, &hdr.pixdim[0]);
        hdr.qform_code = src->hdr.qform_code;
        if (hdr.qform_code == NIFTI_XFORM_UNKNOWN) hdr.qform_code = NIFTI_XFORM_SCANNER_ANAT;
    }
    if (sform) {
        Transform sto_xyz = GetSFormTransform(src->hdr, src->imgfmt);
        memcpy(hdr.srow_x, sto_xyz.m[0], 4 * sizeof(float));
        memcpy(hdr.srow_y, sto_xyz.m[1], 4 * sizeof(float));
        memcpy(hdr.srow_z, sto_xyz.m[2], 4 * sizeof(float));
        hdr.sform_code = src->hdr.sform_code;
        if (hdr.sform_code == NIFTI_XFORM_UNKNOWN) hdr.sform_code = NIFTI_XFORM_SCANNER_ANAT;
    }
    UpdateTransforms();
}

// ---------------------------------------------------------------------------
void Image::CopyDataScaling(const Image* src)
{
    assert(src != NULL);
    hdr.scl_slope = src->hdr.scl_slope;
    hdr.scl_inter = src->hdr.scl_inter;
    hdr.cal_max   = src->hdr.cal_max;
    hdr.cal_min   = src->hdr.cal_min;
}

// ---------------------------------------------------------------------------
void Image::CopyMetaData(const Image* src)
{
    assert(src != NULL);
    filefmt         = src->filefmt;
    compress        = src->compress;
    hdr.intent_p1   = src->hdr.intent_p1;
    hdr.intent_p2   = src->hdr.intent_p2;
    hdr.intent_p3   = src->hdr.intent_p3;
    hdr.intent_code = src->hdr.intent_code;
    memcpy(hdr.intent_name, src->hdr.intent_name, 16 * sizeof(char));
    memcpy(hdr.descrip,     src->hdr.descrip,     80 * sizeof(char));
    memcpy(hdr.magic,       src->hdr.magic,        4 * sizeof(char));
}

// ===========================================================================
// auxiliary functions
// ===========================================================================

// ---------------------------------------------------------------------------
Image::FileFormat StringToFileFormat(const char* str)
{
    if      (strncmp(str, "ANALYZE", 8) == 0) return Image::FILE_FORMAT_ANALYZE;
    else if (strncmp(str, "NIFTI_1", 8) == 0) return Image::FILE_FORMAT_NIFTI_1;
    else if (strncmp(str, "NIFTI_2", 8) == 0) return Image::FILE_FORMAT_NIFTI_2;
    else if (strncmp(str, "NIFTI_A", 8) == 0) return Image::FILE_FORMAT_NIFTI_A;
    else if (strncmp(str, "META",    5) == 0) return Image::FILE_FORMAT_META;
    else                                      return Image::FILE_FORMAT_UNKNOWN;
}

// ---------------------------------------------------------------------------
const char* FileFormatToString(Image::FileFormat fmt)
{
    switch (fmt) {
        case Image::FILE_FORMAT_ANALYZE: return "ANALYZE";
        case Image::FILE_FORMAT_NIFTI_1: return "NIFTI_1";
        case Image::FILE_FORMAT_NIFTI_2: return "NIFTI_2";
        case Image::FILE_FORMAT_META:    return "META";
        default:                         return "UNKNOWN";
    };
}

// ---------------------------------------------------------------------------
Image::FileFormat NiftiTypeToFileFormat(Image::NiftiType type)
{
    return static_cast<Image::FileFormat>(type);
}

// ---------------------------------------------------------------------------
Image::NiftiType FileFormatToNiftiType(Image::FileFormat fmt)
{
    if (fmt < Image::FILE_FORMAT_META) return static_cast<Image::NiftiType>(fmt);
    else                               return Image::NIFTI_TYPE_UNKNOWN;
}

// ---------------------------------------------------------------------------
Image::Format StringToImageFormat(const char* str)
{
    if      (strncmp(str, "DRAMMS", 7) == 0) return Image::FORMAT_DRAMMS;
    else if (strncmp(str, "ITK",    4) == 0) return Image::FORMAT_ITK;
    else if (strncmp(str, "FSL",    4) == 0) return Image::FORMAT_FSL;
    else                                     return Image::FORMAT_UNKNOWN;
}

// ---------------------------------------------------------------------------
const char* ImageFormatToString(Image::Format fmt)
{
    switch (fmt) {
        case Image::FORMAT_DRAMMS: return "DRAMMS";
        case Image::FORMAT_ITK:    return "ITK";
        case Image::FORMAT_FSL:    return "FSL";
        default:                   return "UNKNOWN";
    };
}

// ---------------------------------------------------------------------------
Image::Transform PermuteXY(const Image::Transform T)
{
    Image::Transform t;
    t.m[0][0] = T.m[1][1];
    t.m[0][1] = T.m[1][0];
    t.m[0][2] = T.m[1][2];
    t.m[0][3] = T.m[1][3];
    t.m[1][0] = T.m[0][1];
    t.m[1][1] = T.m[0][0];
    t.m[1][2] = T.m[0][2];
    t.m[1][3] = T.m[0][3];
    t.m[2][0] = T.m[2][1];
    t.m[2][1] = T.m[2][0];
    t.m[2][2] = T.m[2][2];
    t.m[2][3] = T.m[2][3];
    t.m[3][0] = T.m[3][1];
    t.m[3][1] = T.m[3][0];
    t.m[3][2] = T.m[3][2];
    t.m[3][3] = T.m[3][3];
    return t;
}

// ---------------------------------------------------------------------------
Image::Header PermuteXY(const Image::Header& nhdr)
{
    Image::Header dhdr = nhdr;
    dhdr.dim   [1] = nhdr.dim   [2];
    dhdr.dim   [2] = nhdr.dim   [1];
    dhdr.pixdim[1] = nhdr.pixdim[2];
    dhdr.pixdim[2] = nhdr.pixdim[1];
    if (nhdr.qform_code != NIFTI_XFORM_UNKNOWN) {
        mat44 nqto_xyz = nifti_quatern_to_mat44(nhdr.quatern_b,
                                                nhdr.quatern_c,
                                                nhdr.quatern_d,
                                                nhdr.qoffset_x,
                                                nhdr.qoffset_y,
                                                nhdr.qoffset_z,
                                                nhdr.pixdim[1],
                                                nhdr.pixdim[2],
                                                nhdr.pixdim[3],
                                                nhdr.pixdim[0] < 0.0 ? -1.0 : 1.0);
        mat44 dqto_xyz = PermuteXY(nqto_xyz);
        nifti_mat44_to_quatern(dqto_xyz, &dhdr.quatern_b,
                                         &dhdr.quatern_c,
                                         &dhdr.quatern_d,
                                         &dhdr.qoffset_x,
                                         &dhdr.qoffset_y,
                                         &dhdr.qoffset_z,
                                         NULL, NULL, NULL,
                                         &dhdr.pixdim[0]);
    }
    if (nhdr.sform_code != NIFTI_XFORM_UNKNOWN) {
        for (int c = 0; c < 4; c++) {
            if (c == 0) {
                dhdr.srow_x[0] = nhdr.srow_y[1];
                dhdr.srow_y[0] = nhdr.srow_x[1];
                dhdr.srow_z[0] = nhdr.srow_z[1];
            } else if (c == 1) {
                dhdr.srow_x[1] = nhdr.srow_y[0];
                dhdr.srow_y[1] = nhdr.srow_x[0];
                dhdr.srow_z[1] = nhdr.srow_z[0];
            } else {
                dhdr.srow_x[c] = nhdr.srow_y[c];
                dhdr.srow_y[c] = nhdr.srow_x[c];
                dhdr.srow_z[c] = nhdr.srow_z[c];
            }
        }
    }
    return dhdr;
}

// ---------------------------------------------------------------------------
void GetImageOrientation(Image::Header hdr, int& x_orient, int& y_orient, int& z_orient, Image::Format hdrfmt)
{
    nifti_mat44_to_orientation(GetQFormTransform(hdr, hdrfmt), &x_orient, &y_orient, &z_orient);
}

// ---------------------------------------------------------------------------
int CharToOrientation(char c)
{
    switch (c) {
        case 'R': return NIFTI_R2L;
        case 'L': return NIFTI_L2R;
        case 'A': return NIFTI_A2P;
        case 'P': return NIFTI_P2A;
        case 'I': return NIFTI_I2S;
        case 'S': return NIFTI_S2I;
        default:  return 0;
    };
}

// ---------------------------------------------------------------------------
char OrientationToChar(int o)
{
    switch (o) {
        case NIFTI_R2L: return 'R';
        case NIFTI_L2R: return 'L';
        case NIFTI_A2P: return 'A';
        case NIFTI_P2A: return 'P';
        case NIFTI_I2S: return 'I';
        case NIFTI_S2I: return 'S';
        default:        return 'X';
    };
}

// ---------------------------------------------------------------------------
bool StringToOrientation(const char* str, int& x_orient, int& y_orient, int& z_orient)
{
    if (strlen(str) != 3) return false;
    x_orient = CharToOrientation(str[0]);
    y_orient = CharToOrientation(str[1]);
    z_orient = CharToOrientation(str[2]);
    return x_orient != 0 && y_orient != 0 && z_orient != 0;
}

// ---------------------------------------------------------------------------
string OrientationToString(int x_orient, int y_orient, int z_orient)
{
    string str;
    str.reserve(3);
    str += OrientationToChar(x_orient);
    str += OrientationToChar(y_orient);
    str += OrientationToChar(z_orient);
    return str;
}

// ===========================================================================
// initialization
// ===========================================================================

// ---------------------------------------------------------------------------
void Image::Init(Data   img,
                 int    x_dim,
                 int    y_dim,
                 int    z_dim,
                 float  x_pixdim,
                 float  y_pixdim,
                 float  z_pixdim,
                 short  datatype,
                 int    n,
                 Format fmt,
                 bool   ownsimg)
{
    if (x_dim <= 0 || y_dim <= 0 || z_dim <= 0) {
        BASIS_THROW(std::invalid_argument, "Invalid image size: " << x_dim
                << ", " << y_dim << ", " << z_dim);
    }
    if (n <= 0) {
        BASIS_THROW(std::invalid_argument, "Invalid number of components per voxel: " << n);
    }

    // set header information
    int dim[8];
    dim[0] = 5;      // number of dimensions
    dim[1] = x_dim;  // dimension of x axis
    dim[2] = y_dim;  // dimension of y axis
    dim[3] = z_dim;  // dimension of z axis
    dim[4] = 1;      // no time axis
    dim[5] = n;      // number of components per voxel
    dim[6] = 0;      // unused dimension
    dim[7] = 0;      // unused dimension

    Image::Header* hdr = nifti_make_new_header(dim, datatype);
    this->hdr = *hdr;
    free(hdr);

    this->hdr.pixdim[1]  = x_pixdim;
    this->hdr.pixdim[2]  = y_pixdim;
    this->hdr.pixdim[3]  = z_pixdim;
    this->hdr.qform_code = NIFTI_XFORM_SCANNER_ANAT;
    this->hdr.quatern_b  = 0; // default orientation is RAS
    this->hdr.quatern_c  = 0;
    this->hdr.quatern_d  = 1;

    this->hdr.intent_code = ((n > 1) ? NIFTI_INTENT_VECTOR  // vector
                                     : NIFTI_INTENT_NONE);  // scalar
    this->imgfmt   = fmt;
    this->compress = true;
    strncpy(this->hdr.magic, "n+1", 4);

    // set raw image data...
    if (img.uc != NULL) {
        this->img     = img;
        this->ownsimg = ownsimg;
    // ...or allocate memory
    } else {
        this->img = AllocateData(this->hdr.dim[1],
                                 this->hdr.dim[2],
                                 this->hdr.dim[3],
                                 this->hdr.dim[5],
                                 this->hdr.datatype,
                                 this->imgfmt);
        if (this->img.uc == NULL) {
            BASIS_THROW(runtime_error, "Failed to allocate memory for image!");
        }
        this->ownsimg = true;
    }

    // set image transformations
    UpdateTransforms();

    // set image region
    this->region.ox = 0;
    this->region.oy = 0;
    this->region.oz = 0;
    this->region.nx = this->hdr.dim[1];
    this->region.ny = this->hdr.dim[2];
    this->region.nz = this->hdr.dim[3];
}

// ---------------------------------------------------------------------------
Image::Data Image::AllocateData(int x_dim, int y_dim, int z_dim, int n, short datatype, Format fmt)
{
    Data img;
    if      (datatype == DT_UNSIGNED_CHAR) img.uc = allocate_memory<unsigned char> (x_dim, y_dim, z_dim, n, fmt);
    else if (datatype == DT_SIGNED_SHORT)  img.ss = allocate_memory<short>         (x_dim, y_dim, z_dim, n, fmt);
    else if (datatype == DT_UINT16)        img.us = allocate_memory<unsigned short>(x_dim, y_dim, z_dim, n, fmt);
    else if (datatype == DT_SIGNED_INT)    img.si = allocate_memory<int>           (x_dim, y_dim, z_dim, n, fmt);
    else if (datatype == DT_FLOAT)         img.fl = allocate_memory<float>         (x_dim, y_dim, z_dim, n, fmt);
    else BASIS_THROW(invalid_argument, "Unsupported datatype code: " << datatype);
    return img;
}

// ---------------------------------------------------------------------------
void Image::ReleaseData()
{
    if (ownsimg && img.uc != NULL) {
        switch (hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                free_memory(img.uc, region.nx, region.ny, region.nz, imgfmt);
                break;
            case DT_SIGNED_SHORT:
                free_memory(img.ss, region.nx, region.ny, region.nz, imgfmt);
                break;
            case DT_UINT16:
                free_memory(img.us, region.nx, region.ny, region.nz, imgfmt);
                break;
            case DT_SIGNED_INT:
                free_memory(img.si, region.nx, region.ny, region.nz, imgfmt);
                break;
            case DT_FLOAT:
                free_memory(img.fl, region.nx, region.ny, region.nz, imgfmt);
                break;
            default:
                ASSERT(false, "Invalid image datatype: " << hdr.datatype);
                abort();
        }
    }
}

// ===========================================================================
// debugging
// ===========================================================================

// ---------------------------------------------------------------------------
std::ostream& operator<<(ostream& os, const Image::Region& region)
{
    os << "image region with offset [";
    os << region.ox << ", " << region.oy;
    if (region.nz > 0) os << ", " << region.oz;
    os << "] and size " << region.nx << "x" << region.ny;
    if (region.nz > 0) os << "x" << region.nz;
    os << " voxels";
    return os;
}


} // namespace dramms
