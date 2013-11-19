/**
 * @file  image.hxx
 * @brief Definition of template and inline functions declared in image.h.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#pragma once
#ifndef _DRAMMS_COMMON_IMAGE_HXX
#define _DRAMMS_COMMON_IMAGE_HXX


#include <cstdlib> // abort()
#include <math.h>  // floor(), ceil()
#include <climits> // limits of datatypes

#include <basis/assert.h> // ASSERT()

#include "general.h" // round(), clamp()
#include "image.h"


namespace dramms {


// Attention: The order of the inline definitions matters and affects the
//            performance if wrong! Every (inline) method which is called by
//            another method has to be defined before this other method.
//            Otherwise, the compiler will not inline the object code.

// ===========================================================================
// get/set image value
// ===========================================================================

// ---------------------------------------------------------------------------
template <typename TValue>
TValue Image::value(int i, int j, int k, int n) const
{
	const int N = hdr.dim[5] > 0 ? hdr.dim[5] : 1;
    assert(0 <= i && i < region.nx);
    assert(0 <= j && j < region.ny);
    assert(0 <= k && k < region.nz);
    assert(0 <= n && n < N);
    if (imgfmt == FORMAT_DRAMMS) {
        switch (hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                return static_cast<TValue>(img.uc[k][i][j * N + n]);
            case DT_SIGNED_SHORT:
                return static_cast<TValue>(img.ss[k][i][j * N + n]);
            case DT_UINT16:
                return static_cast<TValue>(img.us[k][i][j * N + n]);
            case DT_SIGNED_INT:
                return static_cast<TValue>(img.si[k][i][j * N + n]);
            case DT_FLOAT:
                return static_cast<TValue>(img.fl[k][i][j * N + n]);
            default:
                ASSERT(false, "Image has unknown datatype: " << hdr.datatype);
                abort();
        }
    } else if (imgfmt == FORMAT_ITK || imgfmt == FORMAT_FSL) {
        switch (hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                return static_cast<TValue>(img.uc[k][j][i * N + n]);
            case DT_SIGNED_SHORT:
                return static_cast<TValue>(img.ss[k][j][i * N + n]);
            case DT_UINT16:
                return static_cast<TValue>(img.us[k][j][i * N + n]);
            case DT_SIGNED_INT:
                return static_cast<TValue>(img.si[k][j][i * N + n]);
            case DT_FLOAT:
                return static_cast<TValue>(img.fl[k][j][i * N + n]);
            default:
                ASSERT(false, "Image has unknown datatype: " << hdr.datatype);
                abort();
        }
    } else {
        ASSERT(false, "Unknown image format: " << imgfmt);
        abort();
    }
}

// ---------------------------------------------------------------------------
inline float Image::value(float i, float j, int k, int n) const
{
	if (i < 0 || i > region.nx - 1 ||
            j < 0 || j > region.ny - 1 ||
            k < 0 || k > region.nz - 1) {
        return 0.0f;
    }

    int iFloor = static_cast<int>(floor(i));
    int iCeil  = static_cast<int>(ceil(i));
    int jFloor = static_cast<int>(floor(j));
    int jCeil  = static_cast<int>(ceil(j));

    if (iFloor == iCeil) {
        if (iCeil != region.nx - 1) iCeil ++;
        else if (iFloor > 0)        iFloor--;
    }
    if (jFloor == jCeil) {
        if (jFloor != region.ny - 1) jCeil ++;
        else if (jFloor > 0)         jFloor--;
    }

    float a = i - iFloor; float A = 1.0f - a;
    float b = j - jFloor; float B = 1.0f - b;

    return A * B * value<float>(iFloor, jFloor, k, n) + 
           A * b * value<float>(iFloor, jCeil,  k, n) + 
           a * B * value<float>(iCeil,  jFloor, k, n) + 
           a * b * value<float>(iCeil,  jCeil,  k, n);
}

// ---------------------------------------------------------------------------
inline float Image::value(float i, float j, float k, int n) const
{
	if (i<0 && i>-0.5) i=0.0;
	if (j<0 && j>-0.5) j=0.0;
	if (k<0 && k>-0.5) k=0.0;
	if (i>region.nx-1 && i<region.nx-0.5) i=static_cast<float>(region.nx-1);
	if (j>region.ny-1 && j<region.ny-0.5) j=static_cast<float>(region.ny-1);
	if (k>region.nz-1 && k<region.nz-0.5) k=static_cast<float>(region.nz-1);
	
    if (i < 0 || i > region.nx - 1 ||
            j < 0 || j > region.ny - 1 ||
            k < 0 || k > region.nz - 1) {
        return 0.0f;
    }

    int iFloor = static_cast<int>(floor(i));
    int iCeil  = static_cast<int>(ceil (i));
    int jFloor = static_cast<int>(floor(j));
    int jCeil  = static_cast<int>(ceil (j));
    int kFloor = static_cast<int>(floor(k));
    int kCeil  = static_cast<int>(ceil (k));

    if (iFloor == iCeil) {
        if (iCeil != region.nx - 1) iCeil++;
        else if (iFloor > 0) iFloor--;
    }
    if (jFloor == jCeil) {
        if (jFloor != region.ny - 1) jCeil++;
        else if (jFloor > 0) jFloor--;
    }
    if (kFloor == kCeil) {
        if (kFloor != region.nz - 1) kCeil++;
        else if (kFloor > 0) kFloor--;
    }

    float a = i - iFloor; float A = 1.0f - a;
    float b = j - jFloor; float B = 1.0f - b;
    float c = k - kFloor; float C = 1.0f - c;

    return (C * (value<float>(iFloor, jFloor, kFloor, n) * (A * B) +
                 value<float>(iCeil,  jFloor, kFloor, n) * (a * B) +
                 value<float>(iFloor, jCeil,  kFloor, n) * (A * b) +
                 value<float>(iCeil,  jCeil,  kFloor, n) * (a * b)) +
            c * (value<float>(iFloor, jFloor, kCeil,  n) * (A * B) +
                 value<float>(iCeil,  jFloor, kCeil,  n) * (a * B) +
                 value<float>(iFloor, jCeil,  kCeil,  n) * (A * b) +
                 value<float>(iCeil,  jCeil,  kCeil,  n) * (a * b)));
}

// ---------------------------------------------------------------------------
inline float Image::get(int i, int j, int k, int n) const
{
    return value<float>(i, j, k, n);
}

// ---------------------------------------------------------------------------
template <typename TValue>
void Image::get(int i, int j, int k, std::vector<TValue>& v) const
{
    const int N = hdr.dim[5] > 0 ? hdr.dim[5] : 1;
    v.resize(N);
    for (int n = 0; n < N; n++) v[n] = value<TValue>(i, j, k, n);
}

// ---------------------------------------------------------------------------
template <typename TValue>
void Image::get(float i, float j, int k, std::vector<TValue>& v) const
{
    const int N = hdr.dim[5] > 0 ? hdr.dim[5] : 1;
    v.resize(N);
    for (int n = 0; n < N; n++) v[n] = static_cast<TValue>(value(i, j, k, n));
}

// ---------------------------------------------------------------------------
template <typename TValue>
void Image::get(float i, float j, std::vector<TValue>& v) const
{
    get(i, j, static_cast<int>(0), v);
}

// ---------------------------------------------------------------------------
template <typename TValue>
void Image::get(float i, float j, float k, std::vector<TValue>& v) const
{
    const int N = hdr.dim[5] > 0 ? hdr.dim[5] : 1;
    v.resize(N);
    for (int n = 0; n < N; n++) v[n] = static_cast<TValue>(value(i, j, k, n));
}

// ---------------------------------------------------------------------------
inline void Image::set(int i, int j, int k, int n, float v)
{
    const int N = hdr.dim[5] > 0 ? hdr.dim[5] : 1;
    assert(0 <= i && i < region.nx);
    assert(0 <= j && j < region.ny);
    assert(0 <= k && k < region.nz);
    assert(0 <= n && n < N);
    if (imgfmt == FORMAT_DRAMMS) {
        switch (hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                img.uc[k][i][j * N + n] = static_cast<unsigned char>(round(clamp(v, 0, UCHAR_MAX)));
                break;
            case DT_SIGNED_SHORT:
                img.ss[k][i][j * N + n] = static_cast<short>(round(clamp(v, SHRT_MIN, SHRT_MAX)));
                break;
            case DT_UINT16:
                img.us[k][i][j * N + n] = static_cast<unsigned short>(round(clamp(v, 0, USHRT_MAX)));
                break;
            case DT_SIGNED_INT:
                img.si[k][i][j * N + n] = static_cast<int>(round(clamp(v, INT_MIN, INT_MAX)));
                break;
            case DT_FLOAT:
                img.fl[k][i][j * N + n] = static_cast<float>(v);
                break;
            default:
                ASSERT(false, "Image has unknown datatype: " << hdr.datatype);
                abort();
        }
    } else if (imgfmt == FORMAT_ITK || imgfmt == FORMAT_FSL) {
        switch (hdr.datatype) {
            case DT_UNSIGNED_CHAR:
                img.uc[k][j][i * N + n] = static_cast<unsigned char>(round(clamp(v, 0, UCHAR_MAX)));
                break;
            case DT_SIGNED_SHORT:
                img.ss[k][j][i * N + n] = static_cast<short>(round(clamp(v, SHRT_MIN, SHRT_MAX)));
                break;
            case DT_UINT16:
                img.us[k][j][i * N + n] = static_cast<unsigned short>(round(clamp(v, 0, USHRT_MAX)));
                break;
            case DT_SIGNED_INT:
                img.si[k][j][i * N + n] = static_cast<int>(round(clamp(v, INT_MIN, INT_MAX)));
                break;
            case DT_FLOAT:
                img.fl[k][j][i * N + n] = static_cast<float>(v);
                break;
            default:
                ASSERT(false, "Image has unknown datatype: " << hdr.datatype);
                abort();
        }
    } else {
        ASSERT(false, "Unknown image format: " << imgfmt);
        abort();
    }
}

// ---------------------------------------------------------------------------
inline void Image::set(int i, int j, float v)
{
    set(i, j, 0, 0, v);
}

// ---------------------------------------------------------------------------
inline void Image::set(int i, int j, int k, float v)
{
    set(i, j, k, 0, v);
}

// ---------------------------------------------------------------------------
template <typename TValue>
void Image::set(int i, int j, int k, const std::vector<TValue>& v)
{
    const int N = hdr.dim[5] > 0 ? hdr.dim[5] : 1;
    assert(v.size() <= N);
    for (size_t n = 0; n < v.size(); n++) {
        set(i, j, k, n, static_cast<double>(v[n]));
    }
}

// ---------------------------------------------------------------------------
template <typename TValue>
void Image::set(int i, int j, const std::vector<TValue>& v)
{
    set(i, j, 0, v);
}

// ---------------------------------------------------------------------------
inline const Fvector3d& Image::v3(int i, int j, int k) const
{
    assert(0 <= i && i < region.nx);
    assert(0 <= j && j < region.ny);
    assert(0 <= k && k < region.nz);
    assert(hdr.dim[5] == 3);
    if      (imgfmt == FORMAT_DRAMMS)                      return img.v3[k][i][j];
    else if (imgfmt == FORMAT_ITK || imgfmt == FORMAT_FSL) return img.v3[k][j][i];
    else {
        ASSERT(false, "Unknown image format: " << imgfmt);
        abort();
    }
}

// ===========================================================================
// coordinate transformation
// ===========================================================================

// ---------------------------------------------------------------------------
inline Fvector3d Image::to_xyz(const Fvector3d& ijk) const
{
    Fvector3d xyz;
    xyz.x = ijk.x * qto_xyz.m[0][0] + ijk.y * qto_xyz.m[0][1] + ijk.z * qto_xyz.m[0][2] + qto_xyz.m[0][3];
    xyz.y = ijk.x * qto_xyz.m[1][0] + ijk.y * qto_xyz.m[1][1] + ijk.z * qto_xyz.m[1][2] + qto_xyz.m[1][3];
    xyz.z = ijk.x * qto_xyz.m[2][0] + ijk.y * qto_xyz.m[2][1] + ijk.z * qto_xyz.m[2][2] + qto_xyz.m[2][3];
    return xyz;
}

// ---------------------------------------------------------------------------
inline Fvector3d Image::to_xyz(float i, float j, float k) const
{
    Fvector3d ijk;
    ijk.x = i;
    ijk.y = j;
    ijk.z = k;
    return to_xyz(ijk);
}

// ---------------------------------------------------------------------------
inline Fvector3d Image::to_xyz(const Ivector3d& ijk) const
{
    Fvector3d cijk;
    cijk.x = static_cast<float>(ijk.x);
    cijk.y = static_cast<float>(ijk.y);
    cijk.z = static_cast<float>(ijk.z);
    return to_xyz(cijk);
}

// ---------------------------------------------------------------------------
inline Fvector3d Image::to_xyz(int i, int j, int k) const
{
    Fvector3d cijk;
    cijk.x = static_cast<float>(i);
    cijk.y = static_cast<float>(j);
    cijk.z = static_cast<float>(k);
    return to_xyz(cijk);
}

// ---------------------------------------------------------------------------
inline Fvector3d Image::to_cijk(const Fvector3d& xyz) const
{
    Fvector3d cijk;
    cijk.x = xyz.x * qto_ijk.m[0][0] + xyz.y * qto_ijk.m[0][1] + xyz.z * qto_ijk.m[0][2] + qto_ijk.m[0][3];
    cijk.y = xyz.x * qto_ijk.m[1][0] + xyz.y * qto_ijk.m[1][1] + xyz.z * qto_ijk.m[1][2] + qto_ijk.m[1][3];
    cijk.z = xyz.x * qto_ijk.m[2][0] + xyz.y * qto_ijk.m[2][1] + xyz.z * qto_ijk.m[2][2] + qto_ijk.m[2][3];
    return cijk;
}

// ---------------------------------------------------------------------------
inline Fvector3d Image::to_cijk(float x, float y, float z) const
{
    Fvector3d xyz;
    xyz.x = x;
    xyz.y = y;
    xyz.z = z;
    return to_cijk(xyz);
}

// ---------------------------------------------------------------------------
inline Ivector3d Image::to_ijk(const Fvector3d& xyz) const
{
    Fvector3d cijk = to_cijk(xyz);
    Ivector3d ijk;
    ijk.x = static_cast<int>(round(cijk.x));
    ijk.y = static_cast<int>(round(cijk.y));
    ijk.z = static_cast<int>(round(cijk.z));
    return ijk;
}

// ---------------------------------------------------------------------------
inline Ivector3d Image::to_ijk(float x, float y, float z) const
{
    Fvector3d xyz;
    xyz.x = x;
    xyz.y = y;
    xyz.z = z;
    return to_ijk(xyz);
}


} // namespace dramms


#endif // _DRAMMS_COMMON_IMAGE_HXX
