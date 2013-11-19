/**
 * @file  utilities.cxx
 * @brief Utility imaging functions.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#include <stack>
#include <limits>         // numeric_limits
#include <math.h>         // fabs(), copysignf()

#include <basis/assert.h> // ASSERT()
#include <basis/except.h> // invalid_argument

#include "general.h"      // round()
#include "matrix.h"       // Mat_Inverse()

#include "utilities.h"


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
// update scl_slope and scl_inter entries of image header
static void update_header_scl(Image* image_out, const Image* image_in, float slope, float inter)
{
    // set scl_slope and scl_inter of output image
    // Attention: image_in and image_out may reference the same image!
    if (image_in->hdr.scl_slope == 0.0f) {
        image_out->hdr.scl_slope = 1.0f / slope;
        image_out->hdr.scl_inter = inter / slope;
    } else {
        image_out->hdr.scl_slope = image_in->hdr.scl_slope / slope;
        image_out->hdr.scl_inter = image_in->hdr.scl_inter - inter / slope;
    }
    // if slope is 1 and intercept is 0, set both to 0
    if (fabs(image_out->hdr.scl_slope - 1.0f) < 1.0e-9f && fabs(image_out->hdr.scl_inter) < 1.0e-9f) {
        image_out->hdr.scl_slope = 0.0f;
        image_out->hdr.scl_inter = 0.0f;
    }
}

// ---------------------------------------------------------------------------
// warp image, casting interpolated values to integral type
template <typename T>
static void warp_intensity_image(T***         warpedimage,
                                 const Image* image,
                                 const Image* deffield,
                                 bool         interpolate)
{
    const int x_size_in = image->region.nx;
    const int y_size_in = image->region.ny;
    const int z_size_in = image->region.nz;

    const int x_size_out = deffield->region.nx;
    const int y_size_out = deffield->region.ny;
    const int z_size_out = deffield->region.nz;

    assert(deffield->hdr.dim[5]   == 3);
    assert(deffield->hdr.datatype == DT_FLOAT);

    if (deffield->imgfmt != Image::FORMAT_DRAMMS) {
        BASIS_THROW(runtime_error, "Deformation field must be in DRAMMS format!");
    }

    Fvector3d*** v3 = deffield->img.v3;

    float  x,  y,  z;
    int    ii, jj, kk;

    for (int k = 0; k < z_size_out; k++) {
        for (int i = 0; i < x_size_out; i++) {
            for (int j = 0; j < y_size_out; j++) {
                x = static_cast<float>(i) + v3[k][i][j].x;
                y = static_cast<float>(j) + v3[k][i][j].y;
                z = static_cast<float>(k) + v3[k][i][j].z;

                if (interpolate) {
                    warpedimage->set(i, j, k, image->value(x, y, z));
                } else {
                    ii = static_cast<int>(round(x));
                    jj = static_cast<int>(round(y));
                    kk = static_cast<int>(round(z));

                    if (ii < 0 || ii >= x_size_in ||
                            jj < 0 || jj >= y_size_in ||
                            kk < 0 || kk >= z_size_in) {
                        warpedimage->set(i, j, k, 0);
                    } else {
                        warpedimage->set(i, j, k, image->get(ii, jj, kk));
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
inline
static int grow_region(const Image* image,
                       float        threshold,
                       Image*       regions,
                       Image*       visited,
                       Ivector3d    seed,
                       int          label,
                       bool         withinslice = false)
{
    int              size = 0; // size of region
    stack<Ivector3d> active;   // active boundary voxels
    Ivector3d        voxel;    // current voxel
    Ivector3d        neighbor; // neighboring voxel

    assert(image   != NULL);
    assert(regions != NULL);
    assert(visited != NULL);
    assert(regions->region == image->region);
    assert(visited->region == image->region);
    assert(0 <= seed.x && seed.x < image->region.nx);
    assert(0 <= seed.y && seed.y < image->region.ny);
    assert(0 <= seed.z && seed.z < image->region.nz);
    assert(label > 0);

    active.push(seed);
    visited->set(seed.x, seed.y, seed.z, 1);

    do {
        voxel = active.top();
        active.pop();
        regions->set(voxel.x, voxel.y, voxel.z, label);
        size++;

        int minx = (voxel.x - 1 >= 0               ) ? (voxel.x - 1) : voxel.x;
        int maxx = (voxel.x + 1 <  image->region.nx) ? (voxel.x + 1) : voxel.x;
        int miny = (voxel.y - 1 >= 0               ) ? (voxel.y - 1) : voxel.y;
        int maxy = (voxel.y + 1 <  image->region.ny) ? (voxel.y + 1) : voxel.y;

        int minz = voxel.z;
        int maxz = voxel.z;
        if (!withinslice) {
            minz = (voxel.z - 1 >= 0               ) ? (voxel.z - 1) : voxel.z;
            maxz = (voxel.z + 1 <  image->region.nz) ? (voxel.z + 1) : voxel.z;
        }

        for (neighbor.z = minz; neighbor.z <= maxz; neighbor.z++) {
            for (neighbor.y = miny; neighbor.y <= maxy; neighbor.y++) {
                for (neighbor.x = minx; neighbor.x <= maxx; neighbor.x++) {
                    if (image->get(neighbor.x, neighbor.y, neighbor.z) > threshold &&
                            visited->get(neighbor.x, neighbor.y, neighbor.z) == 0) {
                        active.push(neighbor);
                        visited->set(neighbor.x, neighbor.y, neighbor.z, 1);
                    }
                }
            }
        }
    } while (!active.empty());

    return size;
}

// ===========================================================================
// intensity scaling
// ===========================================================================

// ---------------------------------------------------------------------------
void GetDatatypeRange(short datatype, float& min, float& max)
{
    switch (datatype) {
        case DT_UNSIGNED_CHAR:
            min = static_cast<float>(numeric_limits<unsigned char>::min());
            max = static_cast<float>(numeric_limits<unsigned char>::max());
            break;
        case DT_SIGNED_SHORT:
            min = static_cast<float>(numeric_limits<short>::min());
            max = static_cast<float>(numeric_limits<short>::max());
            break;
        case DT_UINT16:
            min = static_cast<float>(numeric_limits<unsigned short>::min());
            max = static_cast<float>(numeric_limits<unsigned short>::max());
            break;
        case DT_SIGNED_INT:
            min = static_cast<float>(numeric_limits<int>::min());
            max = static_cast<float>(numeric_limits<int>::max());
            break;
        case DT_UINT32:
            min = static_cast<float>(numeric_limits<unsigned int>::min());
            max = static_cast<float>(numeric_limits<unsigned int>::max());
            break;
        case DT_FLOAT:
            min = - numeric_limits<float>::max();
            max = + numeric_limits<float>::max();
            break;
        default:
            ASSERT(false, "Unsupported datatype: " << datatype);
            abort();
    }
}

// ---------------------------------------------------------------------------
bool IsByteImage(const Image* image, bool scale)
{
    if (image->GetNumberOfComponents() > 1) return false;
    float slope = 1.0f;
    float inter = 0.0f;
    if (image->hdr.scl_slope != 0) {
        if (scale) {
            slope = image->hdr.scl_slope;
            inter = image->hdr.scl_inter;
        } else {
            if (image->hdr.scl_slope != 1.0f && image->hdr.scl_inter != 0.0f) return false;
        }
    }
    float min = 256.0f;
    float max = -1.0f;
    float I;
    for (int k = 0; k < image->region.nz; k++) {
        for (int j = 0; j < image->region.ny; j++) {
            for (int i = 0; i < image->region.nx; i++) {
                I = image->get(i, j, k) * slope + inter;
                if (I != floor(I)) return false;
                if (I < min) min = I;
                if (I > max) max = I;
            }
        }
    }
    return min <= max && 0.0f <= min && max <= 255.0f;
}

// ---------------------------------------------------------------------------
void GetDataScaling(const Image* image, float& slope, float& inter, float lmin, float lmax, bool eff)
{
    // if no (valid) output range specified, use scaling function of image header
    if (lmax <= lmin) {
        if (image->hdr.scl_slope == 0.0f) {
            slope = 1.0f;
            inter = 0.0f;
        } else {
            slope = image->hdr.scl_slope;
            inter = image->hdr.scl_inter;
        }
    // otherwise, determine optimal scaling, optionally avoiding 'long tail' phenomenon
    } else {
        // get intensity range of image data
        float min, max;
        GetIntensityRange(image, min, max, false, eff);
        // scaling function to map [lmin, lmax] to [min, max]
        //
        // Note that for the rescaling, the *inverse* of the function
        // "slope * I + inter" is applied! The return values of this function
        // are the data scaling as if it had been applied to the data before
        // it was being saved and in order to recover the "original" intensities,
        // these would need to be scaled by the inverse of this function.
        //
        // If the output datatype is signed, keep intercept fixed to 0 such that
        // zero intensities are preserved and only intensities to the left or right
        // of zero are stretched/shrinked. This rescaling preserves the sign of the
        // intensities. Otherwise, if the output datatype is unsigned, also shift
        // the intensities as needed to not truncate any values.
        if (lmin < 0) {
            slope = std::min(fabs(lmin) / fabs(min),
                             fabs(lmax) / fabs(max));
            inter = 0;
        } else {
            slope = (lmax - lmin) / (max - min);
            inter = lmin - min * slope;
        }
    }
}

/**
 * @brief Rescale image intensities.
 *
 * This function rescales the intensities of the input image using the linear
 * function @c y = @p slope * @c x + @p inter and writes the rescaled intensities
 * to the given output image. If a value @c y is outside the range of values
 * which can be represented by the output data type, these values are truncated
 * and mapped to the closest representable value.
 *
 * @param [out] scaled Image with rescaled intensities.
 * @param [in]  image  Input image.
 * @param [in]  slope  Rescaling slope.
 * @param [in]  inter  Rescaling intercept.
 */
void ScaleImage(Image* scaled, const Image* image, float slope, float inter)
{
    // check arguments
    if (scaled->HasSmallerSizeAs(image)) {
        BASIS_THROW(invalid_argument, "Output image is smaller than input image!");
    }
    const int N = image->GetNumberOfComponents();
    if (scaled->GetNumberOfComponents() < N) {
        BASIS_THROW(invalid_argument, "Output image has less components than input image!");
    }
    if (slope == 0) {
        slope = image->hdr.scl_slope;
        inter = image->hdr.scl_inter;
    }
    if (slope == 0) return; // according to definition of scl_slope of NIfTI-1
    // rescale intensities
    for (int k = 0; k < image->region.nz; k++) {
        for (int j = 0; j < image->region.ny; j++) {
            for (int i = 0; i < image->region.nx; i++) {
                for (int n = 0; n < N; n++) {
                    scaled->set(i, j, k, n, image->get(i, j, k, n) * slope + inter);
                }
            }
        }
    }
    // update scl_slope and scl_inter of output image
    update_header_scl(scaled, image, slope, inter);
}

// ===========================================================================
// copy / conversion
// ===========================================================================

// ---------------------------------------------------------------------------
Image* CastImage(const Image* image, short datatype, bool scale)
{
    // allocate new image
    Image* casted_image = new Image(image->region.nx,
                                    image->region.ny,
                                    image->region.nz,
                                    datatype,
                                    image->hdr.dim[5] > 0
                                        ? image->hdr.dim[5]
                                        : 1,
                                    image->imgfmt);
    if (!casted_image) return NULL;
    // copy image attributes
    casted_image->CopyRegion     (image);
    casted_image->CopyTransform  (image);
    casted_image->CopyDataScaling(image);
    casted_image->CopyMetaData   (image);
    // copy image data and optionally apply scaling
    float slope = 1.0f, inter = 0.0f;
    if (scale) GetDataScaling(image, slope, inter);
    ScaleImage(casted_image, image, slope, inter);
    return casted_image;
}

// ---------------------------------------------------------------------------
Image* CastImage(const Image* image, short datatype, float min, float max, bool smooth, bool eff)
{
    // allocate new image
    Image* casted_image = new Image(image->region.nx,
                                    image->region.ny,
                                    image->region.nz,
                                    datatype,
                                    image->hdr.dim[5] > 0
                                        ? image->hdr.dim[5]
                                        : 1,
                                    image->imgfmt);
    if (!casted_image) return NULL;
    // copy image attributes
    casted_image->CopyRegion     (image);
    casted_image->CopyTransform  (image);
    casted_image->CopyDataScaling(image);
    casted_image->CopyMetaData   (image);
    casted_image->compress = image->compress;
    // determine scaling function
    float slope = 1.0f, inter = 0.0f;
    // default [min, max] intensity range if not specified and output
    // datatype is not a floating point type
    if (min >= max && datatype != DT_FLOAT) {
        if (datatype == image->hdr.datatype) {
            GetIntensityRange(image, min, max);
        } else {
            GetDatatypeRange(datatype, min, max);
        }
    }
    // scale intensities only if a desired output range is specified or
    // the intensities will be rescaled to the maximum output range of
    // the output datatype; otherwise, truncation might occur
    GetDataScaling(image, slope, inter, min, max, eff);
    // rescale intensities of smoothed image
    if (smooth) {
        Fvector3d pixdim;
        Ivector3d radius;
        Fvector3d sigma;
        pixdim.x = image->hdr.pixdim[1];
        pixdim.y = image->hdr.pixdim[2];
        pixdim.z = image->hdr.pixdim[3];
        GetGaussParameters(pixdim, radius, sigma);
        Image* gauss = Gauss(radius, sigma);
        if (gauss == NULL) {
            delete casted_image;
            return NULL;
        }
        Image* smoothed_image = ConvolveImage(image, gauss);
        delete gauss;
        if (smoothed_image == NULL) {
            delete casted_image;
            return NULL;
        }
        for (int k = 0; k < image->region.nz; k++) {
            for (int j = 0; j < image->region.ny; j++) {
                for (int i = 0; i < image->region.nx; i++) {
                    for (int n = 0; n < image->GetNumberOfComponents(); n++) {
                        smoothed_image->set(i, j, k, n, 0.5 * (smoothed_image->get(i, j, k, n) + image->get(i, j, k, n)));
                    }
                }
            }
        }
        ScaleImage(casted_image, smoothed_image, slope, inter);
        delete smoothed_image;
    // otherwise, rescale intensities of image itself only
    } else {
        ScaleImage(casted_image, image, slope, inter);
    }
    return casted_image;
}

// ---------------------------------------------------------------------------
Image* CopyImage(const Image* image)
{
    return CastImage(image, image->hdr.datatype, 1.0f, 0.0f);
}

// ===========================================================================
// image statistics
// ===========================================================================

// ---------------------------------------------------------------------------
void GetIntensityRange(const Image* image, float& min, float& max, bool scale, bool eff)
{
    const int N = image->GetNumberOfComponents();
    // get minimum and maximum intensity values
    min = + numeric_limits<float>::max();
    max = - numeric_limits<float>::max();
    float I;
    for (int k = 0; k < image->region.nz; k++) {
        for (int j = 0; j < image->region.ny; j++) {
            for (int i = 0; i < image->region.nx; i++) {
                for (int n = 0; n < N; n++) {
                    I = image->get(i, j, k, n);
                    if (I < min) min = I;
                    if (I > max) max = I;
                }
            }
        }
    }
    // invalid image...
    if (min > max) {
        min = max = 0.0f;
        return;
    } 
    // adjust intensity range to include "effective" intensities only
    // looking at the histogram of image intensities
    if (max > min && eff) {
        // determine if image is actually a byte image
        const bool byte_image = IsByteImage(image, true);
        // parameters of histogram computation
        const int   nbins = 256;
        // const float hmin  = byte_image ?   0 : min;
        // const float hmax  = byte_image ? 255 : max;
		const float hmin = min;
		const float hmax = max;
        const float width = (hmax - hmin) / (nbins - 1);
        // compute normalized cumulative histogram
        float* histo = ComputeHistogram(image, nbins, false, hmin, hmax);
        if (histo == NULL) {
            BASIS_THROW(runtime_error, "Failed to allocate memory for histogram!");
        }
        histo[0] = 0.0f; // ignore background
        CumulateHistogram(NormalizeHistogram(histo, nbins), nbins);
		
		// Yangming added on 4/18/2013 for dramms-v1.3
		// ----------------------
		int binIndex_threshold = FindAdaptiveThreshold(histo, nbins);
		min = hmin + width * binIndex_threshold; // min used to be 0 before grouping long tails; now min is at least the effective threshold.
		DRAMMS_MSG("effective threshold between background and foreground = " << min << " (" << binIndex_threshold << "-th bin out of " << nbins << ").");
		float prc1  = 0.01 * (1.0-histo[binIndex_threshold]) + histo[binIndex_threshold];
		float prc2  = 0.02 * (1.0-histo[binIndex_threshold]) + histo[binIndex_threshold];
		float prc98 = 0.98 * (1.0-histo[binIndex_threshold]) + histo[binIndex_threshold];
		float prc99 = 0.99 * (1.0-histo[binIndex_threshold]) + histo[binIndex_threshold];        
		// ----------------------

        // determine bins below which certain percentiles of intensities fall
        int prctile1   = 0;
        int prctile2   = 0;
        int prctile98  = nbins - 1;
        int prctile99  = nbins - 1;
        int prctile100 = nbins - 1;
		for (int i = MAX(1,binIndex_threshold); i < prctile100; i++) {
            if (histo[i - 1] <= prc1  && histo[i    ] >= prc1 ) prctile1   = i;
            if (histo[i - 1] <= prc2  && histo[i    ] >= prc2 ) prctile2   = i;
            if (histo[i    ] <= prc98 && histo[i + 1] >= prc98) prctile98  = i + 1;
            if (histo[i    ] <= prc99 && histo[i + 1] >= prc99) prctile99  = i + 1;
            if (histo[i] == histo[nbins - 1])                 prctile100 = i;

        } 
        // destroy no longer needed histogram
        delete [] histo;
        // adjust input range to cover effective range of intensities only
        if (prctile2 > (50.0f / 255.0f) * (nbins - 1)) {
            min = MAX(min, hmin + width * (prctile1 + prctile2) / 5.0);
        }
        if (prctile100 == prctile99 && prctile99 == prctile98) {
            max = hmin + prctile98 * width;
        } else if (prctile100-prctile98 > (50.0f / 255.0f) * (nbins - 1)) {
            max = hmin + width * (9*prctile98 + 3*prctile99 + prctile100) / 13.0;
        }
        // some verbose output
        DRAMMS_MSG("raw intensity range (after smoothing):");
        DRAMMS_MSG("     prctile( 0, 1,  2) = " << hmin                        << ", " << (hmin + prctile1  * width) << ", " << (hmin + prctile2   * width));
        DRAMMS_MSG("     prctile(98,99,100) = " << (hmin + prctile98 * width)  << ", " << (hmin + prctile99 * width) << ", " << (hmin + prctile100 * width));
        DRAMMS_MSG("effective intensity range (after grouping long, flat tails in both ends):");
        DRAMMS_MSG("     min = " << min << ", max = " << max);
    }
    // scale intensities according to scaling function if requested
    if (scale && image->hdr.scl_slope != 0) {
        min = min * image->hdr.scl_slope + image->hdr.scl_inter;
        max = max * image->hdr.scl_slope + image->hdr.scl_inter;
    }
}

// ---------------------------------------------------------------------------
float* ComputeHistogram(const Image* image, int nbins, bool scale, float min, float max)
{
    // actual formula to get bin of voxel:
    // bin = (I * slope + inter) / width
    // where
    // width = ((slope * max + inter) - (slope * min + inter)) / (nbins - 1)
    //       = slope * (max -  min) / (nbins - 1)
    //
    // this can be simplified to
    // bin = I / width + inter / slope
    // where
    // width = (max - min) / (nbins - 1)

    float* histo = new float[nbins];
    if (histo == NULL) return NULL;

    if (min >= max) GetIntensityRange(image, min, max);

    float slope = 1;
    float inter = 0;
    if (scale && image->hdr.scl_slope != 0) {
        slope = image->hdr.scl_slope;
        inter = image->hdr.scl_inter;
    }

    const float width  = (max - min) / (nbins - 1); // bins centered around min and max, respectively
    const float offset = inter / slope;

    for (int i = 0; i < nbins; i++) histo[i] = 0.0f;

    for (int k = 0; k < image->region.nz; k++) {
        for (int j = 0; j < image->region.ny; j++) {
            for (int i = 0; i < image->region.nx; i++) {
                int bin = static_cast<int>(round((image->get(i, j, k)-min) / width + offset));
                if      (bin < 0        ) histo[0  ] += 1; // - 0.5         -> -1    -> 0
                else if (bin > nbins - 1) histo[255] += 1; // (nbins - 0.5) -> nbins -> nbins - 1
                else                      histo[bin] += 1;
            }
        }
    }

    return histo;
}

// ---------------------------------------------------------------------------
float* NormalizeHistogram(float* histo, int nbins)
{
    float sum = 0.0f;
    for (int i = 0; i < nbins; i++) sum += histo[i];
    for (int i = 1; i < nbins; i++) histo[i] /= sum;
    return histo;
}

// ---------------------------------------------------------------------------
float* CumulateHistogram(float* histo, int nbins)
{
    for (int i = 1; i < nbins; i++) histo[i] += histo[i - 1];
    return histo;
}


// ---------------------------------------------------------------------------
// Yangming added on 4/18/2013 for dramms-v1.3
int FindAdaptiveThreshold(float* histo, int nbins)
{
	int ind;
	int ind_bin=12;
	// for (ind = 1; ind < 15; ind++)
		// if (histo[ind]<=0.5) ind_bin = ind; // adpatively set the default threshold (can be as big as 12 in [0,255] range, as long as the accumulative histogram does not excedes 0.5)
	
	//for (ind = 36; ind > 0; ind=ind-4) // the biggest possible threshold is the 28-th bin. // Yangming changed on 11/01/2013. Instead of testing a proper intensity threshold from higher intensity to lower ones, which perfer a relatively high intensity threshold, we chose a more conservative approach --- testing threshold values from low to high intensity, preferring a relatively lower intensity to keep almost all internal structures, at the risk of keeping some moderate level of background level.
	for (ind = 0; ind < 28; ind=ind+4) // the biggest possible threshold is the 28-th bin.
		{
		// printf("testing thresholding at %d-th bin\n", ind);
		// printf("histo(%d)=%f\n", ind, histo[ind]);
		// printf("histo(%d)=%f\n", ind+4, histo[ind+4]);
		// printf("histo(%d)=%f\n", ind+8, histo[ind+8]);
		// printf("histo(%d)=%f\n\n", ind+12, histo[ind+12]);
		
		//if ( ((histo[ind+12]-histo[ind])<0.065) && ((histo[ind+8]-histo[ind])<0.04) && ((histo[ind+4]-histo[ind])<0.02) && (histo[ind]<0.8) )  // Yangming changed on 11/01/2013, because the one used here, although completely and desirably removes background noise, could remove some internal low-intensity structures. Therefore, we chose a more conservative approach --- choose a low threshold that keeps almost all internal structures at the cost of possibly keeping some moderate level of background noise
		if ( ((histo[ind+12]-histo[ind])<0.075) && ((histo[ind+8]-histo[ind])<0.07) && ((histo[ind+4]-histo[ind])<0.025) && (histo[ind]<0.6) )
			{
			//ind_bin = ind+2;
			ind_bin = ind;
			break;
			}
		}
		
	// printf("decide to have threshold at the %d-th bin\n\n", ind_bin);
	return ind_bin;
}


// ===========================================================================
// convolution
// ===========================================================================

// ---------------------------------------------------------------------------
void GetGaussParameters(int size, int& radius, float& sigma)
{
    if (size < 80) {
        radius = 1;
        sigma  = 0.7;
    } else if ((size >= 80) && (size <= 160)) {
        radius = 2;
        sigma  = 1.2;
    } else if ((size >= 160) && (size <= 260)) {
        radius = 1;
        sigma  = 0.4;
    } else {
        radius = 2;
        sigma  = 0.5;
    }
}

// ---------------------------------------------------------------------------
void GetGaussParameters(Fvector2d& pixdim, Ivector2d& radius, Fvector2d& sigma)
{
    float min = std::min(pixdim.x, pixdim.y);

    Fvector2d ratio;
    ratio.x = pixdim.x / min;
    ratio.y = pixdim.y / min;

    radius.x = radius.y = 2;
    sigma .x = sigma .y = 0.6 / std::min(ratio.x, ratio.y);
}

// ---------------------------------------------------------------------------
void GetGaussParameters(Fvector3d& pixdim, Ivector3d& radius, Fvector3d& sigma)
{
    float min = std::min(std::min(pixdim.x, pixdim.y), pixdim.z);

    Fvector3d ratio;
    ratio.x = pixdim.x / min;
    ratio.y = pixdim.y / min;
    ratio.z = pixdim.z / min;

    radius.x = radius.y = 2;
    sigma .x = sigma .y = 0.6 / std::min(ratio.x, ratio.y);

    radius.z = 1;
    sigma .z = 0.25 / ratio.z;
}

// ---------------------------------------------------------------------------
float Gauss(float x, float sigma)
{
    return exp(- pow(x, 2.0f) / (2.0f * pow(sigma, 2.0f))) / (sigma * sqrt(2.0f * G_PI));
}

// ---------------------------------------------------------------------------
Image* Gauss(const Ivector2d& radius, const Fvector2d& sigma, Image::Format fmt)
{
    Ivector2d size;
    size.x = 2 * radius.x + 1;
    size.y = 2 * radius.y + 1;

    Image* gauss = new Image(size.x, size.y, 1, DT_FLOAT, 1, fmt);

    float sum = 0.0f;
    for (int j = 0; j < size.y; j++) {
        for (int i = 0; i < size.x; i++) {
            gauss->set(i, j, Gauss(i - radius.x, sigma.x) * Gauss(j - radius.y, sigma.y));
            sum += gauss->get(i, j);
        }
    }

    for (int j = 0; j < size.y; j++) {
        for (int i = 0; i < size.x; i++) {
            gauss->set(i, j, gauss->get(i, j) / sum);
        }
    }

    return gauss;
}

// ---------------------------------------------------------------------------
Image* Gauss(const Ivector3d& radius, const Fvector3d& sigma, Image::Format fmt)
{
    Ivector3d size;
    size.x = 2 * radius.x + 1;
    size.y = 2 * radius.y + 1;
    size.z = 2 * radius.z + 1;

    Image* gauss = new Image(size.x, size.y, size.z, DT_FLOAT, 1, fmt);

    float sum = 0.0f;
    for (int k = 0; k < size.z; k++) {
        for (int j = 0; j < size.y; j++) {
            for (int i = 0; i < size.x; i++) {
                gauss->set(i, j, k,
                        Gauss(i - radius.x, sigma.x)
                        * Gauss(j - radius.y, sigma.y)
                        * Gauss(k - radius.z, sigma.z));
                sum += gauss->get(i, j, k);
            }
        }
    }

    for (int k = 0; k < size.z; k++) {
        for (int j = 0; j < size.y; j++) {
            for (int i = 0; i < size.x; i++) {
                gauss->set(i, j, k, gauss->get(i, j, k) / sum);
            }
        }
    }

    return gauss;
}

// ---------------------------------------------------------------------------
void ConvolveImage(const Image* image, const Image* kernel, Image* output)
{
    assert(image != NULL);
    assert(kernel != NULL);
    assert(output != NULL);

    if (output->HasSmallerSizeAs(image)) {
        BASIS_THROW(invalid_argument, "Output image too small to store result of convolution!");
    }
    if (kernel->GetNumberOfComponents() > 1) {
        BASIS_THROW(invalid_argument, "Convolution kernel must be scalar!");
    }

    Ivector3d size;   // size of input image
    Ivector3d offset; // 0 if kernel size is odd, and -1 otherwise
    Ivector3d radius; // kernel radius

    size.x = image->region.nx;
    size.y = image->region.ny;
    size.z = image->region.nz;

    radius.x = (kernel->region.nx - 1) / 2;
    radius.y = (kernel->region.ny - 1) / 2;
    radius.z = (kernel->region.nz - 1) / 2;
    if (radius.x < 0) radius.x = 0;
    if (radius.y < 0) radius.y = 0;
    if (radius.z < 0) radius.z = 0;

    offset.x = kernel->region.nx > 0 ? (- (kernel->region.nx - 1) % 2) : 0;
    offset.y = kernel->region.ny > 0 ? (- (kernel->region.ny - 1) % 2) : 0;
    offset.z = kernel->region.nz > 0 ? (- (kernel->region.nz - 1) % 2) : 0;

    const int N = image->GetNumberOfComponents();

    vector<double> v(N);
    for (int k = 0; k < size.z; k++) {
        for (int j = 0; j < size.y; j++) {
            for (int i = 0; i < size.x; i++) {
                for (int n = 0; n < N; n++) {
                    v[n] = 0.0f;
                    for (int kk = - offset.z - radius.z; kk <= radius.z; kk++) {

                        int kkk = k + kk;
                        if      (kkk < 0)       kkk = - kkk - 1;
                        else if (kkk >= size.z) kkk = (size.z - 1) - (kkk - size.z);

                        for (int jj = - offset.y - radius.y; jj <= radius.y; jj++) {

                            int jjj = j + jj;
                            if      (jjj < 0)       jjj = - jjj - 1;
                            else if (jjj >= size.y) jjj = (size.y - 1) - (jjj - size.y);

                            for (int ii = - offset.x - radius.x; ii <= radius.x; ii++) {

                                int iii = i + ii;
                                if      (iii < 0)       iii = - iii - 1;
                                else if (iii >= size.x) iii = (size.x - 1) - (iii - size.x);

                                v[n] += image->get(iii, jjj, kkk, n)
                                        * kernel->get(offset.x + radius.x + ii,
                                                      offset.y + radius.y + jj,
                                                      offset.z + radius.z + kk);
                            }
                        }
                    }
                }
                output->set(i, j, k, v);
            }
        }
    }
}

// ---------------------------------------------------------------------------
Image* ConvolveImage(const Image* image, const Image* kernel)
{
    assert(image != NULL);
    // allocate output image
    Image* output = new Image(image->region.nx,
                              image->region.ny,
                              image->region.nz,
                              DT_FLOAT,
                              image->GetNumberOfComponents(),
                              image->imgfmt);
    if (output == NULL) return NULL;
    output->CopyRegion     (image);
    output->CopyTransform  (image);
    output->CopyDataScaling(image);
    output->CopyMetaData   (image);
    // convolve image
    ConvolveImage(image, kernel, output);
    return output;
}

// ===========================================================================
// image regions
// ===========================================================================

// ---------------------------------------------------------------------------
Image::Region GetForegroundRegion(const Image* image, float threshold)
{
    int min_i = 0;
    int min_j = 0;
    int min_k = 0;
    int max_i = image->region.nx - 1;
    int max_j = image->region.ny - 1;
    int max_k = image->region.nz - 1;

    // minium k index
    for (int k = min_k; k <= max_k; k++) {
        for (int j = min_j; j <= max_j; j++) {
            for (int i = min_i; i <= max_i; i++) {
                if (image->get(i, j, k) > threshold) {
                    min_k = k;
                    // break out of ALL loops
                    j = max_j + 1;
                    k = max_k + 1;
                    break;
                }
            }
        }
    }

    // maximum k index
    for (int k = max_k; k >= min_k; k--) {
        for (int j = min_j; j <= max_j; j++) {
            for (int i = min_i; i <= max_i; i++) {
                if (image->get(i, j, k) > threshold) {
                    max_k = k;
                    // break out of ALL loops
                    j = max_j + 1;
                    k = min_k - 1;
                    break;
                }
            }
        }
    }
 
    // minium j index
    for (int j = min_j; j <= max_j; j++) {
        for (int k = min_k; k <= max_k; k++) {
            for (int i = min_i; i <= max_i; i++) {
                if (image->get(i, j, k) > threshold) {
                    min_j = j;
                    // break out of ALL loops
                    k = max_k + 1;
                    j = max_j + 1;
                    break;
                }
            }
        }
    }

    // maximum j index
    for (int j = max_j; j >= min_j; j--) {
        for (int k = min_k; k <= max_k; k++) {
            for (int i = min_i; i <= max_i; i++) {
                if (image->get(i, j, k) > threshold) {
                    max_j = j;
                    // break out of ALL loops
                    k = max_k + 1;
                    j = min_j - 1;
                    break;
                }
            }
        }
    }

    // minium i index
    for (int i = min_i; i <= max_i; i++) {
        for (int k = min_k; k <= max_k; k++) {
            for (int j = min_j; j <= max_j; j++) {
                if (image->get(i, j, k) > threshold) {
                    min_i = i;
                    // break out of ALL loops
                    k = max_k + 1;
                    i = max_i + 1;
                    break;
                }
            }
        }
    }

    // maximum k index
    for (int i = max_i; i >= min_i; i--) {
        for (int k = min_k; k <= max_k; k++) {
            for (int j = min_j; j <= max_j; j++) {
                if (image->get(i, j, k) > threshold) {
                    max_i = i;
                    // break out of ALL loops
                    k = max_k + 1;
                    i = min_i - 1;
                    break;
                }
            }
        }
    }

    Image::Region region;
    region.ox = image->region.ox + min_i;
    region.oy = image->region.oy + min_j;
    region.oz = image->region.oz + min_k;
    region.nx = max_i - min_i + 1;
    region.ny = max_j - min_j + 1;
    region.nz = max_k - min_k + 1;
    return region;
}

// ---------------------------------------------------------------------------
Image::Region JoinRegions(const Image::Region& region1, const Image::Region& region2)
{
    Image::Region region = region1;

    int maxx  = region .ox + region .nx;
    int maxy  = region .oy + region .ny;
    int maxz  = region .oz + region .nz;
    int maxx2 = region2.ox + region2.nx;
    int maxy2 = region2.oy + region2.ny;
    int maxz2 = region2.oz + region2.nz;

    if (region2.ox < region.ox) region.ox = region2.ox;
    if (region2.oy < region.oy) region.oy = region2.oy;
    if (region2.oz < region.oz) region.oz = region2.oz;

    if (maxx2 > maxx) region.nx = maxx2 - region.ox + 1;
    if (maxy2 > maxy) region.ny = maxy2 - region.oy + 1;
    if (maxz2 > maxz) region.nz = maxz2 - region.oz + 1;

    return region;
}

// ===========================================================================
// image grid manipulation
// ===========================================================================

// ---------------------------------------------------------------------------
bool PermuteAxes(Image* image, int order[3])
{
    assert(image != NULL);
    // memorize image format
    Image::Format fmt = image->imgfmt;
    // convert order to FORMAT_ITK
    if (fmt == Image::FORMAT_DRAMMS) {
        int tmp[3];
        tmp[0] = order[1];
        tmp[1] = order[0];
        tmp[2] = order[2];
        for (int i = 0; i < 3; i++) {
            if      (tmp[i] == 0) order[i] = 1;
            else if (tmp[i] == 1) order[i] = 0;
            else                  order[i] = tmp[i];
        }
    }
    // nothing to do if order is not permuted
    if (order[0] == 0 && order[1] == 1 && order[2] == 2) return true;
    // check permutation
    if (order[0] < 0 || order[0] > 3 ||
            order[1] < 0 || order[1] > 3 ||
            order[2] < 0 || order[2] > 3 ||
            order[0] == order[1] ||
            order[0] == order[2] ||
            order[1] == order[2]) {
        return false;
    }
    // convert image to FORMAT_ITK
    image->SetFormat(Image::FORMAT_ITK);
    // invert order
    int inverse_order[3];
    inverse_order[order[0]] = 0;
    inverse_order[order[1]] = 1;
    inverse_order[order[2]] = 2;
    // permute sizes
    int   in_dim    [3] = {image->region.nx,     image->region.ny,     image->region.nz};
    float in_pixdim [3] = {image->hdr.pixdim[1], image->hdr.pixdim[2], image->hdr.pixdim[3]};
    int   out_dim   [3] = {in_dim   [order[0]], in_dim   [order[1]], in_dim   [order[2]]};
    float out_pixdim[3] = {in_pixdim[order[0]], in_pixdim[order[1]], in_pixdim[order[2]]};
    // permute orientation matrices
    Image::Transform in_qto_xyz = GetQFormTransform(image);
    Image::Transform in_sto_xyz = GetSFormTransform(image);
    Image::Transform out_qto_xyz;
    Image::Transform out_sto_xyz;
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            out_qto_xyz.m[r][c] = in_qto_xyz.m[r][order[c]];
            out_sto_xyz.m[r][c] = in_sto_xyz.m[r][order[c]];
        }
        // origin does not change by a permute
        out_qto_xyz.m[r][3] = in_qto_xyz.m[r][3];
        out_sto_xyz.m[r][3] = in_sto_xyz.m[r][3];
    }
    for (int c = 0; c < 4; c++) {
        out_qto_xyz.m[3][c] = in_qto_xyz.m[3][c];
        out_sto_xyz.m[3][c] = in_sto_xyz.m[3][c];
    }
    // allocate output image
    Image* permuted_image = new Image(out_dim[0],
                                      out_dim[1],
                                      out_dim[2],
                                      out_pixdim[0],
                                      out_pixdim[1],
                                      out_pixdim[2],
                                      image->hdr.datatype,
                                      image->hdr.dim[5] > 0 ? image->hdr.dim[5] : 1,
                                      image->imgfmt);
    if (permuted_image == NULL) {
        image->SetFormat(fmt); // restore image format
        return false;
    }
    permuted_image->CopyDataScaling(image);
    permuted_image->CopyMetaData   (image);
    // copy image region - see ResizeImage()
    int in_offset [3] = {image->region.ox,    image->region.oy,    image->region.oz};
    int out_offset[3] = {in_offset[order[0]], in_offset[order[1]], in_offset[order[2]]};
    permuted_image->region.ox  = out_offset[0];
    permuted_image->region.oy  = out_offset[1];
    permuted_image->region.oz  = out_offset[2];
    permuted_image->hdr.dim[1] = image->hdr.dim[order[0] + 1];
    permuted_image->hdr.dim[2] = image->hdr.dim[order[1] + 1];
    permuted_image->hdr.dim[3] = image->hdr.dim[order[2] + 1];
    // set transformations
    nifti_mat44_to_quatern(out_qto_xyz,
                           &permuted_image->hdr.quatern_b,
                           &permuted_image->hdr.quatern_c,
                           &permuted_image->hdr.quatern_d,
                           &permuted_image->hdr.qoffset_x,
                           &permuted_image->hdr.qoffset_y,
                           &permuted_image->hdr.qoffset_z,
                           NULL, NULL, NULL,
                           &permuted_image->hdr.pixdim[0]);
    if (image->hdr.qform_code == NIFTI_XFORM_UNKNOWN) {
        permuted_image->hdr.qform_code = NIFTI_XFORM_SCANNER_ANAT;
    } else {
        permuted_image->hdr.qform_code = image->hdr.qform_code;
    }
    if (image->hdr.sform_code != NIFTI_XFORM_UNKNOWN) {
        for (int c = 0; c < 4; c++) {
            permuted_image->hdr.srow_x[c] = out_sto_xyz.m[0][c];
            permuted_image->hdr.srow_y[c] = out_sto_xyz.m[1][c];
            permuted_image->hdr.srow_z[c] = out_sto_xyz.m[2][c];
        }
        permuted_image->hdr.sform_code = image->hdr.sform_code;
    }
    // permute image data
    int idx[3];
    vector<float> v;
    for (idx[2] = 0; idx[2] < permuted_image->region.nz; idx[2]++) {
        for (idx[1] = 0; idx[1] < permuted_image->region.ny; idx[1]++) {
            for (idx[0] = 0; idx[0] < permuted_image->region.nx; idx[0]++) {
                image->get(idx[inverse_order[0]],
                           idx[inverse_order[1]],
                           idx[inverse_order[2]], v);
                permuted_image->set(idx[0], idx[1], idx[2], v);
            }
        }
    }
    // swap images and free input image
    Image::Header tmphdr    = image->hdr;
    Image::Region tmpregion = image->region;
    Image::Data   tmpimg    = image->img;
    image->region           = permuted_image->region;
    image->hdr              = permuted_image->hdr;
    image->img              = permuted_image->img;
    permuted_image->region  = tmpregion;
    permuted_image->hdr     = tmphdr;
    permuted_image->img     = tmpimg;
    delete permuted_image;
    // restore image format
    image->SetFormat(fmt);
    return true;
}

// ---------------------------------------------------------------------------
bool FlipAxes(Image* image, bool flip[3], bool update_offset)
{
    assert(image != NULL);
    // memorize image format
    Image::Format fmt = image->imgfmt;
    // convert flip to FORMAT_ITK
    if (fmt == Image::FORMAT_DRAMMS) {
        bool tmp = flip[0];
        flip[0]  = flip[1];
        flip[1]  = tmp;
    }
    // nothing to do if no axis has to be flipped
    if (!flip[0] && !flip[1] && !flip[2]) return true;
    // convert image to FORMAT_ITK
    image->SetFormat(Image::FORMAT_ITK);
    // allocate output image
    Image* flipped_image = new Image(image->region.nx,
                                     image->region.ny,
                                     image->region.nz,
                                     image->hdr.datatype,
                                     image->hdr.dim[5] > 0 ? image->hdr.dim[5] : 1,
                                     image->imgfmt);
    if (flipped_image == NULL) {
        image->SetFormat(fmt); // restore image format
        return false;
    }
    flipped_image->CopyTransform  (image);
    flipped_image->CopyDataScaling(image);
    flipped_image->CopyMetaData   (image);
    flipped_image->CopyRegion     (image);
    // transformation matrix encoding the flips
    Image::Transform flipmat = IdentityTransform();
    for (int i = 0; i < 3; i++) flipmat.m[i][i] = flip[i] ? -1 : 1;
    // change orientation
    Image::Transform in_qto_xyz  = GetQFormTransform(image);
    Image::Transform in_sto_xyz  = GetSFormTransform(image);
    Image::Transform out_qto_xyz = ConcatenateTransforms(flipmat, in_qto_xyz);
    Image::Transform out_sto_xyz = ConcatenateTransforms(flipmat, in_sto_xyz);
    // change offset by transforming the index of the voxel which
    // will end up at voxel position (0,0,0) to physical units
    if (update_offset) {
        int i = flip[0] ? (image->region.nx - 1) : 0;
        int j = flip[1] ? (image->region.ny - 1) : 0;
        int k = flip[2] ? (image->region.nz - 1) : 0;
        for (int r = 0; r < 3; r++) {
            out_qto_xyz.m[r][3] = i * in_qto_xyz.m[r][0] +
                                  j * in_qto_xyz.m[r][1] +
                                  k * in_qto_xyz.m[r][2] +
                                      in_qto_xyz.m[r][3];
        }
        for (int r = 0; r < 3; r++) {
            out_sto_xyz.m[r][3] = i * in_sto_xyz.m[r][0] +
                                  j * in_sto_xyz.m[r][1] +
                                  k * in_sto_xyz.m[r][2] +
                                      in_sto_xyz.m[r][3];
        }
    }
    // update header information of flipped image
    if (flipped_image->hdr.qform_code == NIFTI_XFORM_UNKNOWN) {
        flipped_image->hdr.qform_code = NIFTI_XFORM_SCANNER_ANAT;
    }
    nifti_mat44_to_quatern(out_qto_xyz,
                           &flipped_image->hdr.quatern_b,
                           &flipped_image->hdr.quatern_c,
                           &flipped_image->hdr.quatern_d,
                           &flipped_image->hdr.qoffset_x,
                           &flipped_image->hdr.qoffset_y,
                           &flipped_image->hdr.qoffset_z,
                           NULL, NULL, NULL,
                           &flipped_image->hdr.pixdim[0]);
    if (flipped_image->hdr.sform_code != NIFTI_XFORM_UNKNOWN) {
        for (int c = 0; c < 4; c++) {
            flipped_image->hdr.srow_x[c] = out_sto_xyz.m[0][c];
            flipped_image->hdr.srow_y[c] = out_sto_xyz.m[1][c];
            flipped_image->hdr.srow_z[c] = out_sto_xyz.m[2][c];
        }
    }
    // flip image data
    vector<float> v;
    int d_in_i = flip[0] ? -1 : 1;
    int d_in_j = flip[1] ? -1 : 1;
    int d_in_k = flip[2] ? -1 : 1;
    int in_k = flip[2] ? (image->region.nz - 1) : 0;
    for (int k = 0; k < image->region.nz; k++, in_k += d_in_k) {
        int in_i = flip[0] ? (image->region.nx - 1) : 0;
        for (int i = 0; i < image->region.nx; i++, in_i += d_in_i) {
            int in_j = flip[1] ? (image->region.ny - 1) : 0;
            for (int j = 0; j < image->region.ny; j++, in_j += d_in_j) {
                image->get(in_i, in_j, in_k, v);
                flipped_image->set(i, j, k, v);
            }
        }
    }
    // swap images and free input image
    Image::Header tmphdr    = image->hdr;
    Image::Region tmpregion = image->region;
    Image::Data   tmpimg    = image->img;
    image->region = flipped_image->region;
    image->hdr    = flipped_image->hdr;
    image->img    = flipped_image->img;
    flipped_image->region = tmpregion;
    flipped_image->hdr    = tmphdr;
    flipped_image->img    = tmpimg;
    delete flipped_image;
    // restore image format
    image->SetFormat(fmt);
    return true;
}

// ---------------------------------------------------------------------------
bool OrientImage(Image* image, int orient[3], bool update_offset)
{
    assert(image != NULL);
    // memorize image format
    Image::Format fmt = image->imgfmt;
    // convert orient to FORMAT_ITK
    if (fmt == Image::FORMAT_DRAMMS) {
        int tmp = orient[0];
        orient[0] = orient[1];
        orient[1] = tmp;
    }
    // convert image to FORMAT_ITK
    image->SetFormat(Image::FORMAT_ITK);
    // get current image orientation
    int in_orient[3];
    GetImageOrientation(image->hdr, in_orient[0], in_orient[1], in_orient[2], image->imgfmt);
    // determine required permutations
    int order[3] = {0, 1, 2};
    for (int i = 0; i < 3; i++) {
        int in_dir = 0;
        if      (in_orient[i] == NIFTI_R2L || in_orient[i] == NIFTI_L2R) in_dir = NIFTI_L2R;
        else if (in_orient[i] == NIFTI_A2P || in_orient[i] == NIFTI_P2A) in_dir = NIFTI_P2A;
        else if (in_orient[i] == NIFTI_I2S || in_orient[i] == NIFTI_S2I) in_dir = NIFTI_I2S;
        assert(in_dir != 0);
        for (int j = 0; j < 3; j++) {
            int out_dir = 0;
            if      (orient[j] == NIFTI_R2L || orient[j] == NIFTI_L2R) out_dir = NIFTI_L2R;
            else if (orient[j] == NIFTI_A2P || orient[j] == NIFTI_P2A) out_dir = NIFTI_P2A;
            else if (orient[j] == NIFTI_I2S || orient[j] == NIFTI_S2I) out_dir = NIFTI_I2S;
            if (out_dir == 0) {
                image->SetFormat(fmt); // restore image format
                return false;
            }
            if (out_dir == in_dir) {
                order[j] = i;
                break;
            }
        }
    }
    // determine required flips
    bool flip[3] = {false, false, false};
    for (int i = 0; i < 3; i++) {
        if (orient[i] != in_orient[order[i]]) flip[i] = true;
    }
    // permute vector components and change sign if required
    if (image->hdr.dim[5] > 1 &&
            (order[0] != 0 || order[1] != 1 || order[2] != 2 ||
             flip[0] || flip[1] || flip[2])) {
        if (image->hdr.dim[5] == 3) {
            float d[3];
            for (int k = 0; k < image->region.nz; k++) {
                for (int i = 0; i < image->region.nx; i++) {
                    for (int j = 0; j < image->region.ny; j++) {
                        d[0] = image->img.v3[k][i][j].x;
                        d[1] = image->img.v3[k][i][j].y;
                        d[2] = image->img.v3[k][i][j].z;
                        if (flip[0]) image->img.v3[k][i][j].x = - d[order[0]];
                        else         image->img.v3[k][i][j].x =   d[order[0]];
                        if (flip[1]) image->img.v3[k][i][j].y = - d[order[1]];
                        else         image->img.v3[k][i][j].y =   d[order[1]];
                        if (flip[2]) image->img.v3[k][i][j].z = - d[order[2]];
                        else         image->img.v3[k][i][j].z =   d[order[2]];
                    }
                }
            }
        } else {
            image->SetFormat(fmt); // restore image format
            BASIS_THROW(std::invalid_argument, "OrientImage() function only supports three-dimensional scalar images and deformation fields.");
        }
    }
    // permute and flip axes as required
    DRAMMS_MSG("Change image orientation from "
            << OrientationToString(fmt == Image::FORMAT_DRAMMS ? in_orient[1] : in_orient[0],
                                   fmt == Image::FORMAT_DRAMMS ? in_orient[0] : in_orient[1],
                                   in_orient[2]) << " to "
            << OrientationToString(orient[0], orient[1], orient[2]));
    DRAMMS_MSG("  Permutation: " << order[0] << order[1] << order[2]);
    DRAMMS_MSG("  Flips:       " << (flip[0] ? "y" : "n") << (flip[1] ? "y" : "n") << (flip[2] ? "y" : "n"));
    bool ok = PermuteAxes(image, order) && FlipAxes(image, flip, update_offset);
    // restore image format
    image->SetFormat(fmt);
    return ok;
}

// ---------------------------------------------------------------------------
Image* ResizeImage(const Image* image, Image::Region region, float pad_value, bool subdomain)
{
    // extend region in each invalid dimension
    if (region.nx <= 0) {
        region.ox = 0;
        region.nx = image->hdr.dim[1];
    }
    if (region.ny <= 0) {
        region.oy = 0;
        region.ny = image->hdr.dim[2];
    }
    if (region.nz <= 0) {
        region.oz = 0;
        region.nz = image->hdr.dim[3];
    }
    // allocate resized image
    Image* resized_image = new Image(region.nx, region.ny, region.nz,
                                     image->hdr.datatype,
                                     image->GetNumberOfComponents(),
                                     image->imgfmt);
    if (resized_image == NULL) return NULL;
    // copy image attributes
    resized_image->CopyTransform  (image, true, false); // qform only
    resized_image->CopyDataScaling(image);
    resized_image->CopyMetaData   (image);
    // update image offset
    if (resized_image->hdr.qform_code != NIFTI_XFORM_UNKNOWN) {
        resized_image->hdr.qoffset_x += (region.ox - image->region.ox) * image->hdr.pixdim[1];
        resized_image->hdr.qoffset_y += (region.oy - image->region.oy) * image->hdr.pixdim[2];
        resized_image->hdr.qoffset_z += (region.oz - image->region.oz) * image->hdr.pixdim[3];
    }
    // copy image data within region, padding where necessary
    vector<float> v;
    vector<float> pad_v(image->GetNumberOfComponents(), pad_value);
    for (int k = region.oz; k < (region.oz + region.nz); k++) {
        if (image->region.oz <= k && k < (image->region.oz + image->region.nz)) {
            for (int j = region.oy; j < (region.oy + region.ny); j++) {
                if (image->region.oy <= j && j < (image->region.oy + image->region.ny)) {
                    for (int i = region.ox; i < (region.ox + region.nx); i++) {
                        if (image->region.ox <= i && i < (image->region.ox + image->region.nx)) {
                            image->get(i - image->region.ox, j - image->region.oy, k - image->region.oz, v);
                        } else {
                            v = pad_v;
                        }
                        resized_image->set(i - region.ox, j - region.oy, k - region.oz, v);
                    }
                } else {
                    for (int i = 0; i < region.nx; i++) {
                        resized_image->set(i, j - region.oy, k - region.oz, pad_v);
                    }
                }
            }
        } else {
            for (int j = 0; j < region.ny; j++) {
                for (int i = 0; i < region.nx; i++) {
                    resized_image->set(i, j, k - region.oz, pad_v);
                }
            }
        }
    }
    // memorize position of image region within original image domain
    // if the resize operation is considered to be a sampling of the image
    // data only within the original image domain
    if (subdomain) {
        resized_image->hdr.dim[1] = image->hdr.dim[1];
        resized_image->hdr.dim[2] = image->hdr.dim[2];
        resized_image->hdr.dim[3] = image->hdr.dim[3];
        resized_image->region     = region;
    }
    return resized_image;
}

// ---------------------------------------------------------------------------
Image* ResampleImage(const Image* image, Fvector3d pixdim, Image::Region region)
{
    // check arguments
    assert(image != NULL);
    if (pixdim.x <= 0 || pixdim.y <= 0 || (pixdim.z <= 0 && image->region.nz > 1)) {
        BASIS_THROW(invalid_argument, "Cannot resample image on grid with invalid voxel size ["
                << pixdim.x << ", " << pixdim.y << ", " << pixdim.z << "]");
    }
    // keep image extent (approximately) unchanged if region not specified
    if (region.nx <= 0) {
        region.nx = static_cast<int>(ceil(image->region.nx * image->hdr.pixdim[1] / pixdim.x));
        region.ox = 0;
    }
    if (region.ny <= 0) {
        region.ny = static_cast<int>(ceil(image->region.ny * image->hdr.pixdim[2] / pixdim.y));
        region.oy = 0;
    }
    if (region.nz <= 0) {
        region.nz = static_cast<int>(ceil(image->region.nz * image->hdr.pixdim[3] / pixdim.z));
        region.oz = 0;
    }
    // allocate output image
    Image* resampled_image = new Image(region.nx, region.ny, region.nz,
                                       image->hdr.datatype,
                                       image->GetNumberOfComponents(),
                                       image->imgfmt);
    resampled_image->CopyTransform  (image, true, false); // qform only
    resampled_image->CopyDataScaling(image);
    resampled_image->CopyMetaData   (image);
    // adjust voxel size
    resampled_image->hdr.pixdim[1] = pixdim.x;
    resampled_image->hdr.pixdim[2] = pixdim.y;
    resampled_image->hdr.pixdim[3] = pixdim.z;
    // adjust image offset - after voxel size!
    if (resampled_image->hdr.qform_code != NIFTI_XFORM_UNKNOWN) {
        resampled_image->hdr.qoffset_x += (region.ox - image->region.ox) * image->hdr.pixdim[1];
        resampled_image->hdr.qoffset_y += (region.oy - image->region.oy) * image->hdr.pixdim[2];
        resampled_image->hdr.qoffset_z += (region.oz - image->region.oz) * image->hdr.pixdim[3];
    }
    // resample image data
    vector<float> v;
    for (int k = 0; k < region.nz; k++) {
        for (int j = 0; j < region.ny; j++) {
            for (int i = 0; i < region.nx; i++) {
                image->get(static_cast<float>((region.ox + i) * pixdim.x / image->hdr.pixdim[1]),
                           static_cast<float>((region.oy + j) * pixdim.y / image->hdr.pixdim[2]),
                           static_cast<float>((region.oz + k) * pixdim.z / image->hdr.pixdim[3]), v);
                resampled_image->set(i, j, k, v);
            }
        }
    }
    return resampled_image;
}

// ===========================================================================
// smoothing
// ===========================================================================

// ---------------------------------------------------------------------------
Image* SmoothImage(const Image* image)
{
    assert(image != NULL);
    // create smoothing kernel
    Ivector2d radius;
    Fvector2d sigma;
    GetGaussParameters(image->region.nx, radius.x, sigma.x);
    GetGaussParameters(image->region.ny, radius.y, sigma.y);
    Image* kernel = Gauss(radius, sigma, image->imgfmt);
    if (kernel == NULL) return NULL;
    // convolve image by smoothing kernel
    Image* smoothed_image = ConvolveImage(image, kernel);
    // destroy smoothing kernel
    delete kernel;
    return smoothed_image;
}

// ---------------------------------------------------------------------------
Image* DownsampleImage(const Image* image, int ratio)
{
    assert(image != NULL);
    if (ratio < 1) {
        BASIS_THROW(invalid_argument, "Invalid downsample ratio: " << ratio);
    }
    // // downsample ratio in each dimension
    // if (image->region.nx < ratio || image->region.ny < ratio ||
            // (image->region.nz > 1 && image->region.nz < ratio)) {
        // BASIS_THROW(invalid_argument, "Image size (" <<image->region.nx <<", " << image->region.ny
                // <<", "<< image->region.nz<<") smaller than downsample ratio of " << ratio << " in at least one dimension!");
    // }
    Ivector3d downratio;
    downratio.x = ratio;
    downratio.y = ratio;
    downratio.z = 1;
    if (image->region.nz > ratio) downratio.z = ratio;
    // in case of a padded image, preserve image boundary at original "left-most" voxel
    // such that downsampling the padded image followed by cropping results in the same
    // image grid as downsampling the unpadded image directly. therefore, skip the first
    // offset voxels of the padded image.
    Ivector3d offset;
    if (image->region.ox < 0) offset.x = (- image->region.ox) % downratio.x;
    else                      offset.x = 0;
    if (image->region.oy < 0) offset.y = (- image->region.oy) % downratio.y;
    else                      offset.y = 0;
    if (image->region.oz < 0) offset.z = (- image->region.oz) % downratio.z;
    else                      offset.z = 0;
    // size of downsample image
    Ivector3d downsampledSize;
    downsampledSize.x = static_cast<int>(image->region.nx / downratio.x);
    downsampledSize.y = static_cast<int>(image->region.ny / downratio.y);
    downsampledSize.z = static_cast<int>(image->region.nz / downratio.z);
    // allocate output image
    Image* downsampledImage = new Image(downsampledSize.x,
                                        downsampledSize.y,
                                        downsampledSize.z,
                                        image->hdr.datatype,
                                        image->hdr.dim[5],
                                        image->imgfmt);
    if (!downsampledImage) return NULL;
    downsampledImage->CopyTransform  (image, true, false); // qform only
    downsampledImage->CopyDataScaling(image);
    downsampledImage->CopyMetaData   (image);
    // copy region - i.e., information about original image size before a potential ResizeImage()
    downsampledImage->hdr.dim[1] = static_cast<int>(image->hdr.dim[1] / downratio.x);
    downsampledImage->hdr.dim[2] = static_cast<int>(image->hdr.dim[2] / downratio.y);
    downsampledImage->hdr.dim[3] = static_cast<int>(image->hdr.dim[3] / downratio.z);
    downsampledImage->region.ox  = static_cast<int>(copysignf((fabs(image->region.ox + offset.x) + downratio.x - 1) / downratio.x, image->region.ox));
    downsampledImage->region.oy  = static_cast<int>(copysignf((fabs(image->region.oy + offset.y) + downratio.y - 1) / downratio.y, image->region.oy));
    downsampledImage->region.oz  = static_cast<int>(copysignf((fabs(image->region.oz + offset.z) + downratio.z - 1) / downratio.z, image->region.oz));
    // adjust voxel size
    downsampledImage->hdr.pixdim[1] = image->hdr.pixdim[1] * downratio.x;
    downsampledImage->hdr.pixdim[2] = image->hdr.pixdim[2] * downratio.y;
    downsampledImage->hdr.pixdim[3] = image->hdr.pixdim[3] * downratio.z;
    // adjust image offset - after voxel size!
    if (downsampledImage->hdr.qform_code != NIFTI_XFORM_UNKNOWN) {
        downsampledImage->hdr.qoffset_x += 0.5 * (downsampledImage->hdr.pixdim[1] - image->hdr.pixdim[1]) + offset.x * image->hdr.pixdim[1];
        downsampledImage->hdr.qoffset_y += 0.5 * (downsampledImage->hdr.pixdim[2] - image->hdr.pixdim[2]) + offset.y * image->hdr.pixdim[2];
        downsampledImage->hdr.qoffset_z += 0.5 * (downsampledImage->hdr.pixdim[3] - image->hdr.pixdim[3]) + offset.z * image->hdr.pixdim[3];
    }
    // downsample image
    vector<float> v;
    for (int k = 0; k < downsampledSize.z; k++) {
        for (int j = 0; j < downsampledSize.y; j++) {
            for (int i = 0; i < downsampledSize.x; i++) {
                image->get(offset.x + i * downratio.x, offset.y + j * downratio.y, offset.z + k * downratio.z, v);
                // scale displacements of 2D or 3D DRAMMS deformation field given
                // in voxel units accordingly as the units' scale has changed
                if (image->imgfmt == Image::FORMAT_DRAMMS) {
                    if (v.size() == 2 || v.size() == 3) {
                        v[0] /= downratio.x;
                        v[1] /= downratio.y;
                    }
                    if (v.size() == 3) v[2] /= downratio.z;
                }
                downsampledImage->set(i, j, k, v);
            }
        }
    }
    return downsampledImage;
}

// ---------------------------------------------------------------------------
Image* SmoothAndDownsampleImage(const Image* image, int ratio)
{
    assert(image != NULL);
    // create smoothing kernel
    Ivector2d radius;
    Fvector2d sigma;
    GetGaussParameters(image->region.nx, radius.x, sigma.x);
    GetGaussParameters(image->region.ny, radius.y, sigma.y);
    if (ratio >= 4) {
        sigma.x = 1.5 * sigma.x;
        sigma.y = 1.5 * sigma.y;
    }
    Image* kernel = Gauss(radius, sigma, image->imgfmt);
    if (kernel == NULL) return NULL;
    // convolve image by smoothing kernel
    Image* smoothed_image = new Image(image->region.nx,
                                      image->region.ny,
                                      image->region.nz,
                                      image->hdr.datatype,
                                      image->GetNumberOfComponents(),
                                      image->imgfmt);
    if (smoothed_image == NULL) {
        delete kernel;
        return NULL;
    }
    smoothed_image->CopyRegion     (image);
    smoothed_image->CopyTransform  (image);
    smoothed_image->CopyDataScaling(image);
    smoothed_image->CopyMetaData   (image);
    ConvolveImage(image, kernel, smoothed_image);
    delete kernel;
    // downsample image
    Image* downsampled_image = NULL;
    try {
        downsampled_image = DownsampleImage(smoothed_image, ratio);
    } catch (invalid_argument& e) {
        delete smoothed_image;
        throw e;
    }
    delete smoothed_image;
    return downsampled_image;
}

// ===========================================================================
// operations on transformations
// ===========================================================================

// ---------------------------------------------------------------------------
// inlined into functions following below for better running time
inline Fvector3d TransformPoint(const Image::Transform m, const Fvector3d& p)
{
    Fvector3d x;
    x.x = p.x * m.m[0][0] + p.y * m.m[0][1] + p.z * m.m[0][2] + m.m[0][3];
    x.y = p.x * m.m[1][0] + p.y * m.m[1][1] + p.z * m.m[1][2] + m.m[1][3];
    x.z = p.x * m.m[2][0] + p.y * m.m[2][1] + p.z * m.m[2][2] + m.m[2][3];
    return x;
}

// ---------------------------------------------------------------------------
Image::Transform IdentityTransform()
{
    Image::Transform T;
    memset(T.m, 0, 16 * sizeof(float));
    T.m[0][0] = 1;
    T.m[1][1] = 1;
    T.m[2][2] = 1;
    T.m[3][3] = 1;
    return T;
}

// ---------------------------------------------------------------------------
Image::Transform InvertTransform(const Image::Transform& T)
{
    // use of Mat_Inverse() here as it uses double type for its elements
    // and reproduces the exact same inverse matrix as was used by the
    // initial DRAMMS implementation. otherwise, the nifti_mat44_inverse()
    // function could be used.
    Image::Transform T_inv;
    Matrix* a;
    Matrix* b;
    CreateMatrix(&a, 4, 4);
    CreateMatrix(&b, 4, 4);
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            a->data[r][c] = static_cast<double>(T.m[r][c]);
        }
    }
    Mat_Inverse(a, b);
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            T_inv.m[r][c] = static_cast<float>(b->data[r][c]);
        }
    }
    FreeMatrix(a);
    FreeMatrix(b);
    return T_inv;
}

// ---------------------------------------------------------------------------
Image* InvertTransform(const Image* D, int num_samples, int boundary)
{
    assert(D != NULL);

    // -----------------------------------------------------------------------
    // constants and local variables

    // size of reverse deformation field - will be increased below
    int nx = D->region.nx;
    int ny = D->region.ny;
    int nz = D->region.nz;

    // number of vector components, i.e., dimension of deformation field
    const int N = D->GetNumberOfComponents();

    // spacing of subvoxel grid in voxel units of input deformation field
    float dx = 1.0 / (2.0 * num_samples + 1.0);
    float dy = 1.0 / (2.0 * num_samples + 1.0);
    float dz = 1.0 / (2.0 * num_samples + 1.0);

    float ii, jj, kk;    // sampled subvoxel position
    int   ll, mm, nn;    // voxel position in output deformation field
    int   LL, MM, NN;    // l + 1, m + 1, n + 1
    float b, c, d;       // coefficients for linear distribution of displacement
                         // at sub-voxel position in output deformation field
    float b1, c1, d1;    // 1 - a, 1 - b, 1 - c
    float a;             // combined coefficients for weight normalization
    float w[8];          // weights of displacement at neighboring voxel position
    Fvector3d v;         // interpolated displacement

    // -----------------------------------------------------------------------
    // check input arguments
    if (N == 2) {
        // TODO
        BASIS_THROW(runtime_error, "Estimation of inverse deformation of 2D deformation field not implemented!");
    }

    if (D->hdr.datatype != DT_FLOAT || nx <= 0 || ny <= 0 || nz <= 0 || N < 2 || N > 3 || (nz == 1 && N != 2)) {
        BASIS_THROW(invalid_argument, "Invalid image given for estimation of inverse deformation field!");
    }
    if (D->imgfmt != Image::FORMAT_DRAMMS) {
        BASIS_THROW(invalid_argument, "Deformation field must be in DRAMMS format!");
    }

    if (num_samples < 0) num_samples = 0;

    // -----------------------------------------------------------------------
    // "extend" input deformation field - extrapolation done on-the-fly by
    //                                    interpolate function Image::get()

    // Note: The reason for doing so is, that this yields a smoother reverse
    //       deformation field at the boundary as displacments extrapolated
    //       from the input deformation field can have an impact to the
    //       reverse deformation field inside the original image domain.
    int Nx = boundary + nx + boundary;
    int Ny = boundary + ny + boundary;
    int Nz = boundary + nz + boundary;

    // -----------------------------------------------------------------------
    // allocate output deformation field
    Image* D_inv = new Image(Nx, Ny, Nz, D->hdr.datatype, N, D->imgfmt);
    if (D_inv == NULL) return NULL;
    D_inv->CopyTransform(D);
    D_inv->CopyDataScaling(D);
    D_inv->CopyMetaData(D);

    // -----------------------------------------------------------------------
    // allocate weight image required for normalization
    Image* weights = new Image(Nx, Ny, Nz, DT_FLOAT, 1, D->imgfmt);
    if (weights == NULL) {
        delete D_inv;
        return NULL;
    }

    // -----------------------------------------------------------------------
    // estimate inverse deformation field
    if (N == 2) {
        // TODO
    } else {
        float***     pw = weights->img.fl;
        Fvector3d*** pd = D      ->img.v3;
        Fvector3d*** pr = D_inv  ->img.v3;

        for (int k = 0; k < Nz; k++) {
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    // for each subvoxel position on the grid of size
                    // (2 * num_samples + 1) centered at voxel (i, j, k)...
                    //
                    // Note that this is identical to first upsampling the input
                    // deformation field and then performing the following computations.
                    for (int dk = - num_samples; dk <= num_samples; dk++) {
                        for (int di = - num_samples; di <= num_samples; di++) {
                            for (int dj = - num_samples; dj <= num_samples; dj++) {
                                // current subvoxel position
                                ii = static_cast<float>(i) + static_cast<float>(di) * dx;
                                jj = static_cast<float>(j) + static_cast<float>(dj) * dy;
                                kk = static_cast<float>(k) + static_cast<float>(dk) * dz;
                                // interpolate displacement
                                ll = static_cast<int>(floor(ii)) - boundary;
                                mm = static_cast<int>(floor(jj)) - boundary;
                                nn = static_cast<int>(floor(kk)) - boundary;

                                b = ii - (ll + boundary); b1 = 1.0 - b;
                                c = jj - (mm + boundary); c1 = 1.0 - c;
                                d = kk - (nn + boundary); d1 = 1.0 - d;

                                a = (d1 * ((b1 * c1) + (b * c1) + (b1 * c) + (b * c))
                                        + d * ((b1 * c1) + (b * c1) + (b1 * c) + (b * c)));

                                w[0] = d1 * (b1 * c1) / a; // (ll, mm, nn)
                                w[1] = d1 * (b  * c1) / a; // (LL, mm, nn)
                                w[2] = d1 * (b1 * c ) / a; // (ll, MM, nn)
                                w[3] = d1 * (b  * c ) / a; // (LL, MM, nn)
                                w[4] = d  * (b1 * c1) / a; // (ll, mm, NN)
                                w[5] = d  * (b  * c1) / a; // (LL, mm, NN)
                                w[6] = d  * (b1 * c ) / a; // (ll, MM, NN)
                                w[7] = d  * (b  * c ) / a; // (LL, MM, NN)

                                LL = ll + 1;
                                MM = mm + 1;
                                NN = nn + 1;

                                if      (ll < 0  ) ll = 0;
                                else if (ll >= nx) ll = nx - 1;
                                if      (LL < 0  ) LL = 0;
                                else if (LL >= nx) LL = nx - 1;
                                if      (mm < 0  ) mm = 0;
                                else if (mm >= ny) mm = ny - 1;
                                if      (MM < 0  ) MM = 0;
                                else if (MM >= ny) MM = ny - 1;
                                if      (nn < 0  ) nn = 0;
                                else if (nn >= nz) nn = nz - 1;
                                if      (NN < 0  ) NN = 0;
                                else if (NN >= nz) NN = nz - 1;

                                v.x = w[0] * pd[nn][ll][mm].x +
                                      w[1] * pd[nn][LL][mm].x +
                                      w[2] * pd[nn][ll][MM].x +
                                      w[3] * pd[nn][LL][MM].x +
                                      w[4] * pd[NN][ll][mm].x +
                                      w[5] * pd[NN][LL][mm].x +
                                      w[6] * pd[NN][ll][MM].x +
                                      w[7] * pd[NN][LL][MM].x;

                                v.y = w[0] * pd[nn][ll][mm].y +
                                      w[1] * pd[nn][LL][mm].y +
                                      w[2] * pd[nn][ll][MM].y +
                                      w[3] * pd[nn][LL][MM].y +
                                      w[4] * pd[NN][ll][mm].y +
                                      w[5] * pd[NN][LL][mm].y +
                                      w[6] * pd[NN][ll][MM].y +
                                      w[7] * pd[NN][LL][MM].y;

                                v.z = w[0] * pd[nn][ll][mm].z +
                                      w[1] * pd[nn][LL][mm].z +
                                      w[2] * pd[nn][ll][MM].z +
                                      w[3] * pd[nn][LL][MM].z +
                                      w[4] * pd[NN][ll][mm].z +
                                      w[5] * pd[NN][LL][mm].z +
                                      w[6] * pd[NN][ll][MM].z +
                                      w[7] * pd[NN][LL][MM].z;
                                // add displacement to subvoxel position
                                ii += v.x;
                                jj += v.y;
                                kk += v.z;
                                ll = static_cast<int>(floor(ii));
                                mm = static_cast<int>(floor(jj));
                                nn = static_cast<int>(floor(kk));
                                // and distribute it to surrounding voxels of inverse field
                                if (0 <= ll && ll < Nx - 1 &&
                                        0 <= mm && mm < Ny - 1 &&
                                        0 <= nn && nn < Nz - 1) {
                                    b = ii - ll; b1 = 1.0 - b;
                                    c = jj - mm; c1 = 1.0 - c;
                                    d = kk - nn; d1 = 1.0 - d;

                                    a = (d1 * ((b1 * c1) + (b * c1) + (b1 * c) + (b * c))
                                            + d * ((b1 * c1) + (b * c1) + (b1 * c) + (b * c)));

                                    w[0] = d1 * (b1 * c1) / a; // (ll, mm, nn)
                                    w[1] = d1 * (b  * c1) / a; // (LL, mm, nn)
                                    w[2] = d1 * (b1 * c ) / a; // (ll, MM, nn)
                                    w[3] = d1 * (b  * c ) / a; // (LL, MM, nn)
                                    w[4] = d  * (b1 * c1) / a; // (ll, mm, NN)
                                    w[5] = d  * (b  * c1) / a; // (LL, mm, NN)
                                    w[6] = d  * (b1 * c ) / a; // (ll, MM, NN)
                                    w[7] = d  * (b  * c ) / a; // (LL, MM, NN)

                                    LL = ll + 1;
                                    MM = mm + 1;
                                    NN = nn + 1;

                                    // update weights
                                    pw[nn][ll][mm] += w[0];
                                    pw[nn][LL][mm] += w[1];
                                    pw[nn][ll][MM] += w[2];
                                    pw[nn][LL][MM] += w[3];
                                    pw[NN][ll][mm] += w[4];
                                    pw[NN][LL][mm] += w[5];
                                    pw[NN][ll][MM] += w[6];
                                    pw[NN][LL][MM] += w[7];

                                    // update vector components
                                    pr[nn][ll][mm].x += w[0] * v.x;
                                    pr[nn][LL][mm].x += w[1] * v.x;
                                    pr[nn][ll][MM].x += w[2] * v.x;
                                    pr[nn][LL][MM].x += w[3] * v.x;
                                    pr[NN][ll][mm].x += w[4] * v.x;
                                    pr[NN][LL][mm].x += w[5] * v.x;
                                    pr[NN][ll][MM].x += w[6] * v.x;
                                    pr[NN][LL][MM].x += w[7] * v.x;

                                    pr[nn][ll][mm].y += w[0] * v.y;
                                    pr[nn][LL][mm].y += w[1] * v.y;
                                    pr[nn][ll][MM].y += w[2] * v.y;
                                    pr[nn][LL][MM].y += w[3] * v.y;
                                    pr[NN][ll][mm].y += w[4] * v.y;
                                    pr[NN][LL][mm].y += w[5] * v.y;
                                    pr[NN][ll][MM].y += w[6] * v.y;
                                    pr[NN][LL][MM].y += w[7] * v.y;

                                    pr[nn][ll][mm].z += w[0] * v.z;
                                    pr[nn][LL][mm].z += w[1] * v.z;
                                    pr[nn][ll][MM].z += w[2] * v.z;
                                    pr[nn][LL][MM].z += w[3] * v.z;
                                    pr[NN][ll][mm].z += w[4] * v.z;
                                    pr[NN][LL][mm].z += w[5] * v.z;
                                    pr[NN][ll][MM].z += w[6] * v.z;
                                    pr[NN][LL][MM].z += w[7] * v.z;
                                }
                            }
                        }
                    }
                }
            }
        }
        // normalize and reverse displacements
        for (int k = 0; k < Nz; k++) {
            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < Ny; j++) {
                    const float& w = pw[k][i][j];
                    if (w > 0) {
                        pr[k][i][j].x /= -w;
                        pr[k][i][j].y /= -w;
                        pr[k][i][j].z /= -w;
                    } else {
                        pr[k][i][j].x = 0.0f;
                        pr[k][i][j].y = 0.0f;
                        pr[k][i][j].z = 0.0f;
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // clean up
    delete weights;

    // -----------------------------------------------------------------------
    // remove additional boundary again
    if (boundary > 0) {
        Image* tmp = new Image(D->region.nx,
                               D->region.ny,
                               D->region.nz,
                               D->hdr.datatype,
                               D->hdr.dim[5],
                               D->imgfmt);
        if (tmp == NULL) {
            delete D_inv;
            return NULL;
        }
        tmp->CopyTransform(D);
        tmp->CopyDataScaling(D_inv);
        tmp->CopyMetaData(D_inv);

        if (N == 2) {
            // TODO
        } else {
            Fvector3d*** pr = D_inv->img.v3;
            Fvector3d*** pt = tmp->img.v3;
            for (int k = boundary; k < Nz - boundary; k++) {
                for (int i = boundary; i < Nx - boundary; i++) {
                    for (int j = boundary; j < Ny - boundary; j++) {
                        pt[k - boundary][i - boundary][j - boundary] = pr[k][i][j];
                    }
                }
            }
        }
        delete D_inv;
        D_inv = tmp;
    }

    return D_inv;
}

// ---------------------------------------------------------------------------
Image::Transform ConcatenateTransforms(const Image::Transform& T1, // A -> B
                                       const Image::Transform& T2) // B -> C
{
    Image::Transform T;
    for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
            T.m[r][c] =   T2.m[r][0] * T1.m[0][c]
                        + T2.m[r][1] * T1.m[1][c]
                        + T2.m[r][2] * T1.m[2][c]
                        + T2.m[r][3] * T1.m[3][c];
        }
    }
    return T;
}

// ---------------------------------------------------------------------------
Image* ConcatenateTransforms(const Image::Header&    A,
                             const Image::Transform& T, // A -> B
                             const Image*            D) // B -> B
{
    return ConcatenateTransforms(A, D->hdr, T, D);
}

// ---------------------------------------------------------------------------
Image* ConcatenateTransforms(const Image::Header&    A,
							 const Image::Header&    B,
                             const Image::Transform& T, // A -> B
                             const Image*            D) // B -> C
{
    assert(D != NULL);
    // check input arguments
    if (D->imgfmt != Image::FORMAT_DRAMMS) {
        BASIS_THROW(invalid_argument, "Deformation field must be given in DRAMMS representation!");
    }
    // size of deformation field
    Ivector3d size;
    size.x = D->region.nx;
    size.y = D->region.ny;
    size.z = D->region.nz;
    // allocate output deformation field
    Image* D_out = new Image(size.x, size.y, size.z,
                             D->hdr.datatype,
                             D->GetNumberOfComponents(),
                             D->imgfmt);
    if (D_out == NULL) return NULL;
    D_out->CopyTransform(D);
    D_out->CopyMetaData (D);
    // prepare linear transformation required for composition
    //
    // NOTE: The affine transformation obtained by FLIRT is in world coordinates,
    //       with a conversion based on a so-called sampling matrix. This
    //       sampling matrix has the voxel sizes on the diagonal, but does not
    //       encode any image orientation. Thus, the world coordinate of voxel
    //       (i,j,k) is simply (i*dx,j*dy,k*dz) regardless of the orientation
    //       matrices (qform/sform) of the image. The left-right order of the first
    //       image dimension is further always in radiological order. Therefore,
    //       FSL programs such as FLIRT swap the order of the first dimension
    //       of an image while reading the image data if it is stored in
    //       neurological order (positive determinant of qform/sfrom).
    Image::Transform vox2mm = IdentityTransform();
    vox2mm.m[0][0] = B.pixdim[1];
    vox2mm.m[1][1] = B.pixdim[2];
    vox2mm.m[2][2] = B.pixdim[3];
    if (GetLeftRightOrder(B, D->imgfmt) == Image::LR_NEUROLOGICAL) {
        Image::Transform swap = IdentityTransform();
        if (D->imgfmt == Image::FORMAT_DRAMMS) {
            swap.m[1][1] = -1;
            swap.m[1][3] = B.dim[2] - 1;
        } else {
            swap.m[0][0] = -1;
            swap.m[0][3] = B.dim[1] - 1;
        }
        vox2mm = ConcatenateTransforms(swap, vox2mm);
    }
    Image::Transform mm2vox = IdentityTransform();
    mm2vox.m[0][0] = 1.0f / A.pixdim[1];
    mm2vox.m[1][1] = 1.0f / A.pixdim[2];
    mm2vox.m[2][2] = 1.0f / A.pixdim[3];
    if (GetLeftRightOrder(A, D->imgfmt) == Image::LR_NEUROLOGICAL) {
        Image::Transform swap = IdentityTransform();
        if (D->imgfmt == Image::FORMAT_DRAMMS) {
            swap.m[1][1] = -1;
            swap.m[1][3] = A.dim[2] - 1;
        } else {
            swap.m[0][0] = -1;
            swap.m[0][3] = A.dim[1] - 1;
        }
        mm2vox = ConcatenateTransforms(mm2vox, swap);
    }
    Image::Transform M = InvertTransform(T); // inverse affine transformation
    M = ConcatenateTransforms(vox2mm, M);    // voxel in C to radiological flirt mm
    M = ConcatenateTransforms(M, mm2vox);    // radiological flirt mm to voxel in A
    // compute unified deformation field
    Fvector3d     x, y;
    vector<float> v;
    for (int k = 0; k < size.z; k++) {
        for (int i = 0; i < size.x; i++) {
            for (int j = 0; j < size.y; j++) {
                // get displacement of voxel
                D->get(i, j, k, v);
                if (v.size() < 3) v.resize(3, 0);
                // get voxel coordinates in image C
                x.x = static_cast<float>(i);
                x.y = static_cast<float>(j);
                x.z = static_cast<float>(k);
                // apply displacement to get voxel in image B
                y.x = x.x + v[0];
                y.y = x.y + v[1];
                y.z = x.z + v[2];
                // apply linear transformation to get voxel in image A
                y = TransformPoint(M, y);
                // set output displacement
                v[0] = y.x - x.x;
                v[1] = y.y - x.y;
                if (D_out->GetNumberOfComponents() == 2) {
                    v.resize(2);
                } else {
                    v[2] = y.z - x.z;
                }
                D_out->set(i, j, k, v);
            }
        }
    }
    // return resulting deformation field
    return D_out;
}

// ---------------------------------------------------------------------------
Image* ConcatenateTransforms(const Image* D1, const Image* D2)
{
    const int nx1 = D1->region.nx;
    const int ny1 = D1->region.ny;
    const int nz1 = D1->region.nz;
    const int nx2 = D2->region.nx;
    const int ny2 = D2->region.ny;
    const int nz2 = D2->region.nz;
    const int nu  = D1->hdr.dim[5];

    // check input arguments
    if (D1->imgfmt != Image::FORMAT_DRAMMS) {
        BASIS_THROW(invalid_argument, "Deformation field must be given in DRAMMS representation!");
    }
    if (D2->hdr.dim[5] != nu) {
        BASIS_THROW(invalid_argument, "The dimensions of the displacement vectors do not match!");
    }

    Image* D = new Image(nx2, ny2, nz2, D2->hdr.datatype, nu, D2->imgfmt);
    if (D == NULL) return NULL;
    D->CopyTransform  (D2);
    D->CopyDataScaling(D2);
    D->CopyMetaData   (D2);

    vector<float> v2, v1CCC, v1CCF, v1CFC, v1CFF, v1FCC, v1FCF, v1FFC, v1FFF;
    float ii, jj, kk;
    int iFloor, jFloor, kFloor;
    int iCeil,  jCeil,  kCeil;

    for (int k = 0; k < nz2; k++) {
        for (int j = 0; j < ny2; j++) {
            for (int i = 0; i < nx2; i++) {
                // get input vector of second field
                D2->get(i, j, k, v2);

                // get intermediate position
                ii = i + v2[0];
                jj = j + v2[1];
                kk = k + v2[2];

                // calculate interpolation coefficients
                iFloor = static_cast<int>(fmax(fmin(floor(ii), static_cast<float>(nx1 - 1) ), 0.0f));
                iCeil  = static_cast<int>(fmax(fmin(ceil (ii), static_cast<float>(nx1 - 1) ), 0.0f));
                jFloor = static_cast<int>(fmax(fmin(floor(jj), static_cast<float>(ny1 - 1) ), 0.0f));
                jCeil  = static_cast<int>(fmax(fmin(ceil (jj), static_cast<float>(ny1 - 1) ), 0.0f));
                kFloor = static_cast<int>(fmax(fmin(floor(kk), static_cast<float>(nz1 - 1) ), 0.0f));
                kCeil  = static_cast<int>(fmax(fmin(ceil (kk), static_cast<float>(nz1 - 1) ), 0.0f));

                if ((iFloor == iCeil) && (iFloor != (nx1 - 1))) iCeil  += 1;
                if ((iFloor == iCeil) && (iFloor == (nx1 - 1))) iFloor -= 1;
                if ((jFloor == jCeil) && (jFloor != (ny1 - 1))) jCeil  += 1;
                if ((jFloor == jCeil) && (jFloor == (ny1 - 1))) jFloor -= 1;
                if ((kFloor == kCeil) && (kFloor != (nz1 - 1))) kCeil  += 1;
                if ((kFloor == kCeil) && (kFloor == (nz1 - 1))) kFloor -= 1;

                // get input vectors of first field to interpolate
                D1->get(iCeil,  jCeil,  kCeil,  v1CCC);
                D1->get(iCeil,  jCeil,  kFloor, v1CCF);
                D1->get(iCeil,  jFloor, kCeil,  v1CFC);
                D1->get(iCeil,  jFloor, kFloor, v1CFF);
                D1->get(iFloor, jCeil,  kCeil,  v1FCC);
                D1->get(iFloor, jCeil,  kFloor, v1FCF);
                D1->get(iFloor, jFloor, kCeil,  v1FFC);
                D1->get(iFloor, jFloor, kFloor, v1FFF);

                // add incremental displacement
                for (int n = 0; n < nu; n++) {
                    v2[n] += (ii - iFloor) * (jj - jFloor) * (kk - kFloor) * v1CCC[n]
                           + (ii - iFloor) * (jj - jFloor) * (kCeil  - kk) * v1CCF[n]
                           + (ii - iFloor) * (jCeil  - jj) * (kk - kFloor) * v1CFC[n]
                           + (ii - iFloor) * (jCeil  - jj) * (kCeil  - kk) * v1CFF[n]
                           + (iCeil  - ii) * (jj - jFloor) * (kk - kFloor) * v1FCC[n]
                           + (iCeil  - ii) * (jj - jFloor) * (kCeil  - kk) * v1FCF[n]
                           + (iCeil  - ii) * (jCeil  - jj) * (kk - kFloor) * v1FFC[n]
                           + (iCeil  - ii) * (jCeil  - jj) * (kCeil  - kk) * v1FFF[n];
                }

                // set output vector
                D->set(i, j, k, v2);
            }
        }
    }

    return D;
}

// ---------------------------------------------------------------------------
Image* SubtractTransforms(const Image::Header&    A,
                          const Image*            D, // A -> B unified deformation
                          const Image::Transform& T) // A -> B affine part
{
    assert(D != NULL);
    // check input arguments
    if (D->imgfmt != Image::FORMAT_DRAMMS) {
        BASIS_THROW(invalid_argument, "Deformation field must be given in DRAMMS representation!");
    }
    // size of deformation field
    Ivector3d size;
    size.x = D->region.nx;
    size.y = D->region.ny;
    size.z = D->region.nz;
    // allocate output deformation field
    Image* D_out = new Image(size.x, size.y, size.z,
                             D->hdr.datatype,
                             D->GetNumberOfComponents(),
                             D->imgfmt);
    if (D_out == NULL) return NULL;
    D_out->CopyTransform(D);
    D_out->CopyMetaData (D);
    // prepare linear transformation required for composition
    Image::Transform vox2mm = IdentityTransform();
    vox2mm.m[0][0] = A.pixdim[1];
    vox2mm.m[1][1] = A.pixdim[2];
    vox2mm.m[2][2] = A.pixdim[3];
    if (GetLeftRightOrder(A, D->imgfmt) == Image::LR_NEUROLOGICAL) {
        Image::Transform swap = IdentityTransform();
        if (D->imgfmt == Image::FORMAT_DRAMMS) {
            swap.m[1][1] = -1;
            swap.m[1][3] = A.dim[2] - 1;
        } else {
            swap.m[0][0] = -1;
            swap.m[0][3] = A.dim[1] - 1;
        }
        vox2mm = ConcatenateTransforms(vox2mm, swap);
    }
    Image::Transform mm2vox = IdentityTransform();
    mm2vox.m[0][0] = 1.0f / D->hdr.pixdim[1];
    mm2vox.m[1][1] = 1.0f / D->hdr.pixdim[2];
    mm2vox.m[2][2] = 1.0f / D->hdr.pixdim[3];
    if (GetLeftRightOrder(D->hdr, D->imgfmt) == Image::LR_NEUROLOGICAL) {
        Image::Transform swap = IdentityTransform();
        if (D->imgfmt == Image::FORMAT_DRAMMS) {
            swap.m[1][1] = -1;
            swap.m[1][3] = D->region.ny - 1;
        } else {
            swap.m[0][0] = -1;
            swap.m[0][3] = D->region.nx - 1;
        }
        mm2vox = ConcatenateTransforms(swap, mm2vox);
    }
    Image::Transform M;                   // affine transformation from A -> B
    M = ConcatenateTransforms(vox2mm, T); // voxel in A to radiological flirt mm
    M = ConcatenateTransforms(M, mm2vox); // radiological flirt mm to voxel in B
    // compute unified deformation field
    Fvector3d     x, y;
    vector<float> v;
    for (int k = 0; k < size.z; k++) {
        for (int i = 0; i < size.x; i++) {
            for (int j = 0; j < size.y; j++) {
                // get voxel coordinates in image C
                y.x = static_cast<float>(i);
                y.y = static_cast<float>(j);
                y.z = static_cast<float>(k);
                // apply displacement to get voxel in image A
                D->get(i, j, k, v);
                if (v.size() < 3) v.resize(3, 0);
                x.x = y.x + v[0];
                x.y = y.y + v[1];
                x.z = y.z + v[2];
                // apply linear transformation to get voxel in image B
                x = TransformPoint(M, x);
                // set output displacement for deformation from B -> C
                v[0] = x.x - y.x;
                v[1] = x.y - y.y;
                if (D_out->GetNumberOfComponents() == 2) {
                    v.resize(2);
                } else {
                    v[2] = x.z - y.z;
                }
                D_out->set(i, j, k, v);
            }
        }
    }
    // return resulting deformation field
    return D_out;
}

// ---------------------------------------------------------------------------
Image* SubtractTransforms(const Image::Header&    A,
						  const Image::Header&    B,
                          const Image*            D, // A -> C
                          const Image::Transform& T) // A -> B
{
    assert(D != NULL);
    // check input arguments
    if (D->imgfmt != Image::FORMAT_DRAMMS) {
        BASIS_THROW(invalid_argument, "Deformation field must be given in DRAMMS representation!");
    }
    // size of deformation field
    Ivector3d size;
    size.x = D->region.nx;
    size.y = D->region.ny;
    size.z = D->region.nz;
    // allocate output deformation field
    Image* D_out = new Image(size.x, size.y, size.z,
                             D->hdr.datatype,
                             D->GetNumberOfComponents(),
                             D->imgfmt);
    if (D_out == NULL) return NULL;
    D_out->CopyTransform(D);
    D_out->CopyMetaData (D);
    // prepare linear transformation required for composition
    Image::Transform vox2mm = IdentityTransform();
    vox2mm.m[0][0] = A.pixdim[1];
    vox2mm.m[1][1] = A.pixdim[2];
    vox2mm.m[2][2] = A.pixdim[3];
    if (GetLeftRightOrder(A, D->imgfmt) == Image::LR_NEUROLOGICAL) {
        Image::Transform swap = IdentityTransform();
        if (D->imgfmt == Image::FORMAT_DRAMMS) {
            swap.m[1][1] = -1;
            swap.m[1][3] = A.dim[2] - 1;
        } else {
            swap.m[0][0] = -1;
            swap.m[0][3] = A.dim[1] - 1;
        }
        vox2mm = ConcatenateTransforms(vox2mm, swap);
    }
    Image::Transform mm2vox = IdentityTransform();
    mm2vox.m[0][0] = 1.0f / B.pixdim[1];
    mm2vox.m[1][1] = 1.0f / B.pixdim[2];
    mm2vox.m[2][2] = 1.0f / B.pixdim[3];
    if (GetLeftRightOrder(B, D->imgfmt) == Image::LR_NEUROLOGICAL) {
        Image::Transform swap = IdentityTransform();
        if (D->imgfmt == Image::FORMAT_DRAMMS) {
            swap.m[1][1] = -1;
            swap.m[1][3] = B.dim[2] - 1;
        } else {
            swap.m[0][0] = -1;
            swap.m[0][3] = B.dim[1] - 1;
        }
        mm2vox = ConcatenateTransforms(swap, mm2vox);
    }
    Image::Transform M;                   // affine transformation from A -> B
    M = ConcatenateTransforms(vox2mm, T); // voxel in A to radiological flirt mm
    M = ConcatenateTransforms(M, mm2vox); // radiological flirt mm to voxel in B
    // compute unified deformation field
    Fvector3d     x, y;
    vector<float> v;
    for (int k = 0; k < size.z; k++) {
        for (int i = 0; i < size.x; i++) {
            for (int j = 0; j < size.y; j++) {
                // get voxel coordinates in image C
                y.x = static_cast<float>(i);
                y.y = static_cast<float>(j);
                y.z = static_cast<float>(k);
                // apply displacement to get voxel in image A
                D->get(i, j, k, v);
                if (v.size() < 3) v.resize(3, 0);
                x.x = y.x + v[0];
                x.y = y.y + v[1];
                x.z = y.z + v[2];
                // apply linear transformation to get voxel in image B
                x = TransformPoint(M, x);
                // set output displacement for deformation from B -> C
                v[0] = x.x - y.x;
                v[1] = x.y - y.y;
                if (D_out->GetNumberOfComponents() == 2) {
                    v.resize(2);
                } else {
                    v[2] = x.z - y.z;
                }
                D_out->set(i, j, k, v);
            }
        }
    }
    // return resulting deformation field
    return D_out;
}

// ---------------------------------------------------------------------------
Image* SubtractTransforms(const Image* D1, // A -> C
                          const Image* D2) // A -> B
{
    assert(D1 != NULL);
    assert(D2 != NULL);
    const int N = D1->GetNumberOfComponents();
    if (D2->GetNumberOfComponents() != N || !D2->HasSameSizeAs(D1)) {
        BASIS_THROW(invalid_argument, "Input deformation fields must have identical dimensions!");
    }
    Image* D_out = new Image(D1->region.nx, D1->region.ny, D1->region.nz, D1->hdr.datatype, N, D1->imgfmt);
    if (D_out == NULL) return NULL;
    D_out->CopyTransform(D1);
    D_out->CopyMetaData (D1);
    vector<float> v1(N), v2(N), vo(N);
    for (int k = 0; k < D1->region.nz; k++) {
        for (int j = 0; j < D1->region.ny; j++) {
            for (int i = 0; i < D1->region.nx; i++) {
                D1->get(i, j, k, v1);
                D2->get(i, j, k, v2);
                for (int n = 0; n < N; n++) vo[n] = v1[n] - v2[n];
                D_out->set(i, j, k, vo);
            }
        }
    }
    return D_out;
}

// ---------------------------------------------------------------------------
Image* ApplyTransform(const Image* image, const Image::Transform& T, const Image* reference, bool interpolate)
{
    // check input images
    if (image->hdr.dim[5] > 1) return NULL;
    if (image->imgfmt != reference->imgfmt) {
        BASIS_THROW(invalid_argument, "Subject and template image must be in the same format!");
    }
    // transformed image should be in the template space,
    // with the same datatype of input image
    Image* output_image = new Image(reference->hdr.dim[1],
                                    reference->hdr.dim[2],
                                    reference->hdr.dim[3],
                                    image->hdr.datatype, 1,
                                    reference->imgfmt);
    if (output_image == NULL) return NULL;
    // copy template space attributes
    output_image->CopyTransform(reference);
    // copy other attributes from input image
    output_image->CopyDataScaling(image);
    output_image->CopyMetaData(image);
    output_image->compress = image->compress;
    // prepare linear transformation required for composition
    //
    // NOTE: The affine transformation obtained by FLIRT is in world coordinates,
    //       with a conversion based on a so-called sampling matrix. This
    //       sampling matrix has the voxel sizes on the diagonal, but does not
    //       encode any image orientation. Thus, the world coordinate of voxel
    //       (i,j,k) is simply (i*dx,j*dy,k*dz) regardless of the orientation
    //       matrices (qform/sform) of the image. The left-right order of the first
    //       image dimension is further always in radiological order. Therefore,
    //       FSL programs such as FLIRT swap the order of the first dimension
    //       of an image while reading the image data if it is stored in
    //       neurological order (positive determinant of qform/sfrom).
    Image::Transform vox2mm = IdentityTransform();
    vox2mm.m[0][0] = reference->hdr.pixdim[1];
    vox2mm.m[1][1] = reference->hdr.pixdim[2];
    vox2mm.m[2][2] = reference->hdr.pixdim[3];
    if (GetLeftRightOrder(reference->hdr, reference->imgfmt) == Image::LR_NEUROLOGICAL) {
        Image::Transform swap = IdentityTransform();
        if (reference->imgfmt == Image::FORMAT_DRAMMS) {
            swap.m[1][1] = -1;
            swap.m[1][3] = reference->region.ny - 1;
        } else {
            swap.m[0][0] = -1;
            swap.m[0][3] = reference->region.nx - 1;
        }
        vox2mm = ConcatenateTransforms(swap, vox2mm);
    }
    Image::Transform mm2vox = IdentityTransform();
    mm2vox.m[0][0] = 1.0f / image->hdr.pixdim[1];
    mm2vox.m[1][1] = 1.0f / image->hdr.pixdim[2];
    mm2vox.m[2][2] = 1.0f / image->hdr.pixdim[3];
    if (GetLeftRightOrder(image->hdr, image->imgfmt) == Image::LR_NEUROLOGICAL) {
        Image::Transform swap = IdentityTransform();
        if (image->imgfmt == Image::FORMAT_DRAMMS) {
            swap.m[1][1] = -1;
            swap.m[1][3] = image->region.ny - 1;
        } else {
            swap.m[0][0] = -1;
            swap.m[0][3] = image->region.nx - 1;
        }
        mm2vox = ConcatenateTransforms(mm2vox, swap);
    }
    Image::Transform M = InvertTransform(T); // inverse linear transformation
    M = ConcatenateTransforms(vox2mm, M);    // voxel in reference to radiological flirt mm
    M = ConcatenateTransforms(M, mm2vox);    // radiological flirt mm to voxel in input image
    // interpolate input image at transformed voxels of output image
    Fvector3d x;
    float     value;
    for (int k = 0; k < output_image->region.nz; k++) {
        for (int j = 0; j < output_image->region.ny; j++) {
            for (int i = 0; i < output_image->region.nx; i++) {
                x.x = static_cast<float>(i);
                x.y = static_cast<float>(j);
                x.z = static_cast<float>(k);
                x = TransformPoint(M, x);
                value = 0.0f;
                if (interpolate) {
                    value = image->value(x.x, x.y, x.z);
                } else if (x.x >= 0 && x.x <= (image->region.nx - 1) &&
                           x.y >= 0 && x.y <= (image->region.ny - 1) &&
                           x.z >= 0 && x.z <= (image->region.nz - 1)) {
                    value = image->get(min(image->region.nx-1, static_cast<int>(x.x + 0.5f)),
                                       min(image->region.ny-1, static_cast<int>(x.y + 0.5f)),
                                       min(image->region.nz-1, static_cast<int>(x.z + 0.5f)));
                }
                output_image->set(i, j, k, value);
            }
        }
    }
    return output_image;
}

// ---------------------------------------------------------------------------
Image* ApplyTransform(const Image* image, const Image* deffield, bool interpolate)
{
    // check input images
    if (image->hdr.dim[5] > 1) {
        BASIS_THROW(invalid_argument, "Only scalar images can be warped by this program!");
    }
    if (deffield->hdr.datatype != DT_FLOAT) {
        BASIS_THROW(invalid_argument, "Data type of deformation field must be DT_FLOAT32!")
    }
    if ((image->region.nz <= 1 && deffield->hdr.dim[5] != 2) || (image->region.nz > 1 && deffield->hdr.dim[5] != 3)) {
        BASIS_THROW(invalid_argument, "Dimension of displacement vectors does not match dimension of input image!");
    }
    if (image->imgfmt != deffield->imgfmt) {
        BASIS_THROW(invalid_argument, "Subject and template image must be in the same format!");
    }
    // warped image should be in the template space,
    // with the same datatype of input image
    Image* warpedimage = new Image(deffield->region.nx,
                                   deffield->region.ny,
                                   deffield->region.nz,
                                   image->hdr.datatype, 1,
                                   deffield->imgfmt);
    if (warpedimage == NULL) return NULL;
    // copy template space attributes
    warpedimage->CopyTransform(deffield);
    // copy other attributes from input image
    warpedimage->CopyDataScaling(image);
    warpedimage->CopyMetaData(image);
    warpedimage->compress = image->compress;
    // warp image
    const int x_size_in  = image   ->region.nx;
    const int y_size_in  = image   ->region.ny;
    const int z_size_in  = image   ->region.nz;
    const int x_size_out = deffield->region.nx;
    const int y_size_out = deffield->region.ny;
    const int z_size_out = deffield->region.nz;
    float x, y, z, value;
    vector<float> v;
    for (int k = 0; k < z_size_out; k++) {
        for (int j = 0; j < y_size_out; j++) {
            for (int i = 0; i < x_size_out; i++) {
                deffield->get(i, j, k, v);
                if (v.size() == 2) {
                    v.resize(3);
                    v[2] = 0;
                }
                x = static_cast<float>(i) + v[0];
                y = static_cast<float>(j) + v[1];
                z = static_cast<float>(k) + v[2];
                value = 0.0f;
                if (interpolate) {
                    value = image->value(x, y, z);
                } else if (x >= 0 && x <= (x_size_in - 1) &&
                           y >= 0 && y <= (y_size_in - 1) &&
                           z >= 0 && z <= (z_size_in - 1)) {
                    value = image->get(min(x_size_in, static_cast<int>(x + 0.5f)),
                                       min(y_size_in, static_cast<int>(y + 0.5f)),
                                       min(z_size_in, static_cast<int>(z + 0.5f)));
                }
                warpedimage->set(i, j, k, value);
            }
        }
    }
    return warpedimage;
}

// ---------------------------------------------------------------------------
Image* ConvertTo3DDeformationField(const Image* D, bool copy)
{
    if (D->GetNumberOfComponents() < 2 || D->GetNumberOfComponents() > 3) {
        BASIS_THROW(invalid_argument, "Deformation field must have 2 or 3 vector components!");
    }
    Image* out = new Image(D->region.nx, D->region.ny, D->region.nz,
                           D->hdr.datatype, 3, D->imgfmt);
    if (out == NULL) return NULL;
    out->CopyRegion     (D);
    out->CopyTransform  (D);
    out->CopyDataScaling(D);
    out->CopyMetaData   (D);
    if (D->GetNumberOfComponents() == 3 && !copy) {
        out->img     = D->img;
        out->ownsimg = false;
    } else {
        vector<float> v(2);
        for (int k = 0; k < D->region.nz; k++) {
            for (int j = 0; j < D->region.ny; j++) {
                for (int i = 0; i < D->region.nx; i++) {
                    D->get(i, j, k, v);
                    v.resize(3, 0.0f);
                    out->set(i, j, k, v);
                }
            }
        }
    }
    return out;
}

// ===========================================================================
// connected components
// ===========================================================================

// ---------------------------------------------------------------------------
Image* LabelConnectedComponents(const Image* mask, map<int, int>* sizes)
{
    return LabelConnectedComponents(mask, 0, sizes);
}

// ---------------------------------------------------------------------------
Image* LabelConnectedComponents(const Image* image, float threshold, map<int, int>* sizes)
{
    Ivector3d dim;
    dim.x = image->region.nx;
    dim.y = image->region.ny;
    dim.z = image->region.nz;

    Image* regions = new Image(dim.x, dim.y, dim.z, DT_SIGNED_INT, 1, image->imgfmt);
    if (regions == NULL) return NULL;
    regions->CopyRegion   (image);
    regions->CopyTransform(image);

    Image* visited = new Image(dim.x, dim.y, dim.z, DT_UNSIGNED_CHAR, 1, image->imgfmt);
    if (visited == NULL) {
        delete regions;
        return NULL;
    }

    int label = 1;
    Ivector3d seed;
    for (seed.z = 0; seed.z < dim.z; seed.z++) {
        for (seed.y = 0; seed.y < dim.y; seed.y++) {
            for (seed.x = 0; seed.x < dim.x; seed.x++) {
                if (image->get(seed.x, seed.y, seed.z) > threshold
                        && visited->get(seed.x, seed.y, seed.z) == 0) {
                    int n = grow_region(image, threshold, regions, visited, seed, label, true);
                    if (n == 0) {
                        delete regions;
                        return NULL;
                    }
                    if (sizes) (*sizes)[label] = n;
                    label++;
                }
            }
        }
    }

    delete visited;

    return regions;
}

// ===========================================================================
// mask generation
// ===========================================================================

// ---------------------------------------------------------------------------
Image* GenerateMask(const Image* image, float threshold)
{
    Image* mask = new Image(image->region.nx, image->region.ny, image->region.nz,
                            DT_UNSIGNED_CHAR, 1,
                            image->imgfmt);
    if (mask == NULL) return NULL;
    mask->CopyRegion   (image);
    mask->CopyTransform(image);
    for (int k = 0; k < image->region.nz; k++) {
        for (int j = 0; j < image->region.ny; j++) {
            for (int i = 0; i < image->region.nx; i++) {
                if (image->get(i, j, k) > threshold) {
                    mask->set(i, j, k, 255);
                } else {
                    mask->set(i, j, k, 0);
                }
            }
        }
    }
    return mask;
}

// ---------------------------------------------------------------------------
void RemoveNoiseFromMask(Image* mask)
{
    // label connected components and determine their sizes
    map<int, int> sizes;
    Image* regions = LabelConnectedComponents(mask, &sizes);
    // determine minimum size of connected components
    int       min_region_size = 400;
    const int min_image_size  = min(mask->region.nx, mask->region.ny);
    if      (min_image_size <  65) min_region_size =  80;
    else if (min_image_size < 120) min_region_size = 160;
    else if (min_image_size < 260) min_region_size = 250;
    // eliminate small components which are caused by noise
    for (int k = 0; k < mask->region.nz; k++) {
        for (int j = 0; j < mask->region.ny; j++) {
            for (int i = 0; i < mask->region.nx; i++) {
                int label = regions->value<int>(i, j, k);
                if (label > 0 && sizes[label] < min_region_size) {
                    mask->set(i, j, k, 0);
                }
            }
        }
    }
}

// ===========================================================================
// image registration
// ===========================================================================

// ---------------------------------------------------------------------------
void DetermineControlPointSpacing(const Image* mask,
                                  int& spacing_x, int& spacing_y, int& spacing_z,
                                  Image::Region& region, bool useLessMemory)
{
    // constants
    const int numXY       = 50;
    const int numZ        = 60;
    const int numXYbase   = numXY / 2 - 1;
    const int numZbase    = numZ  / 2 - 1;
    const int numXYbaseSm = numXY / 3 - 1;   
    const int numZbaseSm  = numZ  / 3 - 1;
    const int scaleFactor = 5;  // ideally should be 9, was 9 till 08/15/2012. Changed into 5 to reduce memory usage.

    const int minSpacingX = 3;
    const int minSpacingY = 3;
    const int minSpacingZ = 1;
	
	const int inspacing_x = spacing_x;
	const int inspacing_y = spacing_y;
	const int inspacing_z = spacing_z;

    // compute spacing only if padding is disabled
    if (region == mask->region) {
        if (inspacing_x <= 0) spacing_x = static_cast<int>((region.nx + numXYbase) / numXY);
        if (inspacing_y <= 0) spacing_y = static_cast<int>((region.ny + numXYbase) / numXY);
        if (inspacing_z <= 0) spacing_z = static_cast<int>((region.nz + numZbase)  / numZ);
        if (inspacing_x <= 0 && inspacing_y <= 0 && abs(spacing_x - spacing_y) == 1 && useLessMemory == false ) {
            spacing_x = (spacing_x < spacing_y) ? spacing_x : spacing_y;
            spacing_y = (spacing_x < spacing_y) ? spacing_x : spacing_y;
        }
        if (spacing_x < minSpacingX) spacing_x = minSpacingX;
        if (spacing_y < minSpacingY) spacing_y = minSpacingY;
        if (spacing_z < minSpacingZ) spacing_z = minSpacingZ;
        return;
    }

    // initial image region
    region = mask->region;

    // initial available margin
    int bLeft   = region.nx - 1;
    int bRight  = 0;
    int bTop    = region.ny - 1;
    int bBottom = 0;
    int bFront  = region.nz - 1;
    int bRear   = 0;

    for (int k = 0; k < region.nz; k++) {
        for (int j = 0; j < region.ny; j++) {
            for (int i = 0; i < region.nx; i++) {
                if (mask->get(i, j, k) > 0) {
                    if (i < bLeft)   bLeft   = i;
                    if (i > bRight)  bRight  = i;
                    if (j < bTop)    bTop    = j;
                    if (j > bBottom) bBottom = j;
                    if (k < bFront)  bFront  = k;
                    if (k > bRear)   bRear   = k;
                }
            }
        }
    }

    int marginLeft   = bLeft;
    int marginRight  = region.nx - 1 - bRight;
    int marginTop    = bTop;
    int marginBottom = region.ny - 1 - bBottom;
    int marginFront  = bFront;
    int marginRear   = region.nz - 1 - bRear;

    // initial control point spacing
    if (inspacing_x <= 0) spacing_x = static_cast<int>((region.nx + numXYbase) / numXY);
    if (inspacing_y <= 0) spacing_y = static_cast<int>((region.ny + numXYbase) / numXY);
    if (inspacing_z <= 0) spacing_z = static_cast<int>((region.nz + numZbase)  / numZ);
    if (inspacing_x <= 0 && inspacing_y <= 0 && abs(spacing_x - spacing_y) == 1 && useLessMemory == false ) {
        spacing_x = (spacing_x < spacing_y) ? spacing_x : spacing_y;
        spacing_y = (spacing_x < spacing_y) ? spacing_x : spacing_y;
    }
    if (spacing_x < minSpacingX) spacing_x = minSpacingX;
    if (spacing_y < minSpacingY) spacing_y = minSpacingY;
    if (spacing_z < minSpacingZ) spacing_z = minSpacingZ;

    // iteratively increase image region if required and update spacing of control
    // points until the image region is big enough to fit enough effective control points
    for (;;) {
        // required margin
        int requiredMarginLeftAndRight = scaleFactor * spacing_x;
        int requiredMarginTopAndBottom = scaleFactor * spacing_y;
        int requiredMarginFrontAndRear = scaleFactor * spacing_z;

        // required additional padding
        int paddingLeft   = marginLeft   < requiredMarginLeftAndRight ? requiredMarginLeftAndRight - marginLeft   : 0;
        int paddingRight  = marginRight  < requiredMarginLeftAndRight ? requiredMarginLeftAndRight - marginRight  : 0;
        int paddingTop    = marginTop    < requiredMarginTopAndBottom ? requiredMarginTopAndBottom - marginTop    : 0;
        int paddingBottom = marginBottom < requiredMarginTopAndBottom ? requiredMarginTopAndBottom - marginBottom : 0;
        int paddingFront  = marginFront  < requiredMarginFrontAndRear ? requiredMarginFrontAndRear - marginFront  : 0;
        int paddingRear   = marginRear   < requiredMarginFrontAndRear ? requiredMarginFrontAndRear - marginRear   : 0;

        // done if image region is big enough
        if (paddingLeft == 0 && paddingRight == 0 &&
                paddingTop == 0 && paddingBottom == 0 &&
                paddingFront == 0 && paddingRear == 0) {
            break;
        }

        // update image region
        region.ox -= paddingLeft;
        region.oy -= paddingTop;
        region.oz -= paddingFront;
        region.nx += paddingLeft  + paddingRight;
        region.ny += paddingTop   + paddingBottom;
        region.nz += paddingFront + paddingRear;

        // update available margin
        marginLeft   += paddingLeft;
        marginRight  += paddingRight;
        marginTop    += paddingTop;
        marginBottom += paddingBottom;
        marginFront  += paddingFront;
        marginRear   += paddingRear;

        // update control point spacing
        if (inspacing_x <= 0) spacing_x = max(spacing_x, static_cast<int>( (region.nx+numXYbaseSm) / numXY ));
        if (inspacing_y <= 0) spacing_y = max(spacing_y, static_cast<int>( (region.ny+numXYbaseSm) / numXY ));
        if (inspacing_z <= 0) spacing_z = max(spacing_z, static_cast<int>( (region.nz+numZbaseSm ) / numZ ));
        if (inspacing_x <= 0 && inspacing_y <= 0 && abs(spacing_x - spacing_y) == 1 && useLessMemory == false ) {
            spacing_x = (spacing_x < spacing_y) ? spacing_x : spacing_y;
            spacing_y = (spacing_x < spacing_y) ? spacing_x : spacing_y;
        }
        if (spacing_x < minSpacingX) spacing_x = minSpacingX;
        if (spacing_y < minSpacingY) spacing_y = minSpacingY;
        if (spacing_z < minSpacingZ) spacing_z = minSpacingZ;
    }
}


} // namespace dramms
