/**
 * @file  utilities.h
 * @brief Utility imaging functions.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#pragma once
#ifndef _DRAMMS_COMMON_UTILITIES_H
#define _DRAMMS_COMMON_UTILITIES_H


#include <map>
#include <string>

#include "mvcd.h"  // Fvector2d, Ivector2d
#include "image.h" // Image, nifti1_io.h


namespace dramms {


// ===========================================================================
// intensity scaling
// ===========================================================================

/**
 * @brief Get maximum intensity range of image datatype.
 *
 * @param [in]  datatype Datatype enumeration value, e.g., DT_UNSIGNED_CHAR.
 * @param [out] min      Minimum intensity that can be represented.
 * @param [out] max      Maximum intensity that can be represented.
 */
void GetDatatypeRange(short datatype, float& min, float& max);

/**
 * @brief Determine if image is an actual scalar byte image.
 *
 * An image with only integral intensities in the range [0, 255] is considered
 * an actual byte image no matter in which datatype the intensities are stored.
 *
 * @param [in] image Intensity image.
 * @param [in] scale Whether to apply scaling function before test if intensity
 *                   is of integral type. If @c false, this function only
 *                   returns true if the slope and intercept of the image
 *                   header map each intensity to itself, i.e., either
 *                   @c scl_slope = 0 or @c scl_slope = 1 and @c scl_inter = 0.
 *
 * @returns Whether the given image is a scalar byte image.
 */
bool IsByteImage(const Image* image, bool scale = true);

/**
 * @brief Get data scaling function for image intensities.
 *
 * This function returns the slope and intercept parameters of a linear scaling
 * function which maps the image intensities stored in memory to the desired
 * output range if @p lmin < @p lmax. Moreover, if the @p eff parameter is set
 * to @c true in this case, the histogram of image intensities is computed and
 * only intensities at the boundaries of the histogram (i.e., those corresponding
 * to background and noise) will be mapped outside the specified intensity range.
 * These should be set to the @p lmin or @p lmax intensity when rescaling the image.
 *
 * If the output datatype is signed, i.e., @p lmin < 0, the intercept is always 0
 * such that zero intensities are preserved and only intensities to the left or
 * right of the zero intensity are stretched/shrinked. This rescaling preserves
 * the sign of the intensities. Otherwise, if the output datatype is unsigned, i.e.,
 * @p lmin >= 0, also intensities are also shifted as needed in order to not
 * truncate any values.
 *
 * Otherwise, if @p lmin >= @p lmax, the scaling as defined by the image header
 * is returned.
 *
 * @note Actual byte images, i.e., images whose scaled intensities are all of
 *       integral value and in the range [0, 255], are treated slightly
 *       differently. In this case, the histogram is always computed over the
 *       fixed range [0, 255], no matter what the actual input range of
 *       intensities is. The reason for this is that it is what has been
 *       implemented by the pre-release versions of DRAMMS and the parameters
 *       of the consecutive processing steps were tuned according to this
 *       rescaling of the input intensities.
 *
 * @param [in]  image  Image to be rescaled.
 * @param [out] slope  Scaling slope.
 * @param [out] inter  Scaling intercept.
 * @param [in]  lmin   Minimum value of output range.
 * @param [in]  lmax   Maximum value of output range.
 * @param [in]  eff    Whether to map effective intensities only.
 */
void GetDataScaling(const Image* image, float& slope, float& inter,
                    float lmin = 0.0f, float lmax = 0.0f, bool eff = false);

/**
 * @brief Linearly rescale image intensities.
 *
 * This function rescales the intensities of the input image using the linear
 * function @c y = @p slope * @c x + @p inter and writes the rescaled intensities
 * to the given output image. If a value @c y is outside the range of values
 * which can be represented by the output data type, these values are truncated
 * and mapped to the closest representable value.
 *
 * The scl_slope and scl_inter values of the NIfTI-1 header of the output
 * image are recalculated to reflect any scaling or rescaling that has be
 * applied to the intensity values stored in memory. Hence, by applying this
 * scaling function, one can recover the values of the input image which
 * correspond to the intensities which result from scaling the input values
 * according to the scaling function of the NIfTI-1 header of the input image.
 *
 * @param [out] scaled Image with rescaled intensities.
 * @param [in]  image  Input image.
 * @param [in]  slope  Rescaling slope or 0 if scaling function defined by
 *                     header of input image should be used.
 * @param [in]  inter  Rescaling intercept.
 */
void ScaleImage(Image* scaled, const Image* image, float slope = 0.0f, float inter = 0.0f);

// ===========================================================================
// copy / conversion
// ===========================================================================

/**
 * @brief Cast image to another datatype.
 *
 * This function allocates a new image of the requested datatype and copies
 * the image data from the given image to the new image, casting the image
 * values as necessary.
 *
 * If @p scale is true, the scaling function defined by the scl_slope and
 * scl_inter values of the NIfTI-1 header is applied to the image intensities
 * before casting the values to the output data type. Moreover, if this cast
 * would result in the truncation of intensities, the scaling function is
 * further adjusted beforehand such that the output intensities are rescaled
 * to the maximum intensity range of the output data type.
 *
 * The scl_slope and scl_inter values of the NIfTI-1 header of the output
 * image are recalculated to reflect any scaling or rescaling that has be
 * applied to the intensity values stored in memory. Hence, by applying this
 * scaling function, one can recover the values of the input image which
 * correspond to the intensities which result from scaling the input values
 * according to the scaling function of the NIfTI-1 header of the input image.
 *
 * @param [in] image    Input image.
 * @param [in] datatype NIfTI-1 datatype code of casted image.
 * @param [in] scale    Whether to scale the intensities according to the
 *                      scaling function defined by the scl_slope and scl_inter
 *                      values of the NIfTI-1 header of the input image or
 *                      the maximum range of the output datatype if
 *                      otherwise intensities would be truncated.
 *
 * @returns New image of given datatype with optionally rescaled intensities.
 *          The returned image has to be destroyed by the caller using the
 *          @c delete operator.
 */
Image* CastImage(const Image* image, short datatype, bool scale = true);

/**
 * @brief Cast image to another datatype.
 *
 * This function allocates a new image of the requested datatype and copies
 * the image data from the given image to the new image, scaling and casting
 * the image values as necessary and requested.
 *
 * This function casts an image to the given output data type, applies the
 * scaling function as given by the NIfTI-1 header before, and further
 * rescales the intensities to the range specified by [@p min, @p max] if
 * the intensities would otherwise be truncated by the data type conversion.
 * If either @p min or @p max are outside the range of values representable
 * by the chosen output data type, the image values which are mapped to
 * values outside the range of the data type are truncated.
 *
 * The scl_slope and scl_inter values of the NIfTI-1 header of the output
 * image are recalculated to reflect any scaling or rescaling that has be
 * applied to the intensity values stored in memory. Hence, by applying this
 * scaling function, one can recover the values of the input image which
 * correspond to the intensities which result from scaling the input values
 * according to the scaling function of the NIfTI-1 header of the input image.
 *
 * @param [in] image    Input image.
 * @param [in] datatype NIfTI-1 datatype code of casted image.
 * @param [in] min      Minimum of output intensity range.
 * @param [in] max      Maximum of output intensity range.
 * @param [in] smooth   Whether to average the image with a smoothed version
 *                      of it before applying the scaling function and
 *                      casting the intensities.
 * @param [in] eff      Whether to only rescale intensities within the effective
 *                      range in order to avoid the 'long tail' phenomenon.
 *
 * @returns New image of given datatype with optionally rescaled intensities.
 *          The returned image has to be destroyed by the caller using the
 *          @c delete operator.
 */
Image* CastImage(const Image* image, short datatype, float min, float max,
                 bool smooth = false, bool eff = false);

/**
 * @brief Copy image.
 *
 * @returns New image copy. The returned image has to be
 *          destroyed by the caller using the @c delete operator.
 */
Image* CopyImage(const Image* image);

// ===========================================================================
// image statistics
// ===========================================================================

/**
 * @brief Get intensity range in scalar image.
 *
 * @param [in]  image Intensity image.
 * @param [out] min   Minimum intensity.
 * @param [out] max   Maximum intensity.
 * @param [in]  scale Whether to return scaled intensities using the scaling
 *                    function of the image header.
 * @param [in]  eff   Whether to determine the effective intensity range only
 *                    by considering the histogram of intensities.
 */
void GetIntensityRange(const Image* image, float& min, float& max,
                       bool scale = false, bool eff = false);

/**
 * @brief Compute histogram of intensities.
 *
 * If @p min is greater or equal to @p max, the actual intensity range of the
 * image is determined first and used instead.
 *
 * @param [in] image Intensity image.
 * @param [in] nbins Number of histogram bins.
 * @param [in] scale Whether to compute histogram of scaled intensities using
 *                   the scaling function specified by the image header.
 * @param [in] min   Intensity at center of first histogram bin.
 * @param [in] max   Intensity at center of last histogram bin.
 *
 * @returns Histogram of image intensities or NULL if memory allocation failed.
 */
float* ComputeHistogram(const Image* image, int nbins = 256, bool scale = false,
                        float min = 0.0f, float max = 0.0f);

/**
 * @brief Normalize histogram of intensities.
 *
 * @param [in,out] histo Histogram of image intensities which will be turned
 *                       in-place into a normalized histogram.
 * @param [in]     nbins Number of histogram bins.
 *
 * @returns Normalized histogram. Note that this pointer references
 *          the memory of the input histogram @p histo.
 */
float* NormalizeHistogram(float* histo, int nbins = 256);

/**
 * @brief Compute cumulative histogram of intensities.
 *
 * @param [in,out] histo Histogram of image intensities which will be turned
 *                       in-place into a cumulative histogram.
 * @param [in]     nbins Number of histogram bins.
 *
 * @returns Cumulative histogram. Note that this pointer references
 *          the memory of the input histogram @p histo.
 */
float* CumulateHistogram(float* histo, int nbins = 256);

/**
 * @brief Compute the effective threshold to maximumly remove background.
 *
 * @param [in] histo Histogram of image intensities.
 * @param [in]     nbins Number of histogram bins.
 * @param [out]    the index of bin as the effective threshold.
 *
 * @returns the index of bin as the effective threshold.
 */
int FindAdaptiveThreshold(float* histo, int nbins = 256);

// ===========================================================================
// convolution
// ===========================================================================

/**
 * @brief Get parameters of Gaussian kernel given image size.
 *
 * These Gaussian kernel parameters are used for smoothing an intensity image.
 *
 * @param [in]  size   Image size.
 * @param [out] radius Kernel radius.
 * @param [out] sigma  Sigma value.
 */
void GetGaussParameters(int size, int& radius, float& sigma);

/**
 * @brief Get parameters of Gaussian kernal given voxel size.
 *
 * These Gaussian kernel parameters are used for smoothing an intensity image
 * before (re-)scaling the image intensities.
 *
 * @param [in]  pixdim Voxel size.
 * @param [out] radius Kernel radius.
 * @param [out] sigma  Sigma value.
 */
void GetGaussParameters(Fvector2d& pixdim, Ivector2d& radius, Fvector2d& sigma);

/**
 * @brief Get parameters of Gaussian kernal given voxel size.
 *
 * These Gaussian kernel parameters are used for smoothing an intensity image
 * before (re-)scaling the image intensities.
 *
 * @param [in]  pixdim Voxel size.
 * @param [out] radius Kernel radius.
 * @param [out] sigma  Sigma value.
 */
void GetGaussParameters(Fvector3d& pixdim, Ivector3d& radius, Fvector3d& sigma);

/**
 * @brief Get two-dimensional Gaussian kernel.
 *
 * @param [in] radius Radius of Gauss function in each dimension.
 * @param [in] sigma  Standard deviation of Gauss function in each dimension.
 * @param [in] fmt    Image format used to represent Gauss function.
 *
 * @returns Gaussian kernel image of size (2 * radius.x + 1) x (2 * radius.y + 1)
 *          or NULL if not enough memory available. The returned kernel has to
 *          be destroyed by the caller using the delete operator.
 */
Image* Gauss(const Ivector2d& radius,
             const Fvector2d& sigma,
             Image::Format    fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Get three-dimensional Gaussian kernel.
 *
 * @param [in] radius Radius of Gauss function in each dimension.
 * @param [in] sigma  Standard deviation of Gauss function in each dimension.
 * @param [in] fmt    Image format used to represent Gauss function.
 *
 * @returns Gaussian kernel image of size
 *          (2 * radius.x + 1) x (2 * radius.y + 1) x (2 * radius.z + 1)
 *          or NULL if not enough memory available. The returned kernel has to
 *          be destroyed by the caller using the delete operator.
 */
Image* Gauss(const Ivector3d& radius,
             const Fvector3d& sigma,
             Image::Format    fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Convolve image with given kernel.
 *
 * @param [in]  image  Input image.
 * @param [in]  kernel Convolution kernel. May have less dimensions than
 *                     the input image. For example, a 2D kernel is applied
 *                     to a 3D image, slice by slice.
 * @param [out] output Result of convolution. The output image has to be
 *                     allocated before and the image must be at least as
 *                     big as the input image. Usually, it should have the
 *                     same size as the input image.
 */
void ConvolveImage(const Image* image, const Image* kernel, Image* output);

/**
 * @brief Convolve image with given kernel.
 *
 * @param [in]  image  Input image.
 * @param [in]  kernel Convolution kernel. May have less dimensions than
 *                     the input image. For example, a 2D kernel is applied
 *                     to a 3D image, slice by slice.
 *
 * @returns Result of convolution or NULL if memory allocation failed.
 *          The returned image has to be destroyed using the delete operator.
 */
Image* ConvolveImage(const Image* image, const Image* kernel);

// ===========================================================================
// image regions
// ===========================================================================

/**
 * @brief Determine the foreground region.
 *
 * @param [in] image     Image.
 * @param [in] threshold Upper threshold of background values, i.e., all voxels
 *                       with an intensity greater than this threshold are
 *                       considered to be foreground.
 *
 * @returns Region of foreground.
 */
Image::Region GetForegroundRegion(const Image* image, float threshold = 0.0f);

/**
 * @brief Join image regions.
 *
 * @param [in] region1 First image region.
 * @param [in] region2 Second image region.
 *
 * @returns The union of both regions.
 */
Image::Region JoinRegions(const Image::Region& region1, const Image::Region& region2);

// ===========================================================================
// image grid manipulation
// ===========================================================================

/**
 * @brief Permute image axes.
 *
 * @param [in,out] image Image whose axes are to be permuted.
 * @param [in]     order Array of permuted axes indices. The indices must be
 *                       a rearrangement of the indices 0, 1, 2, where 0 is the
 *                       first axis of the input image, 1 the second axis, and
 *                       2 the third axis.
 *
 * @sa itk::PermuteAxesImageFilter
 *
 * @returns Whether input order was valid and enough temporary memory could
 *          be allocated to perform the permutation.
 */
bool PermuteAxes(Image* image, int order[3]);

/**
 * @brief Flip image axes.
 *
 * @param [in] image         Image whose axes are to be flipped.
 * @param [in] flip          Array of a boolean value for each axis. Each axis for
 *                           which the boolean value is true, is flipped.
 * @param [in] update_offset Whether to change the qoffset such that
 *                           it still corresponds to voxel [0,0,0]
 *                           of the input image, i.e., before any flip.
 *
 * @returns Whether enough temporary memory could be allocated to perform
 *          the flip of the axes.
 */
bool FlipAxes(Image* image, bool flip[3], bool update_offset = true);

/**
 * @brief Orient image.
 *
 * Valid orientation codes are NIFTI_R2L, NIFTI_L2R, NIFTI_A2P, NIFTI_P2A,
 * NIFTI_I2S, and NIFTI_S2I as defined in nifti1_io.h. The input image can
 * be a scalar image or a deformation field. In case of a deformation field,
 * the components of the displacement vectors are reoriented accordingly.
 *
 * @param [in,out] image         Image to be oriented.
 * @param [in]     orient        Array of orientation codes.
 * @param [in]     update_offset Whether to change the qoffset such that
 *                               it still corresponds to voxel [0,0,0]
 *                               of the input image, i.e., before the reorientation.
 *
 * @returns Whether desired orientation was valid and enough temporary memory
 *          could be allocated to perform the reorientation.
 */
bool OrientImage(Image* image, int orient[3], bool update_offset = true);

/**
 * @brief Orient image.
 *
 * Valid orientation codes are NIFTI_R2L, NIFTI_L2R, NIFTI_A2P, NIFTI_P2A,
 * NIFTI_I2S, and NIFTI_S2I as defined in nifti1_io.h. The input image can
 * be a scalar image or a deformation field. In case of a deformation field,
 * the components of the displacement vectors are reoriented accordingly.
 *
 * @param [in,out] image         Image to be oriented.
 * @param [in]     x_orient      Orientation along x axis.
 * @param [in]     y_orient      Orientation along y axis.
 * @param [in]     z_orient      Orientation along z axis.
 * @param [in]     update_offset Whether to change the qoffset such that
 *                               it still corresponds to voxel [0,0,0]
 *                               of the input image, i.e., before the reorientation.
 *
 * @returns Whether desired orientation was valid and enough temporary memory
 *          could be allocated to perform the reorientation.
 */
inline bool OrientImage(Image* image, int x_orient, int y_orient, int z_orient, bool update_offset = true)
{
    int orient[3] = {x_orient, y_orient, z_orient};
    return OrientImage(image, orient, update_offset);
}

/**
 * @brief Downsample image.
 *
 * @param [in] image Input image.
 * @param [in] ratio Downsample ratio.
 *
 * @returns Downsampled image or NULL if not enough memory was available.
 *          The returned image has to be destroyed by the caller using the
 *          delete operator.
 *
 * @throws std::invalid_argument if downsampling ratio is less than 1 or if
 *                               image size is smaller than this ratio in
 *                               at least one dimension.
 */
Image* DownsampleImage(const Image* image, int ratio = 2);

/**
 * @brief Smooth image and then downsample it.
 *
 * This function smoothes the given image using a two-dimensional Gaussian
 * kernel with automatically selected parameters depending on the image size
 * and the downsample ratio. If the input image is three-dimensional, the
 * image is smoothed slice-by-slice. In case of a vector field, the components
 * are smoothed one-by-one.
 */
Image* SmoothAndDownsampleImage(const Image* image, int ratio = 2);

/**
 * @brief Resize image.
 *
 * Resizing the image corresponds to cropping parts of the image and padding
 * other parts with a constant value such that the resulting image is within
 * the bounds specified by the desired image region.
 *
 * If @p subdomain is @c false, the image header is updated in order to reflect
 * the new size of the image and the image @c region attribute is set to the
 * entire image region of the resized image. Otherwise, the image header is
 * not modified and the image @c region attribute set to the specified
 * @p region argument. Thus, if @p subdomain is @c true, it is possible to
 * resize the previously resized image again, specifying the desired region
 * within the original image domain. By not specifying any valid region, the
 * previously image is resized such that the resulting image data is defined on
 * the grid of the original image again. This is for example useful if an image
 * has been resized to concentrate on the actual foreground region only
 * (i.e., saving running time and in particular memory), but the final result
 * shall be defined within the original image domain. Here, by just setting
 * @p subdomain to @c true, the WriteNiftiImage() function will automatically
 * resize the image again such that the written image data is defined on the
 * image grid defined by the image header.
 *
 * Example:
 * @code
 * Image* tmp = NULL;
 * Image* image = ReadNiftiImage("in.nii");
 * // resize image to foreground region only
 * Image::Region roi = GetForegroundRegion(image);
 * tmp = ResizeImage(image, roi, 0, true);
 * delete image;
 * image = tmp;
 * // process image
 * // ...
 * // write image data within original image domain
 * // see documentation of WriteNiftiImage() for details
 * WriteNiftiImage("out.nii", image);
 * delete image;
 * @endcode
 *
 * @param [in] image     Image whose data is resized.
 * @param [in] region    The desired image region within the image bounds
 *                       specified by the header information. If the number
 *                       of voxels in a dimension is zero, the image is
 *                       resized such that the region in this dimension
 *                       corresponds to the image size as specified by the
 *                       image header.
 * @param [in] pad_value The value used to pad the image if necessary.
 * @param [in] subdomain Whether the resized image shall be treated as an image
 *                       defined on a subdomain of the original image domain.
 *                       See description of this function for more information.
 *
 * @returns Resized image or NULL if memory allocation failed.
 */
Image* ResizeImage(const Image*  image,
                   Image::Region region    = Image::Region(),
                   float         pad_value = 0.0f,
                   bool          subdomain = false);

/**
 * @brief Resample image on specified image grid.
 *
 * @param [in] image  Input image.
 * @param [in] pixdim Size of voxels of output image.
 * @param [in] region Image region of output image relative to input image.
 *
 * @returns Resampled image or NULL if memory allocation failed.
 *
 * @throws std::invalid_argument if @p pixdim is invalid.
 */
Image* ResampleImage(const Image*  image,
                     Fvector3d     pixdim,
                     Image::Region region = Image::Region());

// ===========================================================================
// smoothing
// ===========================================================================

/**
 * @brief Smooth image, both scalar or vector field.
 *
 * This function smoothes the given image using a two-dimensional Gaussian
 * kernel with automatically selected parameters depending on the image size.
 * If the input image is three-dimensional, the image is smoothed
 * slice-by-slice. In case of a vector field, the components are smoothed
 * one-by-one.
 *
 * @param [in] image Input image.
 *
 * @returns Smoothed image of type DT_FLOAT or NULL if not enough memory
 *          was available. The returned image has to be destroyed by the
 *          caller using the delete operator. 
 */
Image* SmoothImage(const Image* image);

// ===========================================================================
// transformations
// ===========================================================================

/**
 * @brief Transform point using affine transformation matrix.
 *
 * @param [in] m Affine transformation matrix. The last row of this matrix
 *               is considered to be (0, 0, 0, 1) and thus ignored.
 * @param [in] p Coordinates of point to transform.
 *
 * @returns Coordinates of point after transformation.
 */
Fvector3d TransformPoint(const Image::Transform m, const Fvector3d& p);

/**
 * @brief Get affine identity transform.
 *
 * @returns 4x4 identity matrix.
 */
Image::Transform IdentityTransform();

/**
 * @brief Get inverse transformation.
 *
 * @param [in] T Affine transformation matrix.
 *
 * @returns Inverse affine transformation.
 */
Image::Transform InvertTransform(const Image::Transform& T);

/**
 * @brief Estimate inverse deformation field.
 *
 * This function estimates the inverse deformation field of a two-dimensional
 * or three-dimensional deformation field.
 *
 * @param [in] D           Deformation field in DRAMMS format.
 * @param [in] num_samples Number of samples used around each voxel of the
 *                         input deformation field to obtain a smooth estimate.
 *                         The more samples are used, the higher the input
 *                         deformation field is sampled. Set to less or equal
 *                         to zero to disable subsampling of the input
 *                         deformation field.
 * @param [in] boundary    Number of additional boundary voxels to use to
 *                         estimate inverse deformation field. The more,
 *                         the smoother the resulting deformation field will
 *                         be at the boundary, but the more memory is required
 *                         and the computation time is increased.
 *
 * @returns Inverse deformation field.
 *
 * @throws std::invalid_argument if the input deformation field image is no
 *                               valid deformation field or not in the DRAMMS
 *                               format. Use the Image::SetFormat() method
 *                               to change the format of the image to DRAMMS.
 */
Image* InvertTransform(const Image* D, int num_samples = 2, int boundary = 8);

/**
 * @brief Concatenate two affine transformations.
 *
 * @note The first transformation is applied first, then the second one, i.e.,
 *       the matrix multiplication is carried out in the order T2 * T1, NOT T1 * T2.
 *
 * @param [in] T1 First affine transformation.
 * @param [in] T2 Second affine transformation.
 *
 * @returns Transformation matrix T2 * T1, i.e., [T2 o T1](x) = T2[T1(x)].
 */
Image::Transform ConcatenateTransforms(const Image::Transform& T1, const Image::Transform& T2);

/**
 * @brief Concatenate affine and deformable transformation.
 *
 * @attention All inputs must be given in DRAMMS order of the axes, i.e.,
 *            yxz compared to the NIfTI order xyz. This includes the header
 *            of image A and the affine transformation. Use PermuteXY() to
 *            change the order of the read NIfTI header and the FSL matrix.
 *            The format of the deformation field can be changed to the
 *            DRAMMS format using the Image::SetFormat() method.
 *
 * @note This function requires that the image space B is identical to image space C.
 *
 * @param [in] A Reference image header defining space of image A.
 * @param [in] T Affine transformation from image A to B.
 * @param [in] D Deformable transformation from image B to C.
 *
 * @returns Composition of transformations, i.e., [D o T](x) = D[T(x)],
 *          or NULL if memory allocation failed.
 *
 * @throws std::invalid_argument if the input deformation field image is no
 *                               valid deformation field.
 */
Image* ConcatenateTransforms(const Image::Header&    A,
                             const Image::Transform& T,
                             const Image*            D);

/**
 * @brief Concatenate affine and deformable transformation.
 *
 * @attention All inputs must be given in DRAMMS order of the axes, i.e.,
 *            yxz compared to the NIfTI order xyz. This includes the header
 *            of image A and the affine transformation. Use PermuteXY() to
 *            change the order of the read NIfTI header and the FSL matrix.
 *            The format of the deformation field can be changed to the
 *            DRAMMS format using the Image::SetFormat() method.
 *
 * @param [in] A Reference image header defining space of image A.
 * @param [in] B Reference image header defining space of image B.
 * @param [in] T Affine transformation from image A to B.
 * @param [in] D Deformable transformation from image B to C.
 *
 * @returns Composition of transformations, i.e., [D o T](x) = D[T(x)],
 *          or NULL if memory allocation failed.
 *
 * @throws std::invalid_argument if the input deformation field image is no
 *                               valid deformation field or not in the DRAMMS
 *                               format. Use the Image::SetFormat() method
 *                               to change the format of the image to DRAMMS.
 *                               Further, the dimension of the displacement
 *                               vectors in @p D1 and @p D2 must be identical.
 */
Image* ConcatenateTransforms(const Image::Header&    A,
							 const Image::Header&    B,
                             const Image::Transform& T,
                             const Image*            D);

/**
 * @brief Concatenate two deformations.
 *
 * This function calculates the composition of the deformation fields @p D1
 * and @p D2, i.e., @p D2 o @p D1. Here, @p D1 is usually the deformation obtained
 * when registering image A to image B, and @p D2 the result of registering
 * image B to image C. Then the composition will generate a deformation that warps
 * image A to the space of image C. Note that the two deformation fields can
 * have different size.
 *
 * @param [in] D1 Deformation in DRAMMS format from image space A to B.
 * @param [in] D2 Deformation in DRAMMS format from image space B to C.
 *
 * @returns Composition of deformations, i.e., [D2 o D1](x) = D2[D1(x)],
 *          or NULL if memory allocation failed.
 *
 * @throws std::invalid_argument if the input deformation field image is no
 *                               valid deformation field or not in the DRAMMS
 *                               format. Use the Image::SetFormat() method
 *                               to change the format of the image to DRAMMS.
 *                               Further, the dimension of the displacement
 *                               vectors in @p D1 and @p D2 must be identical.
 */
Image* ConcatenateTransforms(const Image* D1, const Image* D2);

/**
 * @brief Subtract affine from deformable transformation.
 *
 * @attention All inputs must be given in DRAMMS order of the axes, i.e.,
 *            yxz compared to the NIfTI order xyz. This includes the header
 *            of image A and the affine transformation. Use PermuteXY() to
 *            change the order of the read NIfTI header and the FSL matrix.
 *            The format of the deformation field can be changed to the
 *            DRAMMS format using the Image::SetFormat() method.
 *
 * @param [in] A Reference image header defining space of image A.
 * @param [in] D Deformable transformation from image A to C.
 * @param [in] T Affine transformation from image A to B.
 *
 * @returns Deformable transformations from B to C,
 *          or NULL if memory allocation failed.
 *
 * @throws std::invalid_argument if the input deformation field image is no
 *                               valid deformation field.
 */
Image* SubtractTransforms(const Image::Header&    A,
                          const Image*            D,
                          const Image::Transform& T);

/**
 * @brief Subtract affine from deformable transformation.
 *
 * @attention All inputs must be given in DRAMMS order of the axes, i.e.,
 *            yxz compared to the NIfTI order xyz. This includes the header
 *            of image A and the affine transformation. Use PermuteXY() to
 *            change the order of the read NIfTI header and the FSL matrix.
 *            The format of the deformation field can be changed to the
 *            DRAMMS format using the Image::SetFormat() method.
 *
 * @param [in] A Reference image header defining space of image A.
 * @param [in] B Reference image header defining space of image B.
 * @param [in] D Deformable transformation from image A to C.
 * @param [in] T Affine transformation from image A to B.
 *
 * @returns Deformable transformations from B to C,
 *          or NULL if memory allocation failed.
 *
 * @throws std::invalid_argument if the input deformation field image is no
 *                               valid deformation field.
 */
Image* SubtractTransforms(const Image::Header&    A,
						  const Image::Header&    B,
                          const Image*            D,
                          const Image::Transform& T);

/**
 * @brief Subtract deformation from deformable transformation.
 *
 * @param [in] D1 Unified deformable transformation from image A to C.
 * @param [in] D2 Affine transformation from image A to B.
 *
 * @returns Deformable transformation from B to C, or NULL if memory allocation failed.
 *
 * @throws std::invalid_argument if any input deformation field image is no
 *                               valid deformation field or if the dimensions
 *                               of the deformation fields do not match.
 */
Image* SubtractTransforms(const Image* D1, const Image* D2);

/**
 * @brief Linearly transform scalar image.
 *
 * @param [in] image       Subject image.
 * @param [in] T           Linear transformation matrix.
 * @param [in] reference   Reference/template image.
 * @param [in] interpolate Whether to use linear interpolation.
 *                         If false, nearest neighbor interpolation is used.
 *
 * @returns Transformed subject image in the space of the template image or
 *          NULL if memory allocation failed.
 */
Image* ApplyTransform(const Image*            image,
                      const Image::Transform& T,
                      const Image*            reference,
                      bool                    interpolate = true);

/**
 * @brief Apply deformation field on scalar image.
 *
 * @param [in] image       Subject image.
 * @param [in] deffield    Deformation field which warps the subject
 *                         image to the space of the template image.
 *                         Note that the parameters of the template space
 *                         are stored in the header of the deformation field.
 * @param [in] interpolate Whether to use linear interpolation.
 *                         If false, nearest neighbor interpolation is used.
 *
 * @returns Warped subject image in the space of the template image or
 *          NULL if memory allocation failed.
 */
Image* ApplyTransform(const Image* image,
                      const Image* deffield,
                      bool         interpolate = true);

/**
 * @brief Convert a 2D deformation field to a 3D deformation field.
 *
 * This function copies the input deformation field to one with three components
 * per voxel, where the displacement in the third component is set to zero.
 *
 * @param [in] D    Input deformation field (2D or 3D).
 * @param [in] copy Whether to copy the image data even if the input deformation
 *                  field is already three-dimensional.
 *
 * @returns 3D deformation field. If @p copy if @c false and the input deformation
 *          field is three-dimensional, the returned image instance uses the same
 *          image data memory which is not owned by this instance.
 */
Image* ConvertTo3DDeformationField(const Image* D, bool copy = true);

// ===========================================================================
// connected components
// ===========================================================================

/**
 * @brief Label connected components.
 *
 * @param [in]  mask  Mask image.
 * @param [out] sizes Determined sizes of connected components.
 *                    If NULL, the sizes are not returned.
 *
 * @returns Label image of datatype DT_SIGNED_INT where voxels belonging to
 *          a connected component are labeld with a unique label. Voxels
 *          which are background are labeled zero.
 */
Image* LabelConnectedComponents(const Image* mask, std::map<int, int>* sizes = NULL);

/**
 * @brief Label connected components.
 *
 * @param [in]  image     Image.
 * @param [in]  threshold Upper background threshold. Voxels with a value
 *                        exceeding this threshold are considered foreground.
 * @param [out] sizes     Determined sizes of connected components.
 *                        If NULL, the sizes are not returned.
 *
 * @returns Label image of datatype DT_SIGNED_INT where voxels belonging to
 *          a connected component are labeld with a unique label. Voxels
 *          which are background are labeled zero.
 */
Image* LabelConnectedComponents(const Image* image, float threshold, std::map<int, int>* sizes = NULL);

// ===========================================================================
// mask generation
// ===========================================================================

/**
 * @brief Generate mask from intensity image.
 *
 * @param [in] image     Intensity image.
 * @param [in] threshold Upper background threshold. Voxels with a value
 *                       exceeding this threshold are considered foreground.
 *
 * @returns Generated mask of type DT_UNSIGNED_CHAR with zero values at
 *          background voxels and 255 at foreground voxels or NULL if memory
 *          allocation failed.
 */
Image* GenerateMask(const Image* image, float threshold = 0.0f);

/**
 * @brief Remove noise from mask.
 *
 * This function determine the size of the connected components and removes
 * those connected components from the mask which are too small. The minimum
 * size of the connected components is hereby determined based on the size of
 * the mask image.
 *
 * @param [in, out] mask Mask image.
 */
void RemoveNoiseFromMask(Image* mask);

// ===========================================================================
// image registration
// ===========================================================================

/**
 * @brief Determine control point spacing and required image size.
 *
 * If the input spacing for a dimension is already set to a positive value,
 * it will not be modified by this function. Otherwise, this function will
 * determine the optimal spacing for the control points in this dimension.
 * Similarly, it the input region is set to the region of the @p mask image,
 * the image region is not increased by this function. Otherwise, the
 * image region is increased until enough control points fit into it.
 *
 * @param [in]     mask          A mask for foreground defined in template space.
 * @param [in,out] spacing_x     Space between neighboring control points in x.
 * @param [in,out] spacing_y     Space between neighboring control points in y.
 * @param [in,out] spacing_z     Space between neighboring control points in z.
 * @param [out]    region        Required image region. This region is either identical
 *                               to the image region of the mask if no padding of the
 *                               image domain is required or a region bigger than the
 *                               original image domain.
 * @param [in]     useLessMemory Whether to choose a spacing which requires less memory.
 */
void DetermineControlPointSpacing(const Image*   mask,
                                  int&           spacing_x,
                                  int&           spacing_y,
                                  int&           spacing_z,
                                  Image::Region& region,
								  bool 			 useLessMemory);


} // namespace dramms


#endif // _DRAMMS_COMMON_UTILITIES_H
