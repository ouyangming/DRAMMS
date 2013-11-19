/**
 * @file  imageio.h
 * @brief Image I/O functions.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#pragma once
#ifndef _DRAMMS_COMMON_IMAGEIO_H
#define _DRAMMS_COMMON_IMAGEIO_H


#include <string>

#include "image.h"


namespace dramms {


// ===========================================================================
// image - either one of the supported formats below
// ===========================================================================

/**
 * @brief Check whether a given image file can be read by ReadImage().
 *
 * @param [in] filename Name of image file.
 *
 * @returns Whether ReadImage() can read the given image file.
 */
bool CanReadImage(const char* filename);

/**
 * @brief Read image from file.
 *
 * The supported file formats are ANALYZE 7.5, NIfTI-1, and the Meta image format.
 *
 * @param [in] filename Name of image file.
 * @param [in] fmt      Image format.
 *
 * @returns Image of datatype as defined by the image header or NULL on error.
 *          The returned image has to be destroyed using the @c delete operator.
 *
 * @sa ReadNiftiImage(const char*, Image::Format)
 * @sa ReadMetaImage(const char*, Image::Format)
 */
Image* ReadImage(const char* filename, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Read image from file.
 *
 * The supported file formats are ANALYZE 7.5, NIfTI-1, and the Meta image format.
 *
 * See documentation of the CastImage() function for details on the rescaling
 * of intensity values that is optionally applied by this function when reading
 * the input image.
 *
 * @param [in] filename Name of image file.
 * @param [in] datatype NIfTI-1 datatype code, e.g., DT_UNSIGNED_CHAR or
 *                      NIFTI_TYPE_UINT8. The image data, which may be of
 *                      other datatype, is casted to the desired output
 *                      datatype using CastImage().
 * @param [in] scale    See @p scale parameter of CastImage() function.
 * @param [in] fmt      Image format.
 *
 * @returns Image or NULL on error. The returned image has to be destroyed
 *          by the caller using the @c delete operator.
 *
 * @sa ReadNiftiImage(const char*, short, bool, Image::Format)
 * @sa ReadMetaImage(const char*, short, bool, Image::Format)
 * @sa CastImage(const Image*, short, bool)
 */
Image* ReadImage(const char*   filename,
                 short         datatype,
                 bool          scale = false,
                 Image::Format fmt   = Image::FORMAT_DRAMMS);

/**
 * @brief Read image sequence from file.
 *
 * The supported file formats are ANALYZE 7.5, and NIfTI-1.
 *
 * @param [in] filename Name of image file.
 * @param [in] fmt      Image format.
 *
 * @returns Images of datatype as defined by the image header or empty sequence on error.
 *          The returned image has to be destroyed using the @c delete operator.
 *
 * @sa ReadNiftiSequence(const char*, Image::Format)
 */
Sequence ReadSequence(const char* filename, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Read image sequence from file.
 *
 * The supported file formats are ANALYZE 7.5, and NIfTI-1.
 *
 * See documentation of the CastImage() function for details on the rescaling
 * of intensity values that is optionally applied by this function when reading
 * the input image.
 *
 * @param [in] filename Name of image file.
 * @param [in] datatype NIfTI-1 datatype code, e.g., DT_UNSIGNED_CHAR or
 *                      NIFTI_TYPE_UINT8. The image data, which may be of
 *                      other datatype, is casted to the desired output
 *                      datatype using CastImage().
 * @param [in] scale    See @p scale parameter of CastImage() function.
 * @param [in] fmt      Image format.
 *
 * @returns Image sequence or empty sequence on error. The returned images have
 *          to be destroyed by the caller using the @c delete operator.
 *
 * @sa ReadNiftiSequence(const char*, short, bool, Image::Format)
 * @sa CastImage(const Image*, short, bool)
 */
Sequence ReadSequence(const char*   filename,
                      short         datatype,
                      bool          scale = false,
                      Image::Format fmt   = Image::FORMAT_DRAMMS);

/**
 * @brief Read deformation field from NIfTI-1, ANALYZE 7.5, or Meta image file.
 *
 * @param [in] filename Name of image file.
 * @param [in] fmt      Image format.
 *
 * @returns Deformation field or NULL on error.
 *
 * @sa ReadNiftiDeformationField()
 * @sa ReadMetaDeformationField()
 */
Image* ReadDeformationField(const char*   filename,
                            Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Write image to file.
 *
 * @param [in] filename Name of image file.
 * @param [in] image    Image to write.
 * @param [in] fmt      Format of image file. If Image::FILE_FORMAT_UNKNOWN,
 *                      the default, the format of the @c filefmt attribute of
 *                      the @p image is used.
 *
 * @returns Whether the image was written to file.
 *
 * @sa WriteImage(const char*, const Image*, const Image*, Image::FileFormat)
 */
bool WriteImage(const char*       filename,
                const Image*      image,
                Image::FileFormat fmt = Image::FILE_FORMAT_UNKNOWN);

/**
 * @brief Write image to file.
 *
 * The image format used to write the image to a file depends on the value of
 * the @c filefmt attribute of the @p image instance which is usually set by
 * the ReadImage*() functions to the format of the file from which an image
 * has been read from before it was processed. Optionally, a different file
 * format can be specified instead using the @p fmt parameter.
 *
 * If the image to write is only defined within a certain region of the
 * image domain as defined by the image header, the image is first
 * resized by this function such that the image data is defined on the
 * image domain specified by the image header before the image is written
 * to disk. Therefore, the ResizeImage() function is called with only the
 * image as argument. See also documentation of ResizeImage().
 *
 * If the image shall actually be written with a different size (image domain),
 * the @c subdomain argument of the ResizeImage() function has to be set to
 * @c false when resizing the image, or the image @c region attribute be set as
 * follows before writing the image to disk:
 * @code
 * image->region.ox = image->region.oy = image->region.oz = 0;
 * image->region.nx = image->hdr.dim[1];
 * image->region.ny = image->hdr.dim[2];
 * image->region.nz = image->hdr.dim[3];
 * @endcode
 *
 * @param [in] filename Name of image file.
 * @param [in] image    Image to write.
 * @param [in] mask     Mask image. Zero is written at voxels outside mask.
 * @param [in] fmt      Format of image file. If Image::FILE_FORMAT_UNKNOWN,
 *                      the default, the format of the @c filefmt attribute of
 *                      the @p image is used.
 *
 * @returns Whether the image was written to file.
 *
 * @sa WriteNiftiImage()
 * @sa WriteMetaImage()
 */
bool WriteImage(const char*       filename,
                const Image*      image,
                const Image*      mask,
                Image::FileFormat fmt = Image::FILE_FORMAT_UNKNOWN);

/**
 * @brief Write image sequence to file.
 *
 * @param [in] filename Name of image file.
 * @param [in] images   Image sequence to write.
 * @param [in] fmt      Format of image file. If Image::FILE_FORMAT_UNKNOWN,
 *                      the default, the format of the @c filefmt attribute of
 *                      the @p image is used.
 *
 * @returns Whether the image was written to file.
 *
 * @sa WriteSequence(const char*, const ConstSequence, const Image*, Image::FileFormat)
 */
bool WriteSequence(const char*         filename,
                   const ConstSequence images,
                   Image::FileFormat   fmt = Image::FILE_FORMAT_UNKNOWN);

/**
 * @brief Write image sequence to file.
 *
 * The image format used to write the images to a file depends on the value of
 * the @c filefmt attribute of the first @p image instance which is usually set
 * by the ReadImage*() functions to the format of the file from which an image
 * has been read from before it was processed. Optionally, a different file
 * format can be specified instead using the @p fmt parameter.
 *
 * If a image to write is only defined within a certain region of the
 * image domain as defined by the image header, the image is first
 * resized by this function such that the image data is defined on the
 * image domain specified by the image header before the image is written
 * to disk. Therefore, the ResizeImage() function is called with only the
 * image as argument. See also documentation of ResizeImage().
 *
 * If the image shall actually be written with a different size (image domain),
 * the @c subdomain argument of the ResizeImage() function has to be set to
 * @c false when resizing the image, or the image @c region attribute be set as
 * follows before writing the image to disk:
 * @code
 * image->region.ox = image->region.oy = image->region.oz = 0;
 * image->region.nx = image->hdr.dim[1];
 * image->region.ny = image->hdr.dim[2];
 * image->region.nz = image->hdr.dim[3];
 * @endcode
 *
 * @param [in] filename Name of image file.
 * @param [in] images   Image sequence to write.
 * @param [in] mask     Mask image. Zero is written at voxels outside mask.
 * @param [in] fmt      Format of image file. If Image::FILE_FORMAT_UNKNOWN,
 *                      the default, the format of the @c filefmt attribute of
 *                      the @p image is used.
 *
 * @returns Whether the image sequence was written to file.
 *
 * @sa WriteNiftiSequence()
 */
bool WriteSequence(const char*         filename,
                   const ConstSequence image,
                   const Image*        mask,
                   Image::FileFormat   fmt = Image::FILE_FORMAT_UNKNOWN);

/**
 * @brief Read data from raw image data file.
 *
 * This function should rarely be used. Instead, consider the use of either
 * ReadImage() or one of the specialized functions for a particular image
 * file format, i.e., ReadNiftiImage(), or ReadMetaImage().
 *
 * @param [in]      filename Name of image file.
 * @param [in, out] image    Image structure with allocated image buffer.
 *                           The size of the image buffer must be set, i.e.,
 *                           the image header information and image region
 *                           must be set properly before calling this function.
 *
 * @returns Whether the image data was read successfully.
 *
 * @sa ReadImage()
 * @sa ReadNiftiImage()
 * @sa ReadMetaImage()
 */
bool ReadRawImage(const char* filename, Image* image);

/**
 * @brief Write raw image data.
 *
 * If the image to write is only defined within a certain region of the
 * image domain as defined by the image header, the image is first
 * resized by this function such that the image data is defined on the
 * image domain specified by the image header before the image is written
 * to disk. Therefore, the ResizeImage() function is called with only the
 * image as argument. See also documentation of ResizeImage().
 *
 * If the image shall actually be written with a different size (image domain),
 * the @c subdomain argument of the ResizeImage() function has to be set to
 * @c false when resizing the image, or the image @c region attribute be set as
 * follows before writing the image to disk:
 * @code
 * image->region.ox = image->region.oy = image->region.oz = 0;
 * image->region.nx = image->hdr.dim[1];
 * image->region.ny = image->hdr.dim[2];
 * image->region.nz = image->hdr.dim[3];
 * @endcode
 *
 * @param [in] filename Name of image file.
 * @param [in] image    Image to write.
 * @param [in] mask     Mask image. Zero is written at voxels outside mask.
 *
 * @returns Whether the image was written to file.
 *
 * @sa WriteNiftiImage()
 * @sa WriteMetaImage()
 */
bool WriteRawImage(const char* filename, const Image* image, const Image* mask = NULL);

// ===========================================================================
// transformation - either affine or deformable
// ===========================================================================

/**
 * @brief Read affine or deformable transformation from file.
 *
 * @param [in]  filename          Name of transformation file.
 * @param [out] deformation_field Read deformation field or NULL if input file
 *                                could not be read, the memory allocation failed,
 *                                or the input file is an affine transformation instead.
 * @param [out] transform         Read affine transformation or an invalid
 *                                transformation, i.e., all elements zero,
 *                                if input file could not be read or the input
 *                                file is a deformation field instead.
 * @param [in]  fmt               Format of the transformation. In case of
 *                                DRAMMS format, the first two dimensions are
 *                                swapped after reading the transformation.
 *
 * @returns Whether either @p deformation_field or @p transform was read.
 */
bool ReadTransform(const char*       filename,
                   Image*&           deformation_field,
                   Image::Transform& transform,
                   Image::Format     fmt = Image::FORMAT_DRAMMS);

// ===========================================================================
// ANALYZE 7.5 or NIfTI-1 image
// ===========================================================================

/**
 * @brief Check whether a given image file is in the ANALYZE 7.5 or NIfTI-1 format.
 *
 * This is a slightly modified version of the is_nifti_file() function of
 * the nifticlib library, where only the output parameter for the read header
 * was added. Moreover, no debug messages are printed by this function.
 *
 * @param [in]  filename Image file name optionally without extension.
 * @param [out] hdr      Read NIfTI-1 or ANALYZE 7.5 header (in NIfTI-1 format) if not NULL.
 *
 * @retval  0 if file looks like ANALYZE 7.5 [checks sizeof_hdr field == 348].
 * @retval  1 if file marked as NIfTI-1 (header+data in 1 file).
 * @retval  2 if file marked as NIfTI-1 (header+data in 2 files).
 * @retval -1 if it cannot tell, file does not exist, etc.
 */
int IsNiftiImage(const char* filename, nifti_1_header *hdr = NULL);

/**
 * @brief Check whether a given image file is a deformation field in the
 *        ANALYZE 7.5 or NIfTI-1 format.
 *
 * This is a slightly modified version of the is_nifti_file() function of
 * the nifticlib library, where only the output parameter for the read header
 * was added. Moreover, no debug messages are printed by this function.
 *
 * @param [in]  filename Image file name optionally without extension.
 * @param [out] hdr      Read NIfTI-1 or ANALYZE 7.5 header (in NIfTI format) if not NULL.
 *
 * @retval  0 if file is a deformation field and looks like ANALYZE 7.5 [checks sizeof_hdr field == 348].
 * @retval  1 if file is a deformation field marked as NIfTI-1 (header+data in 1 file).
 * @retval  2 if file is a deformation field marked as NIfTI-1 (header+data in 2 files).
 * @retval -1 if image is not a deformation field, the format is not known,
 *            file does not exist, etc.
 */
int IsNiftiDeformationField(const char* filename, nifti_1_header *hdr = NULL);

/**
 * @brief Check whether a given image file can be read by ReadNiftiImage().
 *
 * @param [in] filename Name of image file.
 *
 * @returns Whether ReadNiftiImage() can read the given image file.
 */
bool CanReadNiftiImage(const char* filename);

/**
 * @brief Read NIfTI-1 or ANALYZE 7.5 image.
 *
 * @param [in] filename Name of image file.
 * @param [in] fmt      Image format.
 *
 * @returns Image of datatype as defined by the NIfTI-1 header or NULL on error.
 *          The returned image has to be destroyed using the @c delete operator.
 *
 * @sa ReadNiftiImage(const char*, short, bool, Image::Format)
 */
Image* ReadNiftiImage(const char* filename, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Read NIfTI-1 or ANALYZE 7.5 image.
 *
 * See documentation of the CastImage() function for details on the rescaling
 * of intensity values that is optionally applied by this function when reading
 * the input image.
 *
 * @param [in] filename Name of image file.
 * @param [in] datatype NIfTI-1 datatype code, e.g., DT_UNSIGNED_CHAR or
 *                      NIFTI_TYPE_UINT8. The image data, which may be of
 *                      other datatype, is casted to the desired output
 *                      datatype using CastImage().
 * @param [in] scale    See @p scale parameter of CastImage() function.
 * @param [in] fmt      Image format.
 *
 * @returns Image or NULL on error. The returned image has to be destroyed
 *          by the caller using the @c delete operator.
 *
 * @sa ReadNiftiImage(const char*, Image::Format)
 * @sa CastImage(const Image*, short, bool)
 */
Image* ReadNiftiImage(const char*   filename,
                      short         datatype,
                      bool          scale = false,
                      Image::Format fmt   = Image::FORMAT_DRAMMS);

/**
 * @brief Read NIfTI-1 or ANALYZE 7.5 image sequence.
 *
 * @param [in] filename Name of image file.
 * @param [in] fmt      Image format.
 *
 * @returns Image of datatype as defined by the NIfTI-1 header or NULL on error.
 *          The returned images have to be destroyed using the @c delete operator.
 *
 * @sa ReadNiftiSequence(const char*, short, bool, Image::Format)
 */
Sequence ReadNiftiSequence(const char* filename, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Read NIfTI-1 or ANALYZE 7.5 image.
 *
 * See documentation of the CastImage() function for details on the rescaling
 * of intensity values that is optionally applied by this function when reading
 * the input image.
 *
 * @param [in] filename Name of image file.
 * @param [in] datatype NIfTI-1 datatype code, e.g., DT_UNSIGNED_CHAR or
 *                      NIFTI_TYPE_UINT8. The image data, which may be of
 *                      other datatype, is casted to the desired output
 *                      datatype using CastImage().
 * @param [in] scale    See @p scale parameter of CastImage() function.
 * @param [in] fmt      Image format.
 *
 * @returns Image or NULL on error. The returned image has to be destroyed
 *          by the caller using the @c delete operator.
 *
 * @sa ReadNiftiSequence(const char*, Image::Format)
 * @sa CastImage(const Image*, short, bool)
 */
Sequence ReadNiftiSequence(const char*   filename,
                           short         datatype,
                           bool          scale = false,
                           Image::Format fmt   = Image::FORMAT_DRAMMS);

/**
 * @brief Read deformation field from NIfTI-1 file.
 *
 * @param [in] filename Name of image file.
 * @param [in] fmt      Image format.
 *
 * @returns Deformation field or NULL on error.
 */
Image* ReadNiftiDeformationField(const char*   filename,
                                 Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Write ANALYZE 7.5 or NIfTI-1 image.
 *
 * @param [in] filename Name of image file.
 * @param [in] image    Image to write.
 * @param [in] type     NIfTI-1 type of image file.
 *
 * @returns Whether the image was written to file.
 *
 * @sa WriteNiftiImage(const char*, const Image*, const Image*, Image::NiftiType)
 */
bool WriteNiftiImage(const char*      filename,
                     const Image*     image,
                     Image::NiftiType type = Image::NIFTI_TYPE_UNKNOWN);

/**
 * @brief Write ANALYZE 7.5 or NIfTI-1 image.
 *
 * If the image to write is only defined within a certain region of the
 * image domain as defined by the image header, the image is first
 * resized by this function such that the image data is defined on the
 * image domain specified by the image header before the image is written
 * to disk. Therefore, the ResizeImage() function is called with only the
 * image as argument. See also documentation of ResizeImage().
 *
 * If the image shall actually be written with a different size (image domain),
 * the @c subdomain argument of the ResizeImage() function has to be set to
 * @c false when resizing the image, or the image @c region attribute be set as
 * follows before writing the image to disk:
 * @code
 * image->region.ox = image->region.oy = image->region.oz = 0;
 * image->region.nx = image->hdr.dim[1];
 * image->region.ny = image->hdr.dim[2];
 * image->region.nz = image->hdr.dim[3];
 * @endcode
 *
 * @param [in] filename Name of image file.
 * @param [in] image    Image to write.
 * @param [in] mask     Mask image. Zero is written at voxels outside mask.
 * @param [in] type     NIfTI-1 type of image file. If set to
 *                      Image::NIFTI_TYPE_UNKNOWN, the default, the NIfTI-1 type
 *                      is derived first from the specified file name extension
 *                      if provided and secondly from the @c filefmt attribute
 *                      of the @p image. Otherwise, the specified NIfTI-1 type
 *                      is used instead.
 *
 * @returns Whether the image was written to file.
 */
bool WriteNiftiImage(const char*      filename,
                     const Image*     image,
                     const Image*     mask,
                     Image::NiftiType type = Image::NIFTI_TYPE_UNKNOWN);

/**
 * @brief Write ANALYZE 7.5 or NIfTI-1 image sequence.
 *
 * @param [in] filename Name of image file.
 * @param [in] images   Image sequence to write.
 * @param [in] type     NIfTI-1 type of image file.
 *
 * @returns Whether the image sequence was written to file.
 *
 * @sa WriteNiftiSequence(const char*, const ConstSequence, const Image*, Image::NiftiType)
 */
bool WriteNiftiSequence(const char*         filename,
                        const ConstSequence images,
                        Image::NiftiType    type = Image::NIFTI_TYPE_UNKNOWN);

/**
 * @brief Write ANALYZE 7.5 or NIfTI-1 image sequence.
 *
 * If an image to write is only defined within a certain region of the
 * image domain as defined by the image header, the image is first
 * resized by this function such that the image data is defined on the
 * image domain specified by the image header before the image is written
 * to disk. Therefore, the ResizeImage() function is called with only the
 * image as argument. See also documentation of ResizeImage().
 *
 * If the image shall actually be written with a different size (image domain),
 * the @c subdomain argument of the ResizeImage() function has to be set to
 * @c false when resizing the image, or the image @c region attribute be set as
 * follows before writing the image to disk:
 * @code
 * image->region.ox = image->region.oy = image->region.oz = 0;
 * image->region.nx = image->hdr.dim[1];
 * image->region.ny = image->hdr.dim[2];
 * image->region.nz = image->hdr.dim[3];
 * @endcode
 *
 * @param [in] filename Name of image file.
 * @param [in] images   Image sequence to write.
 * @param [in] mask     Mask image. Zero is written at voxels outside mask.
 * @param [in] type     NIfTI-1 type of image file. If set to
 *                      Image::NIFTI_TYPE_UNKNOWN, the default, the NIfTI-1 type
 *                      is derived first from the specified file name extension
 *                      if provided and secondly from the @c filefmt attribute
 *                      of the @p image. Otherwise, the specified NIfTI-1 type
 *                      is used instead.
 *
 * @returns Whether the image sequence was written to file.
 */
bool WriteNiftiSequence(const char*         filename,
                        const ConstSequence images,
                        const Image*        mask,
                        Image::NiftiType    type = Image::NIFTI_TYPE_UNKNOWN);

// ===========================================================================
// Meta image
// ===========================================================================

/**
 * @brief Check whether a given image file is in the Meta image format.
 *
 * @param [in]  filename Image file name optionally without extension.
 * @param [out] hdr      Read Meta image header in NIfTI format if not NULL.
 *
 * @retval  1 if file is a Meta image.
 * @retval -1 if it cannot tell, file does not exist, etc.
 */
int IsMetaImage(const char* filename, nifti_1_header* hdr = NULL);

/**
 * @brief Check whether a given image file is a deformation field in the Meta image format.
 *
 * @param [in]  filename Image file name optionally without extension.
 * @param [out] hdr      Read Meta image header in NIfTI format if not NULL.
 *
 * @retval  1 if file is a Meta image representing a deformation field.
 * @retval -1 if file is no deformation field, it cannot tell, file does not exist, etc.
 */
int IsMetaDeformationField(const char* filename, nifti_1_header* hdr = NULL);

/**
 * @brief Check whether a given image file can be read by ReadMetaImage().
 *
 * @param [in] filename Name of image file.
 *
 * @returns Whether ReadMetaImage() can read the given Meta image file.
 */
bool CanReadMetaImage(const char* filename);

/**
 * @brief Read Meta image from file.
 *
 * @param [in] filename Name of image file.
 * @param [in] fmt      Image format.
 *
 * @returns Image of datatype as defined by the image header or NULL on error.
 *          The returned image has to be destroyed using the @c delete operator.
 *
 * @sa ReadMetaImage(const char*, short, bool, Image::Format)
 */
Image* ReadMetaImage(const char* filename, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Read Meta image from file.
 *
 * See documentation of the CastImage() function for details on the rescaling
 * of intensity values that is optionally applied by this function when reading
 * the input image.
 *
 * @param [in] filename Name of image file.
 * @param [in] datatype NIfTI-1 datatype code, e.g., DT_UNSIGNED_CHAR or
 *                      NIFTI_TYPE_UINT8. The image data, which may be of
 *                      other datatype, is casted to the desired output
 *                      datatype using CastImage().
 * @param [in] scale    See @p scale parameter of CastImage() function.
 * @param [in] fmt      Image format.
 *
 * @returns Image or NULL on error. The returned image has to be destroyed
 *          by the caller using the @c delete operator.
 *
 * @sa ReadMetaImage(const char*, Image::Format)
 * @sa CastImage(const Image*, short, bool)
 */
Image* ReadMetaImage(const char*   filename,
                     short         datatype,
                     bool          scale = false,
                     Image::Format fmt   = Image::FORMAT_DRAMMS);

/**
 * @brief Read deformation field from Meta image file.
 *
 * @param [in] filename Name of image file.
 * @param [in] fmt      Image format.
 *
 * @returns Deformation field or NULL on error.
 */
Image* ReadMetaDeformationField(const char*   filename,
                                Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Write Meta image.
 *
 * If the image to write is only defined within a certain region of the
 * image domain as defined by the image header, the image is first
 * resized by this function such that the image data is defined on the
 * image domain specified by the image header before the image is written
 * to disk. Therefore, the ResizeImage() function is called with only the
 * image as argument. See also documentation of ResizeImage().
 *
 * If the image shall actually be written with a different size (image domain),
 * the @c subdomain argument of the ResizeImage() function has to be set to
 * @c false when resizing the image, or the image @c region attribute be set as
 * follows before writing the image to disk:
 * @code
 * image->region.ox = image->region.oy = image->region.oz = 0;
 * image->region.nx = image->hdr.dim[1];
 * image->region.ny = image->hdr.dim[2];
 * image->region.nz = image->hdr.dim[3];
 * @endcode
 *
 * @param [in] filename Name of image file.
 * @param [in] image    Image to write.
 * @param [in] mask     Mask image. Zero is written at voxels outside mask.
 *
 * @returns Whether the image was written to file.
 *
 * @sa WriteRawImage()
 */
bool WriteMetaImage(const char* filename, const Image* image, const Image* mask = NULL);

// ===========================================================================
// affine transformation
// ===========================================================================

/**
 * @brief Read affine transformation matrix from file.
 *
 * @param [in]  filename  Name of transformation file.
 * @param [out] transform Affine transformation matrix.
 * @param [in] fmt       Format of the affine transformation matrix. In case of
 *                       DRAMMS, the first two columns and rows of the matrix
 *                       are swapped after reading the matrix from file.
 *
 * @retval -1 Failed to open file for reading.
 * @retval  0 Invalid transformation file.
 * @retval  1 Successfully read affine transformation.
 */
int ReadAffineTransform(const char* filename, Image::Transform& transform, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Write affine transformation matrix to file.
 *
 * @param [in] filename  Name of transformation file.
 * @param [in] transform Affine transformation matrix.
 * @param [in] fmt       Format of the affine transformation matrix. In case of
 *                       DRAMMS, the first two columns and rows of the matrix
 *                       are swapped back before writing the matrix to file.
 *
 * @returns Whether the affine transformation was written to the specified file.
 */
bool WriteAffineTransform(const char* filename, Image::Transform transform, Image::Format fmt = Image::FORMAT_DRAMMS);


} // namespace dramms


#endif // _DRAMMS_COMMON_IMAGEIO_H
