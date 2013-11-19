/**
 * @file  image.h
 * @brief Definition of image type and declaration of image functions.
 *
 * Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
 * See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#pragma once
#ifndef _SBIA_DRAMMS_COMMON_IMAGE_H
#define _SBIA_DRAMMS_COMMON_IMAGE_H


#include <iostream>
#include <vector>
#include <nifti1_io.h> // Image::Header, mat44

#include "mvcd.h" // Fvector3d


namespace dramms {


/**
 * @class Image
 * @brief Data structure used to represent an image.
 *
 * The layout of the image data in @c img is identical for scalar images,
 * for @c FORMAT_ITK, @c FORMAT_FSL, and @c FORMAT_DRAMMS, i.e., given voxel
 * coordinates (i, j, k) and orienation RAS, i runs from left to right,
 * j from posterior to anterior, and k from inferior to superior, the data
 * associated with this voxel position is stored in img.uc[k][j][i]
 * (if the image datatype is @c DT_UNSIGNED_CHAR). However, the meaning of
 * i and j (or x and y, respectively) is exchanged in the associated NIfTI
 * header in case of @c FORMAT_DRAMMS, i.e., the first dimension refers to
 * the y axis and the second dimension to the x axis. This allows code such
 * as the following to access the image data:
 * @code
 * Image image;
 * for (int k = 0; k < image.hdr.dim[3]; k++) {
 *     for (int i = 0; i < image.hdr.dim[1]; i++) {
 *         for (int j = 0; j < image.hdr.dim[2]; j++) {
 *             image.img.uc[k][i][j] = 0;
 *         }
 *     }
 * }
 * @endcode
 * Using the function PermuteXY(const Image::Header&), the format of the
 * NIfTI-1 header can be converted from one to the other.
 *
 * @note Looking at above code and given the swap of the header information
 *       related to the x and y axes, it appears as if the image axes were
 *       permuted to yxz. In fact, however, the image data in memory yet
 *       corresponds to the usual xyz order, only the header information is
 *       swapped.
 *
 * If the datatype of the image is not known, the set() and get() methods can
 * be used to set or get the image data of a scalar image as float. The
 * templated value() method, on the other side, can be used to access the image
 * data given any other type of the return value. Note that if the image data
 * is actually stored as a different datatype, a static cast is performed.
 * Moreover, these methods take the image format into consideration when
 * accessing the image data, i.e., in case of @c FORMAT_DRAMMS, the i (or x)
 * argument refers to the j (or y) coordinate of the @c FORMAT_ITK, and
 * the j (or y) argument to the i (or x) coordinate, respectively. Thus,
 * the following two loops result in accessing the same voxels in each
 * iteration:
 * @code
 * if (image.fmt == FORMAT_DRAMMS) {
 *     for (int k = 0; k < image.hdr.dim[3]; k++) {
 *         for (int i = 0; i < image.hdr.dim[1]; i++) {
 *             for (int j = 0; j < image.hdr.dim[2]; j++) {
 *                 float I = image.get(i, j, k);
 *             }
 *         }
 *     }
 * } else {
 *     for (int k = 0; k < image.hdr.dim[3]; k++) {
 *         for (int j = 0; j < image.hdr.dim[2]; j++) {
 *             for (int i = 0; i < image.hdr.dim[1]; i++) {
 *                 float I = image.get(i, j, k);
 *             }
 *         }
 *     }
 * }
 * @endcode
 *
 * As the image data in the memory referred to by @c img can be optionally
 * cropped to a smaller extent as the one specified by the image header,
 * the image @c region should be used instead of the information about the
 * original image size in @c hdr.dim, i.e.,
 * @code
 * if (image.fmt == FORMAT_DRAMMS) {
 *     for (int k = 0; k < image.region.nz; k++) {
 *         for (int i = 0; i < image.region.nx; i++) {
 *             for (int j = 0; j < image.region.ny; j++) {
 *                 unsigned char I = image.img.uc[k][i][j];
 *             }
 *         }
 *     }
 * } else {
 *     for (int k = 0; k < image.region.nz; k++) {
 *         for (int j = 0; j < image.region.ny; j++) {
 *             for (int i = 0; i < image.region.nx; i++) {
 *                 unsigned char I = image.img.uc[k][j][i];
 *             }
 *         }
 *     }
 * }
 * @endcode
 * When writing the image to file, the image data has to be padded with zeros
 * (or any other reasonable value) if the image region is smaller than the
 * image size specified by the image header.
 */
class Image
{
    // -----------------------------------------------------------------------
    // image related types
public:

    /**
     * @brief Represents an image region in three dimensions.
     *
     * An image region is defined by an offset in voxel units and the size of the
     * image region. It is always relative to the entire image volume as specified
     * by the image's header information.
     *
     * Example:
     * @code
     * Image::Region region;
     * @endcode
     *
     * Constraints on image region:
     * @code
     * 0 <= region.ox <= image.region.nx
     * 0 <= region.oy <= image.region.ny
     * 0 <= region.oz <= image.region.nz
     * @endcode
     *
     * @note If these constraints on the image region are not met,
     *       the image region is considered invalid. Functions
     *       can return an invalid region to indicate an error.
     *
     * Image region of entire image:
     * @code
     * region.ox = 0;
     * region.oy = 0;
     * region.oz = 0;
     * region.nx = image.hdr.dim[1];
     * region.ny = image.hdr.dim[2];
     * region.nz = image.hdr.dim[3];
     * @endcode
     */
    struct Region
    {
        int ox; ///< Offset of image region in x direction.
        int oy; ///< Offset of image region in y direction.
        int oz; ///< Offset of image region in z direction.
        int nx; ///< Size of image region in x direction.
        int ny; ///< Size of image region in y direction.
        int nz; ///< Size of image region in z direction.

        Region() : ox(0), oy(0), oz(0), nx(0), ny(0), nz(0) {}

        bool operator==(const Region& rhs)
        {
            return ox == rhs.ox && oy == rhs.oy && oz == rhs.oz &&
                   nx == rhs.nx && ny == rhs.ny && nz == rhs.nz;
        }

        bool operator!=(const Region& rhs)
        {
            return !(*this == rhs);
        }
    };

    /**
     * @brief Type used for pointer to image data.
     */
    union Data
    {
        unsigned char***  uc; // unsigned char
        short***          ss; // signed short
        unsigned short*** us; // unsigned short
        int***            si; // signed integer
        float***          fl; // float
        Fvector2d***      v2; // 2D float vector
        Fvector3d***      v3; // 3D float vector

        Data() : uc(NULL) {}
        Data(unsigned char***  p) : uc(p) {}
        Data(short***          p) : ss(p) {}
        Data(unsigned short*** p) : us(p) {}
        Data(int***            p) : si(p) {}
        Data(float***          p) : fl(p) {}
        Data(Fvector2d***      p) : v2(p) {}
        Data(Fvector3d***      p) : v3(p) {}
    };

    /**
     * @brief Format of image representation.
     */
    enum Format {
        /**
         * @brief Unknown format.
         */
        FORMAT_UNKNOWN,

        /**
         * @brief Format used by ITK.
         *
         * In case of a vector field, the displacement vector components are
         * given in physical units and the vector components at each voxel are
         * stored consecutively, i.e., xyzxyzxyz...
         */
        FORMAT_ITK,

        /**
         * @brief Image format used by FSL.
         *
         * The format of FSL is similar to the @c FORMAT_ITK. The only
         * difference is the storage of vector fields such as the deformation
         * field output by FNIRT, in particular. Here, the vector component
         * images are concatenated which results in an image which has the
         * deformation field components stored as xxx...yyy...zzz...
         *
         * @note As long as the image is represented in memory only, i.e.,
         *       by an instance of Image, this format is identical to
         *       FORMAT_ITK. Otherwise, accessing the vectors of a deformation
         *       field would be very inefficient and confusing. Only when
         *       an image is read from disk or written to disk, the
         *       rearrangement of the vector components is applied.
         *
         * @sa http://www.fmrib.ox.ac.uk/fsl/fnirt/index.html
         */
        FORMAT_FSL,

        /**
         * @brief Image format used by DRAMMS (and initially DRAMMS as well).
         *
         * In this format, the x and y axis are exchanged, without, however,
         * also exchanging the image data (!), i.e., only the information in
         * the image header is changed such that it appears as if the x and y
         * axis where permuted. It is important to note that this permutation
         * only is applied to the header information, not the image data itself.
         *
         * To iterate over an image given in this format, the following code
         * is usually used:
         * @code
         * for (int k = 0; k < image->region.nz; k++) {
         *     for (int i = 0; i < image->region.nx; i++) {
         *         for (int j = 0; j < image->region.ny; j++) {
         *             unsigned char I = image->img.uc[k][i][j];
         *         }
         *     }
         * }
         * @endcode
         *
         * For vector fields, on the other side, it is important to not that
         * the first two vector components are as a result of the different
         * interpretation of the first two image dimensions permuted. Hence,
         * in order to convert to the standard format used by most image file
         * formats and imaging libraries, the information in the header has to
         * be adjusted (see PermuteXY(const Image::Header&)) and in case of a
         * vector field, the first two components of each vector have to be
         * swapped.
         *
         * Moreover, displacements are given in voxel units instead of physical
         * units, i.e., the displacmente has to be multiplied by the voxel size
         * (pixdim field of NIfTI-1 header) in order to get the corresponding
         * displacement of the spatial coordinates.
         */
        FORMAT_DRAMMS
    };

    /**
     * @brief Type used to represent 4x4 coordinate transformation matrices.
     */
    typedef mat44 Transform;

    /**
     * @brief Left-right order of image data.
     */
    enum LROrder {
        LR_UNKNOWN      =  0, ///< left-right order unknown
        LR_RADIOLOGICAL = -1, ///< radiological order as in ANALYZE, det(s/qform) < 0
        LR_NEUROLOGICAL =  1  ///< neurological order as in NIfTI, det(s/qform) >= 0
    };

    /**
     * @brief Type used to represent image header.
     */
    typedef nifti_1_header Header;

    /**
     * @brief File format of NIfTI or ANALYZE image.
     *
     * Note that the enumeration values DO NOT correspond to the integer values used
     * for the nifti_type field of the nifti_image data structure as defined in
     * the nifti1_io.h header file. Instead, they are offset by 1.
     */
    enum NiftiType {
        NIFTI_TYPE_UNKNOWN,                            ///< Unknown NIfTI-1 type.
        NIFTI_TYPE_ANALYZE = NIFTI_FTYPE_ANALYZE  + 1, ///< ANALYZE 7.5 .hdr(.gz)/.img(.gz).
        NIFTI_TYPE_NIFTI_1 = NIFTI_FTYPE_NIFTI1_1 + 1, ///< NIfTI-1 .nii(.gz) file.
        NIFTI_TYPE_NIFTI_2 = NIFTI_FTYPE_NIFTI1_2 + 1, ///< NIfTI-1 .hdr(.gz)/.img(.gz) image pair.
        NIFTI_TYPE_NIFTI_A = NIFTI_FTYPE_ASCII    + 1  ///< NIfTI-1 .nia(.gz) file.
    };

    /**
     * @brief Format of file from which image has been read or in which image
     *        should be written to disk, respectively.
     *
     * Whether the image data was compressed or shall be compressed, respectively,
     * is controlled by the @c compress attribute of the Image.
     *
     * Note that the mapping of NiftiType to FileFormat is 1-to-1 for all defined
     * NIfTI-1 types. The FileFormat enumeration, however, contains some more values.
     */
    enum FileFormat {
        FILE_FORMAT_UNKNOWN = NIFTI_TYPE_UNKNOWN, ///< Unknown file format.
        FILE_FORMAT_ANALYZE = NIFTI_TYPE_ANALYZE, ///< ANALYZE 7.5 .hdr(.gz)/.img(.gz).
        FILE_FORMAT_NIFTI_1 = NIFTI_TYPE_NIFTI_1, ///< NIfTI-1 .nii(.gz).
        FILE_FORMAT_NIFTI_2 = NIFTI_TYPE_NIFTI_2, ///< NIfTI-1 .hdr(.gz)/.img(.gz).
        FILE_FORMAT_NIFTI_A = NIFTI_TYPE_NIFTI_A, ///< NIfTI-1 .nia(.gz).
        FILE_FORMAT_META                          ///< Meta image .mhd/.raw(.gz).
    };

    // -----------------------------------------------------------------------
    // construction / destruction
public:

    /**
     * @brief Default constructor.
     *
     * This constructor does not allocate any image data and initializes the
     * NIfTI header to zero.
     */
    Image();

    /**
     * @brief Construct new image and allocate image buffer.
     *
     * Constructs a new image with voxel size 1x1x1 mm^3, image origin at
     * (0, 0, 0), and orientation LPS, i.e., transverse flipped.
     *
     * @param [in] x_dim    Size along x axis.
     * @param [in] y_dim    Size along y axis.
     * @param [in] z_dim    Size along z axis.
     * @param [in] datatype NIfTI-1 datatype code.
     * @param [in] n        Number of components per voxel.
     * @param [in] fmt      Image format.
     *
     * @throws std::invalid_argument if any of the input arguments has an
     *                               invalid value such as a dimension of zero,
     *                               for example.
     */
    Image(int    x_dim,
          int    y_dim,
          int    z_dim,
          short  datatype = DT_UNSIGNED_CHAR,
          int    n        = 1,
          Format fmt      = Image::FORMAT_DRAMMS);

    /**
     * @brief Construct new image and allocate image buffer.
     *
     * Constructs a new image with image origin at (0, 0, 0), and orientation
     * LPS, i.e., transverse flipped.
     *
     * @param [in] x_dim    Size along x axis.
     * @param [in] y_dim    Size along y axis.
     * @param [in] z_dim    Size along z axis.
     * @param [in] x_pixdim Voxel size along x axis.
     * @param [in] y_pixdim Voxel size along y axis.
     * @param [in] z_pixdim Voxel size along z axis.
     * @param [in] datatype NIfTI-1 datatype code.
     * @param [in] n        Number of components per voxel.
     * @param [in] fmt      Image format.
     *
     * @throws std::invalid_argument if any of the input arguments has an
     *                               invalid value such as a dimension of zero,
     *                               for example.
     */
    Image(int    x_dim,
          int    y_dim,
          int    z_dim,
          float  x_pixdim,
          float  y_pixdim,
          float  z_pixdim,
          short  datatype = DT_UNSIGNED_CHAR,
          int    n        = 1,
          Format fmt      = Image::FORMAT_DRAMMS);

    /**
     * @brief Construct new image given already allocated image buffer.
     *
     * Constructs a new image using the already pre-allocated image buffer
     * with voxel size 1x1x1 mm^3, image origin at (0, 0, 0), and orientation
     * LPS, i.e., transverse flipped.
     *
     * @param [in] img      Pre-allocated image buffer.
     * @param [in] x_dim    Size along x axis.
     * @param [in] y_dim    Size along y axis.
     * @param [in] z_dim    Size along z axis.
     * @param [in] datatype NIfTI-1 datatype code.
     * @param [in] n        Number of components per voxel.
     * @param [in] fmt      Image format.
     * @param [in] owns_img Whether the newly constructed image shall take
     *                      over ownership of the image data. If true, this
     *                      image instance will delete the image data upon
     *                      construction.
     *
     * @throws std::invalid_argument if any of the input arguments has an
     *                               invalid value such as a dimension of zero,
     *                               for example.
     */
    Image(Data   img,
          int    x_dim,
          int    y_dim,
          int    z_dim,
          short  datatype = DT_UNSIGNED_CHAR,
          int    n        = 1,
          Format fmt      = Image::FORMAT_DRAMMS,
          bool   owns_img = true);

    /**
     * @brief Construct new image given already allocated image buffer.
     *
     * Constructs a new image using the already pre-allocated image buffer
     * with voxel size 1x1x1 mm^3, image origin at (0, 0, 0), and orientation
     * LPS, i.e., transverse flipped.
     *
     * @param [in] img      Pre-allocated image buffer.
     * @param [in] x_dim    Size along x axis.
     * @param [in] y_dim    Size along y axis.
     * @param [in] z_dim    Size along z axis.
     * @param [in] x_pixdim Voxel size along x axis.
     * @param [in] y_pixdim Voxel size along y axis.
     * @param [in] z_pixdim Voxel size along z axis.
     * @param [in] datatype NIfTI-1 datatype code.
     * @param [in] n        Number of components per voxel.
     * @param [in] fmt      Image format.
     * @param [in] owns_img Whether the newly constructed image shall take
     *                      over ownership of the image data. If true, this
     *                      image instance will delete the image data upon
     *                      construction.
     *
     * @throws std::invalid_argument if any of the input arguments has an
     *                               invalid value such as a dimension of zero,
     *                               for example.
     */
    Image(Data   img,
          int    x_dim,
          int    y_dim,
          int    z_dim,
          float  x_pixdim,
          float  y_pixdim,
          float  z_pixdim,
          short  datatype = DT_UNSIGNED_CHAR,
          int    n        = 1,
          Format fmt      = Image::FORMAT_DRAMMS,
          bool   owns_img = true);

    /**
     * @brief Destructor.
     *
     * Frees the image buffer if this Image instance is the owner of this
     * buffer, i.e., the attribute owns_img is true.
     */
    ~Image();

    // -----------------------------------------------------------------------
    // set/get voxel data
public:

    /**
     * @brief Get image intensity at voxel position.
     *
     * @code
     * Image image;
     * image.value<unsigned char>(i, j, k);
     * image.value<float>(i, j, k);
     * @endcode
     *
     * @param [in] i Voxel index in first dimension.
     * @param [in] j Voxel index in second dimension.
     * @param [in] k Voxel index in third dimension.
     *               Set to zero if image is two-dimensional.
     * @param [in] n Vector component. Set to zero if voxel data is scalar.
     *
     * @returns Voxel value casted to specified template type.
     */
    template <typename TValue>
    TValue value(int i, int j, int k = 0, int n = 0) const;

    /**
     * @brief Get image intensity at voxel position.
     *
     * @param [in] i Voxel index in first dimension.
     * @param [in] j Voxel index in second dimension.
     * @param [in] k Voxel index in third dimension.
     *               Set to zero if image is two-dimensional.
     * @param [in] n Vector component. Set to zero if voxel data is a scalar.
     *
     * @returns Voxel value casted to float.
     */
    float get(int i, int j, int k = 0, int n = 0) const;

    /**
     * @brief Get image components at voxel position.
     *
     * @param [in]  i Voxel index in first dimension.
     * @param [in]  j Voxel index in second dimension.
     * @param [out] v Vector components casted to TValue.
     */
    template <typename TValue>
    void get(int i, int j, std::vector<TValue>& v) const;

    /**
     * @brief Get image components at voxel position.
     *
     * @param [in]  i Voxel index in first dimension.
     * @param [in]  j Voxel index in second dimension.
     * @param [in]  k Voxel index in third dimension.
     * @param [out] v Vector components casted to TValue.
     */
    template <typename TValue>
    void get(int i, int j, int k, std::vector<TValue>& v) const;

    /**
     * @brief Get interpolated image components at sub-voxel position.
     *
     * This method performs a bilinear interpolation for each vector component.
     *
     * @param [in]  i Continuous voxel index in first dimension.
     * @param [in]  j Continuous voxel index in second dimension.
     * @param [out] v Interpolated vector components.
     */
    template <typename TValue>
    void get(float i, float j, std::vector<TValue>& v) const;

    /**
     * @brief Get interpolated image components at sub-voxel position.
     *
     * This method performs a bilinear interpolation for each vector component
     * within the specified slice of a three-dimensional volume.
     *
     * @param [in]  i Continous voxel index in first dimension.
     * @param [in]  j Continous voxel index in second dimension.
     * @param [in]  k Slice index. Set to 0 for two-dimensional images.
     * @param [out] v Vector components casted to TValue.
     */
    template <typename TValue>
    void get(float i, float j, int k, std::vector<TValue>& v) const;

    /**
     * @brief Get interpolated image components at sub-voxel position.
     *
     * This method performs trilinear interpolation for each vector component.
     *
     * @param [in]  i Continous voxel index in first dimension.
     * @param [in]  j Continous voxel index in second dimension.
     * @param [in]  k Continous voxel index in third dimension.
     * @param [out] v Vector components casted to TValue.
     */
    template <typename TValue>
    void get(float i, float j, float k, std::vector<TValue>& v) const;

    /**
     * @brief Set image intensity at given voxel position.
     *
     * @param [in] i Voxel index in first dimension.
     * @param [in] j Voxel index in second dimension.
     * @param [in] v Voxel intensity. Will be casted to image datatype.
     */
    void set(int i, int j, float v);

    /**
     * @brief Set image intensity at given voxel position.
     *
     * @param [in] i Voxel index in first dimension.
     * @param [in] j Voxel index in second dimension.
     * @param [in] k Voxel index in third dimension.
     * @param [in] v Voxel intensity. Will be casted to image datatype.
     */
    void set(int i, int j, int k, float v);

    /**
     * @brief Set image intensity at given voxel position.
     *
     * @param [in] i Voxel index in first dimension.
     * @param [in] j Voxel index in second dimension.
     * @param [in] k Voxel index in third dimension.
     * @param [in] n Vector component. Set to zero if voxel data is scalar.
     * @param [in] v Voxel intensity. Will be casted to image datatype.
     */
    void set(int i, int j, int k, int n, float v);

    /**
     * @brief Set image components at given voxel position.
     *
     * @param [in] i Voxel index in first dimension.
     * @param [in] j Voxel index in second dimension.
     * @param [in] v Vector components. Will be casted to image datatype.
     */
    template <typename TValue>
    void set(int i, int j, const std::vector<TValue>& v);

    /**
     * @brief Set image components at given voxel position.
     *
     * @param [in] i Voxel index in first dimension.
     * @param [in] j Voxel index in second dimension.
     * @param [in] k Voxel index in third dimension.
     * @param [in] v Vector components. Will be casted to image datatype.
     */
    template <typename TValue>
    void set(int i, int j, int k, const std::vector<TValue>& v);

    /**
     * @brief Get image intensity at sub-voxel position.
     *
     * This function interpolates the intensity using bilinear interpolation
     * of the intensities in the slice of the given @p z coordinate.
     *
     * @param [in] i Continuous voxel index in first dimension.
     * @param [in] j Continuous voxel index in second dimension.
     * @param [in] k Voxel index in third dimension.
     *               Set to zero if image is two-dimensional.
     * @param [in] n Vector component. Set to zero if voxel data is scalar.
     *
     * @returns Interpolated intensity value.
     */
    float value(float i, float j, int k = 0, int n = 0) const;

    /**
     * @brief Get image intensity at sub-voxel position.
     *
     * This function interpolates the intensity using trilinear interpolation.
     *
     * @param [in] i Continuous voxel index in first dimension.
     * @param [in] j Continuous voxel index in second dimension.
     * @param [in] k Continuous voxel index in third dimension.
     * @param [in] n Vector component. Set to zero if voxel data is scalar.
     *
     * @returns Interpolated intensity value.
     */
    float value(float i, float j, float k, int n = 0) const;

    /**
     * @brief Get displacement vector at voxel position.
     *
     * @param [in] i Voxel index in first dimension.
     * @param [in] j Voxel index in second dimension.
     * @param [in] k Voxel index in third dimension.
     *
     * @returns Displacement vector of deformation field.
     */
    const Fvector3d& v3(int i, int j, int k) const;

    // -----------------------------------------------------------------------
    // sub-image
public:

    /**
     * @brief Get Image instance representing a single slice of an image volume.
     *
     * @param [in] k    Slice index.
     * @param [in] copy Whether to copy the image data.
     *
     * @returns Two-dimensional Image instance representing the specified slice.
     *          The returned instance has to be destroyed by the caller using
     *          the delete operator. Note, however, that if @p copy is false,
     *          the image data itself remains in the possession of this image
     *          instance. The returned slice image should in this case only be
     *          used temporarily.
     */
    Image* GetSlice(int k, bool copy = true);

    /**
     * @brief Get Image instance representing a single slice of an image volume.
     *
     * @param [in] k    Slice index.
     * @param [in] copy Whether to copy the image data.
     *
     * @returns Two-dimensional Image instance representing the specified slice.
     *          The returned instance has to be destroyed by the caller using
     *          the delete operator. Note, however, that if @p copy is false,
     *          the image data itself remains in the possession of this image
     *          instance. The returned slice image should in this case only be
     *          used temporarily.
     */
    const Image* GetSlice(int k, bool copy = true) const;

    // -----------------------------------------------------------------------
    // format of image representation
public:

    /**
     * @brief Set image format.
     *
     * This method performs the necessary conversions such that this Image
     * instance has the specified format afterwards.
     *
     * @param [in] fmt Desired image format.
     *
     * @returns Reference to this image instance.
     */
    void SetFormat(const Image::Format fmt);

    // -----------------------------------------------------------------------
    // queries
public:

    /**
     * @brief Get number of image components, i.e., 1 for scalar fields and
     *        > 1 for vector fields.
     *
     * @returns Number of image components.
     */
    int GetNumberOfComponents() const;

    /**
     * @brief Query whether this image has the same size as another image.
     *
     * This method takes into consideration the possible different image format.
     *
     * @param [in] rhs Other image.
     *
     * @returns Whether this image has the same size as the other image.
     */
    bool HasSameSizeAs(const Image* rhs) const;

    /**
     * @brief Query whether this image has a smaller size than another image.
     *
     * This method takes into consideration the possible different image format.
     *
     * @param [in] rhs Other image.
     *
     * @returns Whether this image has smaller size as the other image, i.e.,
     *          the size of this image in each dimension is less than the
     *          one of the corresponding dimension of the other image.
     */
    bool HasSmallerSizeAs(const Image* rhs) const;

    /**
     * @brief Query whether this image has a smaller or same size than another image.
     *
     * This method takes into consideration the possible different image format.
     *
     * @param [in] rhs Other image.
     *
     * @returns Whether this image has smaller or same size as the other image,
     *          i.e., the size of this image in each dimension is less than or
     *          equal to the one of the corresponding dimension of the other image.
     */
    bool HasSmallerOrSameSizeAs(const Image* rhs) const;

    // -----------------------------------------------------------------------
    // coordinate transformation
public:

    /**
     * @brief Set image orientation.
     *
     * Valid orientation codes are NIFTI_R2L, NIFTI_L2R, NIFTI_A2P, NIFTI_P2A,
     * NIFTI_I2S, and NIFTI_S2I as defined in nifti1_io.h.
     * This method also sets the qform_code to @c NIFTI_XFORM_SCANNER_ANAT.
     *
     * @attention If the image is in DRAMMS format, the @p x_orient corresponds to
     *            the @p y_orient in NIfTI and vice versa!
     *
     * @param [in] x_orient Orientation along x axis.
     * @param [in] y_orient Orientation along y axis.
     * @param [in] z_orient Orientation along z axis.
     */
    void SetOrientation(int x_orient, int y_orient, int z_orient);

    /**
     * @brief Update coordinate transformations.
     *
     * After the offset and quaternion representation of the qform
     * transformation of the image header has been modified, this method must
     * be called to update the 4x4 matrices qto_xyz and qto_ijk which are
     * used by the methods to_xyz() and to_ijk().
     */
    void UpdateTransforms();

    /**
     * @brief Transform voxel coordinates to spatial coordinates.
     *
     * This method uses the qform transformation to transform the given
     * voxel coordinates in the image space to spatial coordinates.
     *
     * @param [in] ijk Voxel coordinates.
     *
     * @returns Spatial coordinates.
     *
     * @sa UpdateTransforms()
     */
    Fvector3d to_xyz(const Ivector3d& ijk) const;

    /**
     * @brief Transform voxel coordinates to spatial coordinates.
     *
     * @param [in] i Voxel coorindate along x axis.
     * @param [in] j Voxel coorindate along y axis.
     * @param [in] k Voxel coorindate along z axis.
     *
     * @returns Spatial coordinates.
     *
     * @sa UpdateTransforms()
     */
    Fvector3d to_xyz(int i, int j, int k) const;

    /**
     * @brief Transform voxel coordinates to spatial coordinates.
     *
     * This method uses the qform transformation to transform the given
     * voxel coordinates in the image space to spatial coordinates.
     *
     * @param [in] ijk Continuous voxel coordinates.
     *
     * @returns Spatial coordinates.
     *
     * @sa UpdateTransforms()
     */
    Fvector3d to_xyz(const Fvector3d& ijk) const;

    /**
     * @brief Transform voxel coordinates to spatial coordinates.
     *
     * @param [in] i Continuous voxel coorindate along x axis.
     * @param [in] j Continuous voxel coorindate along y axis.
     * @param [in] k Continuous voxel coorindate along z axis.
     *
     * @returns Spatial coordinates.
     *
     * @sa UpdateTransforms()
     */
    Fvector3d to_xyz(float i, float j, float k) const;

    /**
     * @brief Transform spatial coordinates to continuous voxel coordinates.
     *
     * This method uses the inverse qform transformation to transform the given
     * spatial coordinates to voxel coordinates in the image space.
     *
     * @param [in] xyz Spatial coordinates.
     *
     * @returns Continuous voxel coordinates.
     *
     * @sa UpdateTransforms()
     */
    Fvector3d to_cijk(const Fvector3d& xyz) const;

    /**
     * @brief Transform spatial coordinates to continuous voxel coordinates.
     *
     * @param [in] x Spatial coorindate along x axis.
     * @param [in] y Spatial coorindate along y axis.
     * @param [in] z Spatial coorindate along z axis.
     *
     * @returns Continuous voxel coordinates.
     *
     * @sa UpdateTransforms()
     */
    Fvector3d to_cijk(float x, float y, float z) const;

    /**
     * @brief Transform spatial coordinates to voxel coordinates.
     *
     * This method uses the inverse qform transformation to transform the given
     * spatial coordinates to voxel coordinates in the image space.
     *
     * @param [in] xyz Spatial coordinates.
     *
     * @returns Voxel coordinates.
     *
     * @sa UpdateTransforms()
     */
    Ivector3d to_ijk(const Fvector3d& xyz) const;

    /**
     * @brief Transform spatial coordinates to voxel coordinates.
     *
     * @param [in] x Spatial coorindate along x axis.
     * @param [in] y Spatial coorindate along y axis.
     * @param [in] z Spatial coorindate along z axis.
     *
     * @returns Voxel coordinates.
     *
     * @sa UpdateTransforms()
     */
    Ivector3d to_ijk(float x, float y, float z) const;

    // -----------------------------------------------------------------------
    // copy image attributes
public:

    /**
     * @brief Copy image size and region.
     *
     * This function copies the image size as specified by the image header
     * and the actual image region from the specified source image. As the image
     * data of this image is, however, not resized by this function, the size
     * of the image region of the @p src image and the one of this image must
     * match. This function is in particular useful to copy the image information
     * required to properly recover the original image domain after a ResizeImage()
     * operation.
     *
     * @param [in] src Source image.
     */
    void CopyRegion(const Image* src);

    /**
     * @brief Copy image transformation, i.e., qform and sform.
     *
     * This function copies all NIfTI header fields related to the transformation
     * of voxel positions to world coordinates from the given image, i.e.,
     *
     * - pixdim
     * - xyzt_units
     * - qform_code
     * - quatern_b
     * - quatern_c
     * - quatern_d
     * - qoffset_x
     * - qoffset_y
     * - qoffset_z
     * - sform_code
     * - srow_x
     * - srow_y
     * - srow_z
     *
     * @note If no qform and/or sform transformation is set, this function will set
     *       the qform and/or sform code to NIFTI_XFORM_SCANNER_ANAT and use the
     *       default transformation as returned by GetQFormTransform() and
     *       GetSFormTransform(), respectively.
     *
     * @param [in] src   Source image.
     * @param [in] qform Copy qform transformation.
     *                   Note that this includes the fields pixdim and xyzt_units.
     * @param [in] sform Copy sform transformation.
     */
    void CopyTransform(const Image* src, bool qform = true,
                                         bool sform = true);

    /**
     * @brief Copy data scaling function.
     *
     * This function copies the NIfTI header fields related to the intensity data
     * from the given image, i.e.,
     *
     * - scl_slope
     * - scl_inter
     * - cal_max
     * - cal_min
     *
     * @param [in] src Source image.
     */
    void CopyDataScaling(const Image* src);

    /**
     * @brief Copy meta data.
     *
     * This function copies the NIfTI-1 meta data from the given image, i.e.,
     *
     * - intent_p1
     * - intent_p2
     * - intent_p3
     * - intent_code
     * - intent_name
     * - descrip
     * - magic
     *
     * Moreover, this function copies the @c filefmt and @c compress attributes
     * of the given Image instance. Hence, it copies all the information related
     * to the format of the file an image was read from and image was read from
     * and will be stored in by WriteNiftiImage().
     *
     * @param [in] src Source image.
     */
    void CopyMetaData(const Image* src);

    // -----------------------------------------------------------------------
    // initialization
protected:

    /**
     * @brief Initialize image during construction.
     *
     * This method is called by constructors to initialize the image instance.
     *
     * Initializes a new image with image origin at (0, 0, 0), and orientation
     * LPS, i.e., transverse flipped.
     *
     * @param [in] img      Pre-allocated image buffer. If NULL, this method
     *                      allocates a new image buffer of the given size and
     *                      and datatype.
     * @param [in] x_dim    Size along x axis.
     * @param [in] y_dim    Size along y axis.
     * @param [in] z_dim    Size along z axis.
     * @param [in] x_pixdim Voxel size along x axis.
     * @param [in] y_pixdim Voxel size along y axis.
     * @param [in] z_pixdim Voxel size along z axis.
     * @param [in] datatype NIfTI-1 datatype code.
     * @param [in] n        Number of components per voxel.
     * @param [in] fmt      Image format.
     * @param [in] owns_img Whether the newly constructed image shall take
     *                      over ownership of the image data. If true, this
     *                      image instance will delete the image data upon
     *                      construction.
     *
     * @throws std::invalid_argument if any of the input arguments has an
     *                               invalid value such as a dimension of zero,
     *                               for example.
     */
    void Init(Data   img,
              int    x_dim,
              int    y_dim,
              int    z_dim,
              float  x_pixdim,
              float  y_pixdim,
              float  z_pixdim,
              short  datatype,
              int    n,
              Format fmt,
              bool   owns_img);

    /**
     * @brief Allocate image data.
     */
    Data AllocateData(int x_dim, int y_dim, int z_dim, int n, short datatype, Format fmt);

    /**
     * @brief Release image data.
     */
    void ReleaseData();

    // -----------------------------------------------------------------------
    // members
public:

    Header     hdr;      ///< Header of image in NIfTI-1 format.
    Data       img;      ///< Image data.
    FileFormat filefmt;  ///< File format of image on disk.
    Format     imgfmt;   ///< Format of image representation.
    Region     region;   ///< Region of image data within image extent specified by @c hdr.
    bool       ownsimg;  ///< Whether this image owns the image data.
    bool       compress; ///< Whether to write image data to compressed image file.

    // The following members have to be updated whenever the image header
    // is changed by calling the method UpdateTransform() and may not be
    // modified directly.
    Transform qto_xyz;    ///< Transforms voxel coordinates to spatial
                          ///< coordinates using the qform transform.
    Transform qto_ijk;    ///< Transforms spatial coordinates to voxel
                          ///< coordinates using the qform transform.
}; // class Image

// ===========================================================================
// image sequence
// ===========================================================================

/**
 * @class Sequence
 * @brief Container for an (ordered) collection of images.
 */
class Sequence : public std::vector<Image*>
{
public:

  /**
   * @brief Default constructor
   */
  Sequence() {}

  /**
   * @brief Copy constructor.
   */
  Sequence(const Sequence &images) {
    for (size_t i = 0; i < images.size(); i++) push_back(images[i]);
  }

  /**
   * @brief Destructor
   */
  virtual ~Sequence() {}

  /**
   * @brief Get i-th volume of image sequence.
   */
  Image *GetVolume(int i)
  {
    return at(i);
  }

  /**
   * @brief Get i-th volume of image sequence.
   */
  const Image *GetVolume(int i) const
  {
    return at(i);
  }

  /**
   * @brief Delete images in sequence and clear list.
   */
  void Clear()
  {
    for (iterator it = begin(); it != end(); ++it) delete *it;
    clear();
  }

}; // class Sequence

/**
 * @class ConstSequence
 * @brief Container for an (ordered) collection of constant images.
 */
class ConstSequence : public std::vector<const Image*>
{
public:

  /**
   * @brief Default constructor
   */
  ConstSequence() {}

  /**
   * @brief Copy constructor.
   */
  ConstSequence(const ConstSequence &images) {
    for (size_t i = 0; i < images.size(); i++) push_back(images[i]);
  }

  /**
   * @brief Copy constructor.
   */
  ConstSequence(const Sequence &images) {
    for (size_t i = 0; i < images.size(); i++) push_back(images[i]);
  }

  /**
   * @brief Destructor
   */
  virtual ~ConstSequence() {}

  /**
   * @brief Get i-th volume of image sequence.
   */
  const Image *GetVolume(int i) const
  {
    return at(i);
  }

  /**
   * @brief Delete images in sequence and clear list.
   */
  void Clear()
  {
    for (iterator it = begin(); it != end(); ++it) delete *it;
    clear();
  }

}; // class ConstSequence

// ===========================================================================
// basic image utility functions
// ===========================================================================

/**
 * @brief Print image region to output stream.
 *
 * @param [in] os     Output stream.
 * @param [in] region Image region.
 *
 * @returns Output stream.
 */
std::ostream& operator<<(std::ostream& os, const Image::Region& region);

/**
 * @brief Convert string to file format enumeration value.
 *
 * @param [in] str Format string, e.g., NIFTI_1, ANALYZE, META,...
 *
 * @returns Corresponding enumeration value or Image::FILE_FORMAT_UNKNOWN.
 */
Image::FileFormat StringToFileFormat(const char* str);

/**
 * @brief Convert Image::NiftiType enumeration value to corresponding Image::FileFormat value.
 *
 * @param [in] type NIfTI-1 type enumeration value.
 *
 * @returns Corresponding Image::FileFormat enumeration value or Image::FILE_FORMAT_UNKNOWN.
 */
Image::FileFormat NiftiTypeToFileFormat(Image::NiftiType type);

/**
 * @brief Convert Image::FileFormat enumeration value to corresponding Image::NiftiType value.
 *
 * @param [in] fmt File format enumeration value.
 *
 * @returns Corresponding Image::NiftiType enumeration value or Image::NIFTI_TYPE_UNKNOWN.
 */
Image::NiftiType FileFormatToNiftiType(Image::FileFormat fmt);

/**
 * @brief Convert file format enumeration value to string.
 *
 * @param [in] fmt File format enumeration value.
 *
 * @returns File format string.
 */
const char* FileFormatToString(Image::FileFormat fmt);

/**
 * @brief Convert string to image format enumeration value.
 *
 * @param [in] str Format string, i.e., DRAMMS, ITK, or FSL.
 *
 * @returns Corresponding enumeration value or Image::FORMAT_UNKNOWN.
 */
Image::Format StringToImageFormat(const char* str);

/**
 * @brief Convert image format enumeration value to string.
 *
 * @param [in] fmt Image format enumeration value.
 *
 * @returns Image format string.
 */
const char* ImageFormatToString(Image::Format fmt);

/**
 * @brief Permute x and y axes in image transformation.
 *
 * @param [in] T Image transformation matrix.
 *
 * @returns Image transformation matrix with x and y axes swapped.
 */
Image::Transform PermuteXY(const Image::Transform T);

/**
 * @brief Permute x and y axes in image header.
 *
 * @param [in] nhdr Image header.
 *
 * @returns Image header with x and y axes permuted.
 *
 * @sa Image::FORMAT_DRAMMS
 */
Image::Header PermuteXY(const Image::Header& nhdr);

/**
 * @brief Get 4x4 qform transformation matrix given image header.
 *
 * @param [in] hdr Image header.
 * @param [in] fmt If DRAMMS format, the default, the x and y axes in the
 *                 header fields are considered to be exchanged.
 *
 * @returns The 4x4 transformation matrix.
 */
Image::Transform GetQFormTransform(Image::Header hdr, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Get 4x4 qform transformation matrix of image.
 *
 * @param [in] image Image.
 *
 * @returns The 4x4 transformation matrix.
 */
inline Image::Transform GetQFormTransform(const Image* image)
{
    return GetQFormTransform(image->hdr, image->imgfmt);
}

/**
 * @brief Get 4x4 sform transformation matrix given image header.
 *
 * @param [in] hdr Image header.
 * @param [in] fmt If DRAMMS format, the default, the x and y axes in the
 *                 header fields are considered to be exchanged.
 *
 * @returns The 4x4 transformation matrix.
 */
Image::Transform GetSFormTransform(Image::Header hdr, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Get 4x4 sform transformation matrix of image.
 *
 * @param [in] image Image.
 *
 * @returns The 4x4 transformation matrix.
 */
inline Image::Transform GetSFormTransform(const Image* image)
{
    return GetSFormTransform(image->hdr, image->imgfmt);
}

/**
 * @brief Get left-right order of image data.
 *
 * @param [in] hdr Image header.
 * @param [in] fmt If DRAMMS format, the default, the x and y axes in the
 *                 header fields are considered to be exchanged.
 *
 * @returns Left-right order of image, i.e., either LR_NEUROLOGICAL or LR_RADIOLOGICAL.
 */
Image::LROrder GetLeftRightOrder(Image::Header hdr, Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Get string naming left-right order of image data.
 *
 * @param [in] order Enumeration value.
 *
 * @returns String corresponding to enumeration value of left-right order of image data.
 */
std::string LROrderToString(Image::LROrder order);

/**
 * @brief Get the image orientation codes.
 *
 * @param [in]  hdr      Image header.
 * @param [out] orient_x Orientation of the first axis.
 * @param [out] orient_y Orientation of the second axis.
 * @param [out] orient_z Orientation of the third axis.
 * @param [in]  fmt      If DRAMMS format, the default, the x and y axes in the
 *                       header fields are considered to be exchanged.
 */
void GetImageOrientation(Image::Header hdr,
                         int& orient_x, int& orient_y, int& orient_z,
                         Image::Format fmt = Image::FORMAT_DRAMMS);

/**
 * @brief Get orientation code from orientation letter.
 *
 * @param [in] c Orientation letter, i.e., L(eft), R(right), P(osterior),
 *               A(anterior), S(uperior), I(inferior).
 *
 * @returns Orientation code, i.e., NIFTI_R2L, NIFTI_L2R, NIFTI_A2P,
 *          NIFTI_P2A, NIFTI_I2S, or NIFTI_S2I.
 */
int CharToOrientation(char c);

/**
 * @brief Get orientation letter from orientation code.
 *
 * @param [in] o Orientation code, i.e., NIFTI_R2L, NIFTI_L2R, NIFTI_A2P,
 *               NIFTI_P2A, NIFTI_I2S, or NIFTI_S2I.
 *
 * @returns Orientation letter, i.e., L(eft), R(right), P(osterior),
 *          A(anterior), S(uperior), I(inferior).
 */
char OrientationToChar(int o);

/**
 * @brief Get orientation from three letter code.
 *
 * @param [in]  str      Three letter orientation code.
 * @param [out] x_orient Orientation code of first dimension.
 * @param [out] y_orient Orientation code of second dimension.
 * @param [out] z_orient Orientation code of third dimension.
 *
 * @returns Whether the three letter code was valid or not.
 */
bool StringToOrientation(const char* str, int& x_orient, int& y_orient, int& z_orient);

/**
 * @brief Get orientation from three letter code.
 *
 * @param [in] x_orient Orientation code of first dimension.
 * @param [in] y_orient Orientation code of second dimension.
 * @param [in] z_orient Orientation code of third dimension.
 *
 * @returns Three letter orientation code.
 */
std::string OrientationToString(int x_orient, int y_orient, int z_orient);


} // namespace dramms


// definition of inline functions
#include "image.hxx"


#endif // _DRAMMS_COMMON_IMAGE_H
