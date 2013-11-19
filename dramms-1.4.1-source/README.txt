
  Section of Biomedical Image Analysis
  Department of Radiology
  University of Pennsylvania
  3600 Market Street, Suite 380
  Philadelphia, PA 19104

  Web:   http://www.rad.upenn.edu/sbia/
  Email: sbia-software at uphs.upenn.edu

  Copyright (c) 2012, 2013 University of Pennsylvania. All rights reserved.
  See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.



INTRODUCTION
============

  Deformable Registration via Attribute Matching and Mutual-Saliency Weighting (DRAMMS)
  is a software package designed for 2D-to-2D and 3D-to-3D image registration tasks.

  Some typical applications of DRAMMS include:
  - Cross-subject registration of the same organ (can be brain, breast, cardiac, etc)
  - Mono- and Multi-modality registration (MRI, CT, histology)
  - Longitudinal registration (pediatric brain growth, cancer development, mouse brain development, etc)
  - Registration under missing correspondences (e.g., vascular lesions, tumors, histological cuts)

  DRAMMS is implemented as a Unix command-line tool. It is fully automatic and easy to use.
  Users input two images and DRAMMS will output the registered image along with the associated deformation.
  No need for pre-segmentation of structures, prior knowledge, or human initialization/interventions.



PACKAGE OVERVIEW
================

  Source Package
  --------------

  - BasisProject.cmake   Meta-data used by BASIS to configure the project.
  - CMakeLists.txt       Root CMake configuration file.
  - build/               Contains CMake configuration file for bundle build.
  - config/              Package configuration files.
  - doc/                 Software documentation such as the software manual.
  - src/                 Source code files.
  - test/                Implementation of software tests and corresponding data.

  - AUTHORS.txt          A list of the people who contributed to this software.
  - COPYING.txt          The copyright and license notices.
  - ChangeLog.txt        History of changes and releases.
  - INSTALL.txt          Build and installation instructions.
  - README.txt           This readme file.


  Binary Package
  --------------

  Please refer to the INSTALL file for details on where the built executables
  and libraries, the auxiliary data, and the documentation files are installed.



LICENSING
=========

  See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.



INSTALLATION
============

  See http://www.rad.upenn.edu/sbia/software/dramms/installation.html or INSTALL file.



DOCUMENTATION
=============

  See the software manual [3] for details on the software including a demonstration
  of how to apply the software tools provided by this package.



UTILITY PROGRAMS
================

  ATTENTION: The following documentation of the auxiliary executables was mainly used
             during the design of the DRAMMS software. Please see the software manual [3]
             and the help output of each built executable for more up-to-date help.

  ConvertImage
  ------------

  This program is used to change the file format, datatype, and/or image orientation
  of any supported scalar image. Moreover, it can be used to rescale the image
  intensities, e.g., to the maximum range of image intensities supported by the datatype.

  - Change image orientation:
      ConvertImage -O LPS image.nii
      ConvertImage -O LPS image.nii image_lps.nii

  - Change file format from .nii to .hdr/.img:
      ConvertImage image.nii image.hdr

  - Change datatype to DT_UNSIGNED_CHAR:
      ConvertImage -t DT_UNSIGNED_CHAR image.nii image_byte.nii

  - Change datatype to DT_UNSIGNED_CHAR and rescale intensities to [0, 255]:
      ConvertImage -t uchar -S image.nii image_byte.nii
  
  - Change datatype to DT_SIGNED_SHORT and reorient image to RAI:
      ConvertImage -t short -O RAI image.nii image_byte.nii

  - Read raw image data and write NIfTI-1 image:
      ConvertImage -d 128,128,64 -p 1,1,2 image.raw image.nii


  ImageOperation
  --------------

  This program performs common operations on a scalar image. It reads in one image,
  performs one or more operations, and writes the resulting image to the specified
  output file. If no output file is specified, it overwrites the input image.

  - Sample image intensities on specified image grid.

    Image (sub-)grid specified by number of voxels and voxel size:

      ImageOperation -d 64,64,32 image.nii out.nii
      ImageOperation -p 3.5,3.5,6 image.nii out.nii
      ImageOperation -d 64,64,32 -p 3.5,3.5,6 image.nii out.nii

      Note: Performs trilinear interpolation.

      Note: If only -d or -p is specified, the extent of the image grid
            corresponds to the same extent as the one of the input image.
            Otherwise, if both -d and -p are specified, the extent of the
            image grid may vary. Both, the input and output image grids are
            in this case overlayed such that the center of the image grids
            coincides.

    Image sub-grid specified by bounding box of foreground mask:

      ImageOperation -m fgmask.nii image.nii subimage.nii

    Image sub-grid specified by first and last voxel indices:

      ImageOperation -x 20,100 image.nii subimage.nii
      ImageOperation -y 50,75  image.nii subimage.nii
      ImageOperation -z 40,60  image.nii subimage.nii
      ImageOperation -x 20,100 -y 50,75 -z 40,60 image.nii


  ConvertDeformation
  ------------------

  1. Read deformation field given in specified format (default DRAMMS) either from
     one single image file or three separate scalar component images. Each of the
     input images may be given in one of the following file formats:

       - NIfTI-1
       - ANALYZE 7.5
       - MetaImage
       - raw image file

     In case of a raw image file, the -d and -p options have to be given as these
     specify the image and voxel size, respectively.

  2. Optionally, change representation of deformation field (e.g., DRAMMS -> ITK).

  3. Write deformation field either to a single image file or to separate
     scalar component images (note that only selected components may be written only).
     If requested, write magnitude image of deformation field.

  Examples:

    Read DRAMMS deformation field and write component image(s):
      ConvertDeformation def.nii -x defx.nii
      ConvertDeformation def.nii -y defy.nii
      ConvertDeformation def.nii -z defz.nii
      ConvertDeformation def.nii -x defx.nii -y defy.nii -z defz.nii

    Read ITK deformation field and write component image(s):
      ConvertDeformation def.nii -f ITK -x defx.nii -y defy.nii -z defz.nii

    Read DRAMMS deformation field and write magnitude image:
      ConvertDeformation def.nii -m defmag.nii

    Read ITK deformation field and save it in DRAMMS format:
      ConvertDeformation def.nii -f ITK,DRAMMS -o outdef.nii
      ConvertDeformation -f ITK,DRAMMS def.nii
      ConvertDeformation defx.nii defy.nii defz.nii -f ITK,DRAMMS -o outdef.nii
      ConvertDeformation defx.nii defy.nii defz.nii -f ITK,DRAMMS
                         -x outx.nii -y outy.nii -z outz.nii

    Read component images and write single deformation field file:
      ConvertDeformation defx.nii defy.nii defz.nii -o def.nii

    Read deformation field from raw data file and write component image(s):
      ConvertDeformation -d 128,128,64 -p 1,1,1 def.raw
                         -x outx.nii -y outy.nii -z outz.nii
      ConvertDeformation -d 128,128,64 -p 1,1,1 -f ITK def.img
                         -x outx.nii -z outz.nii


  TransformationOperation
  -----------------------

  This program reads in a deformation field, processes it, and writes the modified
  deformation field to the specified output file.  If no second positional argument,
  i.e., the name of the output file is specified, it overwrites the input file.

  - Print certain quantities of deformation field such as the maximum displacement
    in x, y, and z.

      TransformationOperation def.nii

  - Print certain quantities (e.g. displacement) of one selected (sub-)voxel.

      TransformationOperation -c 120.3,100,43.78 def.nii

  - Sample displacements on specified image grid.

    Image (sub-)grid specified by number of voxels and voxel size:

      TransformationOperation -d 64,64,32 def.nii outdef.nii
      TransformationOperation -p 3.5,3.5,6 def.nii outdef.nii
      TransformationOperation -d 64,64,32 -p 3.5,3.5,6 def.nii outdef.nii

      Note: Performs trilinear interpolation.

      Note: If only -d or -p is specified, the extent of the image grid
            corresponds to the same extent as the one of the input image.
            Otherwise, if both -d and -p are specified, the extent of the
            image grid may vary. Both, the input and output image grids are
            in this case overlayed such that the center of the image grids
            coincides.

    Image sub-grid specified by first and last voxel indices:

      TransformationOperation -x 20,90 def.nii subdef.nii
      TransformationOperation -y 50,75 def.nii subdef.nii
      TransformationOperation -z 40,60 def.nii subdef.nii
      TransformationOperation -x 20,90 -y 50,75 -z 40,60 def.nii

  - Eliminate displacements within ROI.

    Image ROI specified by foreground mask, i.e., eliminates displacements at
    voxels outside the mask (value == 0):

      TransformationOperation -m fgmask.nii def.nii outdef.nii

    Image ROI specified by background mask, i.e., eliminates displacements at
    voxels inside the mask (value != 0):

      TransformationOperation -b bgmask.nii def.nii outdef.nii

    Image ROI specified by first and last voxel index:

      TransformationOperation -X 20,90 def.nii subdef.nii
      TransformationOperation -Y 50,75 def.nii subdef.nii
      TransformationOperation -Z 40,60 def.nii subdef.nii
      TransformationOperation -X 20,90 -Y 50,75 -Z 40,60 def.nii

  - Invert transformation.

    TransformationOperation -r AtoB.xfm BtoA.xfm
    TransformationOperation -r defAtoB.nii defBtoA.nii

  - Smooth deformation field.

    TransformationOperation -s def.nii
    TransformationOperation -s def.nii outdef.nii

  - Clamp displacements.

    TransformationOperation -t 10 def.nii

  More than one operation can be applied in one program execution.
  The order in which the operations are applied corresponds to the order in which
  the corresponding option flags are given on the command-line.

    TransformationOperation -s -d 128,128,64 -t 3.5 -r def.nii outdef.nii

  This command does

    1. Smooth the input deformation field
    2. Sample the smoothed deformation field on the image grid of size
       128x128x64 voxels using the voxel size of the input field.
    3. Clamp displacements in each dimension which exceed the value of 3.5.
    4. Invert the resulting deformation field.

  and write the final deformation field to the file outdef.nii.


  CombineTransformations
  ----------------------

  This program requires two input files and writes one output file. Both in input
  files represent a transformation, either a deformable transformation represented by
  a displacement vector field, or an affine transformation represented by a 4x4 matrix.

  - Add two deformation fields:
      CombineTransformations -a defA.nii defB.nii defA+B.nii

  - Subtract two deformation fields:
      CombineTransformations -s defA.nii defB.nii defA-B.nii

  - Calculate mean deformation:
      CombineTransformations -m defA.nii defB.nii defmean.nii

  - Concatenate transformations:
      CombineTransformations [-c] defAtoB.nii defBtoC.nii defAtoC.nii
      CombineTransformations [-c] -r B.nii defAtoB.nii matBtoC.xfm defAtoC.nii
      CombineTransformations [-c] -r A.nii matAtoB.xfm defBtoC.nii defAtoC.nii
      CombineTransformations [-c] matAtoB.xfm matBtoC.xfm matAtoC.xfm



REFERENCES
==========

  [1] Yangming Ou, Aristeidis Sotiras, Nikos Paragios, Christos Davatzikos.
      DRAMMS: Deformable registration via attribute matching and mutual-saliency weighting.
      Medical Image Analysis 15(4): 622-639 (2011).

  [2] Yangming Ou, Christos Davatzikos.
      DRAMMS: Deformable Registration via Attribute Matching and Mutual-Saliency Weighting.
      IPMI 2009: 50-62.

  [3] http://www.rad.upenn.edu/sbia/software/dramms/manual.html
