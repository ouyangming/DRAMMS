##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of software testing.
#
# Unit tests test only single functions, classes, and modules of the software.
# System tests, on the other side, execute the programs of this package with
# given test input data and compare the output of the programs to the expected
# results. For the definition of system tests which perform image regression
# tests, i.e., compare the images written by the test to given baseline images,
# the basis.testdriver executable target imported from BASIS can be used to
# execute the executable to test and perform the regression tests.
#
# Example:
# @code
# basis_add_test (
#  COMMAND basis.testdriver -v
#    --intensity-tolerance 2
#    --compare "t1_bc.nii.gz"
#              "${PROJECT_TESTING_DIR}/baseline/t1_bc.nii.gz"
#    --
#    $<TARGET_FILE:biascorrect>
#      "${PROJECT_TESTING_DIR}/input/t1.nii.gz"
#      t1_bc.nii.gz
# )
# @endcode
#
# If the tests are implemented in C/C++ as supported by CMake's
# create_test_sourcelist() function, a custom test driver can be generated
# using the basis_add_test_driver() function and this test driver then be used
# to execute the tests and perform the regression tests.
#
# Example:
# @code
# basis_add_test_driver (testdriver test1.cxx test2.cxx test3.cxx)
#
# basis_add_test (
#   COMMAND testdriver -v
#     --intensity-tolerance 2
#     --compare "t1_bc.nii.gz"
#               "${PROJECT_TESTING_DIR}/baseline/t1_bc.nii.gz"
#     --
#     biascorrect # test1.cxx defines int biascorrect(int argc, char* argv[])
#      "${PROJECT_TESTING_DIR}/input/t1.nii.gz"
#      t1_bc.nii.gz
# )
#
# # ...
# @endcode
#
# @sa http://www.cmake.org/cmake/help/cmake-2-8-docs.html#command:create_test_sourcelist
#
# The built test executables can be found in @c TESTING_RUNTIME_DIR.
# The test output shall be written to @c TESTING_OUTPUT_DIR. In most cases,
# a subdirectory for each test should be created under this output directory.
# Therefore, the basis_add_test() function creates such subdirectory and sets
# the working directory of the test to this output directory.
#
# Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.<br />
# See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

# ============================================================================
# internal tests
# ============================================================================

if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/internal")
  add_subdirectory ("internal")
endif ()

# ============================================================================
# public tests
# ============================================================================

set (INPUT_DIR    "${CMAKE_CURRENT_SOURCE_DIR}/input")
set (BASELINE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/baseline")

if (APPLE)
  # these very liberal tolerances seem to be required due to the
  # affine registration that is done prior to the deformable one;
  # only after increasing them, the test passed on my
  # MacBook Pro Retina (early 2013) with OS X 10.8.4 and FSL 5.0.1
  # -schuha
  set (WITH_AFFINE_INTENSITY_TOLERANCE          5)
  set (WITH_AFFINE_MAX_NUMBER_OF_DIFFERENCES 1000)
else ()
  set (WITH_AFFINE_INTENSITY_TOLERANCE          5)
  set (WITH_AFFINE_MAX_NUMBER_OF_DIFFERENCES  500)
endif ()


basis_add_test(
  dramms_default_noaffine
  COMMAND
    basis.testdriver
      --compare "H08c_to_H09c_default_noaffine.nii.gz"
                "${BASELINE_DIR}/H08c_to_H09c_default_noaffine.nii.gz"
      --intensity-tolerance 0.01
      --
      $<TARGET_FILE:dramms> -v -v
        -a 0
        -S ${INPUT_DIR}/H08c.nii.gz
        -T ${INPUT_DIR}/H09c.nii.gz
        -O H08c_to_H09c_default_noaffine.nii.gz
        -D def_H08c_to_H09c_noaffine.nii.gz
)


basis_add_test(
  dramms_default_affine
  COMMAND
    basis.testdriver
      --intensity-tolerance       ${WITH_AFFINE_INTENSITY_TOLERANCE}
      --max-number-of-differences ${WITH_AFFINE_MAX_NUMBER_OF_DIFFERENCES}
      --compare "H08c_to_H09c_default_affine.nii.gz"
                "${BASELINE_DIR}/H08c_to_H09c_default_affine.nii.gz"
      --
      $<TARGET_FILE:dramms> -v -v
        -S "${INPUT_DIR}/H08c.nii.gz"
        -T "${INPUT_DIR}/H09c.nii.gz"
        -O H08c_to_H09c_default_affine.nii.gz
)


basis_add_test(
  dramms_warping
  COMMAND
    basis.testdriver
      --compare "warped_H08c_to_H09c.nii.gz"
                "${BASELINE_DIR}/warped_H08c_to_H09c.nii.gz"
      --
      $<TARGET_FILE:dramms> -v -v
        -a 0
        -g 0.2 
        -x 4
        -y 5
        -S "${INPUT_DIR}/H08c.nii.gz"
        -T "${INPUT_DIR}/H09c.nii.gz"
        -L "${INPUT_DIR}/H08c.nii.gz"
        -W warped_H08c_to_H09c.nii.gz
        -r 0
)


basis_add_test(
  dramms_jacobian
  COMMAND
    basis.testdriver
      --intensity-tolerance 0.01
      --compare jac_H08c_to_H09c.nii.gz
                "${BASELINE_DIR}/jac_H08c_to_H09c.nii.gz"
      --
      $<TARGET_FILE:dramms> -v -v
        -a 0
        -g 0.5
        -S "${INPUT_DIR}/H08c.nii.gz"
        -T "${INPUT_DIR}/H09c.nii.gz"
        -J jac_H08c_to_H09c.nii.gz
)


basis_add_test(
  dramms_ravens
  COMMAND
    basis.testdriver
      --compare "ravens_H08c_to_H09c.nii.gz"
                "${BASELINE_DIR}/ravens_H08c_to_H09c_1.nii.gz"
      --
      $<TARGET_FILE:dramms> -v -v
        -a 0
        -g 0.3
        -n 4
        -k 5
        -t 30
        -S "${INPUT_DIR}/H08c.nii.gz"
        -T "${INPUT_DIR}/H09c.nii.gz"
        -L "${INPUT_DIR}/H08c_tissueforravens.nii.gz"
        -l 1
        -R ravens_H08c_to_H09c.nii.gz
)


basis_add_test(
  dramms_initdef
  COMMAND
    basis.testdriver
      --intensity-tolerance       3
      --max-number-of-differences 20
      --compare "H08c_to_H09c_with_initdef.nii.gz"
                "${BASELINE_DIR}/H08c_to_H09c_with_initdef.nii.gz"
      --
      $<TARGET_FILE:dramms> -v -v
        -a 0
        -S "${INPUT_DIR}/H08c.nii.gz"
        -T "${INPUT_DIR}/H09c.nii.gz"
        -d "${INPUT_DIR}/def_H08c_to_H09c_noaffine.nii.gz"
        -O H08c_to_H09c_with_initdef.nii.gz
)


basis_add_test(
  dramms_brainregwarpjac
  COMMAND
    basis.testdriver
      --intensity-tolerance       ${WITH_AFFINE_INTENSITY_TOLERANCE}
      --max-number-of-differences ${WITH_AFFINE_MAX_NUMBER_OF_DIFFERENCES}
      --compare "na01_to_002.nii.gz"
                "${BASELINE_DIR}/na01_RPI_to_002_RPI.nii.gz"
      --compare "na01_to_002_nn.nii.gz"
                "${BASELINE_DIR}/na01_RPI_to_002_RPI_nn.nii.gz"
      --intensity-tolerance       0.1
      --compare "jac_na01_to_002.nii.gz"
                "${BASELINE_DIR}/jac_na01_RPI_to_002_RPI.nii.gz"
      --
      $<TARGET_FILE:dramms> -v -v
        -g 0.4 
        -k 4
        -n 6
        -z 2
        -t 15
        -S "${INPUT_DIR}/na01_ds_analyze.hdr"
        -T "${INPUT_DIR}/002_ds_analyze.hdr"
        -O na01_to_002.nii.gz
        -J jac_na01_to_002.nii.gz
        -W na01_to_002_nn.nii.gz
        -L "${INPUT_DIR}/na01_ds_analyze.hdr"
        -r 0
)


basis_add_test(
  dramms_3dimagewith2slices
  COMMAND
    basis.testdriver
      --intensity-tolerance         3
      --max-number-of-differences  20
      --compare "cardiac_T05_to_T01.img"
                "${BASELINE_DIR}/cardiac_T05_to_T01.img"
      --
      $<TARGET_FILE:dramms> -v -v
        -S "${INPUT_DIR}/cardiacT05c.hdr"
        -T "${INPUT_DIR}/cardiacT01c.hdr"
        -O cardiac_T05_to_T01.hdr
)

    
basis_set_tests_properties(
  dramms_default_noaffine
  dramms_default_affine
  dramms_warping
  dramms_jacobian
  dramms_ravens
  dramms_initdef
  dramms_brainregwarpjac
  dramms_3dimagewith2slices
  PROPERTIES
    LABELS dramms
)
