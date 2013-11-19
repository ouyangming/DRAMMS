##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of software.
#
# Copyright (c) 2011, 2012 University of Pennsylvania. All rights reserved.
# See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

# ----------------------------------------------------------------------------
## @brief Add auxiliary executable for internal use by dramms script only.
macro (dramms_add_advanced_tool EXENAME)
  basis_add_executable(${EXENAME}.cxx LIBEXEC ${ARGN})
  basis_target_link_libraries(${EXENAME} common)
endmacro ()

# ----------------------------------------------------------------------------
## @brief Add auxiliary executable also for use by end users.
macro (dramms_add_user_tool EXENAME NAME_SUFFIX)
  basis_add_executable(${EXENAME}.cxx ${ARGN})
  basis_target_link_libraries(${EXENAME} common)
  basis_set_target_properties (${EXENAME} PROPERTIES OUTPUT_NAME "dramms-${NAME_SUFFIX}")
endmacro ()

# ----------------------------------------------------------------------------
# auxiliary executables
dramms_add_advanced_tool(CheckTwoImages)
dramms_add_advanced_tool(GenerateGaborListFile)
dramms_add_advanced_tool(GenerateImageListFile)
dramms_add_advanced_tool(CalculateImageSimilarity)

# ----------------------------------------------------------------------------
# tools to combine files
dramms_add_advanced_tool(CombineFeatureLists)
dramms_add_user_tool(CombineTransforms combine)

# ----------------------------------------------------------------------------
# basic processing tools
dramms_add_advanced_tool(CalculateGaborTextures)
dramms_add_advanced_tool(MatchHistograms)

dramms_add_user_tool(ApplyTransform       warp)
dramms_add_user_tool(CalculateJacobianMap jacobian)
dramms_add_user_tool(CalculateRavensMaps  ravens)
dramms_add_user_tool(ImageOperation       imgop)
dramms_add_user_tool(TransformOperation   defop)

# ----------------------------------------------------------------------------
# conversion tools
dramms_add_user_tool(ConvertDeformation convert)
dramms_add_advanced_tool(ConvertImage)

# ----------------------------------------------------------------------------
# deformable registration
dramms_add_advanced_tool(Deform2D)
dramms_add_advanced_tool(Deform3D DiscreteOptimizationUtilities.cxx GeneralUtilities.cxx)
basis_target_link_libraries(Deform3D "${DRAMMSFastPD_LIBRARY}")
