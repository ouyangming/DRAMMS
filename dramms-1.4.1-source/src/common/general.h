/**
 * @file  general.h
 * @brief Common constants and macros.
 *
 * Copyright (c) 2011 University of Pennsylvania. All rights reserved.
 * See https://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
 *
 * Contact: SBIA Group <sbia-software at uphs.upenn.edu>
 */

#pragma once
#ifndef _DRAMMS_GENERAL_H
#define _DRAMMS_GENERAL_H

#include <math.h>   // floor(), ceil()
#include <iostream> // cout, cerr, endl
#include <sstream>  // istringstream
#include <string>   // getline(), string
#include <vector>


namespace dramms {


// ===========================================================================
// constants
// ===========================================================================

#define YYES     1
#define NNO      0
#define TRUE     1
#define FALSE    0
#define OOKK     1


#define E_PI     3.1415926536
#define G_PI     3.1415926
#define INFINITE 1000000000

#define MaxOfUC  255 

// in case we need to label brain tissue or structure. 
// avoid use if possible.
#define BG    0
#define CSF   10
#define VN    50
#define GM    150
#define WM    250
#define Z_THICKINK    1.5 

// ===========================================================================
// debug / verbose output
// ===========================================================================


/**
 * @brief Enumeration of global debug level used by library functions.
 *
 * Depending on the debug level, the library function may or may not print
 * errors, warnings, and additional verbose output to STDERR or STDOUT,
 * respectively.
 */
enum DebugLevel
{
    DEBUG_QUIET, ///< Do not print anything.
    DEBUG_ERROR, ///< Print only errors.
    DEBUG_WARN,  ///< Print errors and warnings.
    DEBUG_ALL    ///< Print all messages.
};

extern DebugLevel g_debug; ///< Global debug level used by library functions.

/// Print warning to @c stderr if debug level is set to @c DEBUG_LEVEL_WARN or above.
#define DRAMMS_WARN(msg) \
    { \
        if (::dramms::g_debug >= ::dramms::DEBUG_WARN) { \
            ::std::cerr << msg << ::std::endl; \
        } \
    }

/// Print error message to @c stderr if debug level is set to @c DEBUG_LEVEL_ERROR or above.
#define DRAMMS_ERROR(msg) \
    { \
        if (::dramms::g_debug >= ::dramms::DEBUG_ERROR) { \
            ::std::cerr << msg << ::std::endl; \
        } \
    }

/// Print verbose message to @c stdout if debug level is set to @c DEBUG_LEVEL_ALL.
#define DRAMMS_MSG(msg) \
    { \
        if (::dramms::g_debug >= ::dramms::DEBUG_ALL) { \
            ::std::cout << msg << ::std::endl; \
        } \
    }

// ===========================================================================
// macros
// ===========================================================================

#define MIN(x,y)  ((x) < (y) ? (x) : (y))
#define MAX(x,y)  ((x) > (y) ? (x) : (y))


// ===========================================================================
// numeric
// ===========================================================================

/**
 * @brief Clamp floating-point value.
 *
 * @param [in] v   Value.
 * @param [in] min Minimum allowed value.
 * @param [in] max Maximum allowed value.
 *
 * @returns Value if within specified range, @p min if @p v is smaller than
 *          @p min, and @p max if @p v is greater than @p max.
 */
inline
float clamp(float v, float min, float max)
{
    return (v < min ? min : (v > max ? max : v));
}

/**
 * @brief Round floating-point value.
 */
inline
float round(float v)
{
    return (v > 0.0f) ? floor(v + 0.5f) : ceil(v - 0.5f);
}

// ===========================================================================
// string manipulation
// ===========================================================================

/**
 * @brief Split string at delimiters.
 *
 * @param [in] text  Input string.
 * @param [in] delim Delimiter.
 *
 * @returns Delimited parts of string, excluding the delimiter.
 */
inline
std::vector<std::string> split(const std::string& text, char delim = '\n')
{
    std::vector<std::string> tokens;
    std::istringstream iss(text);
    std::string token;
    while (std::getline(iss, token, delim)) {
        tokens.push_back(token);
    }
    return tokens;
}


} // namespace dramms


#endif // _DRAMMS_GENERAL_H
