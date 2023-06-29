/*
  This file is part of the SC Library, version 3.
  The SC Library provides support for parallel scientific applications.

  Copyright (C) 2019-2023 individual authors, Scott Aiton

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice,
  this list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.
*/

/* We provide string functions to not make iniparser depend on libsc. */

#ifndef INISTRING_H
#define INISTRING_H

/* we set the GNU feature test macro before including anything */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE             /**< Enable additional functionality. */
#endif
#include <stdio.h>

#ifdef __cplusplus
extern              "C"
{
#if 0
}
#endif
#endif

/** Provide a string copy function.
 * \param [out] dest    Buffer of length at least \a size.
 *                      On output, not touched if NULL or \a size == 0.
 * \param [in] size     Allocation length of \a dest.
 * \param [in] src      Null-terminated string.
 * \return              Equivalent to \ref
 *                      sc3_snprintf (dest, size, "%s", src).
 */
void                ini_strcopy (char *dest, size_t size, const char *src);

/** Wrap the system snprintf function, allowing for truncation.
 * The snprintf function may truncate the string written to the specified length.
 * In some cases, compilers warn when this may occur.
 * For the usage in sc3, this is permitted behavior and we avoid the warning.
 * \param [out] str     Buffer of length at least \a size.
 *                      On output, not touched if NULL or \a size == 0.
 *                      Otherwise, "" on snprintf error or the proper result.
 * \param [in] size     Allocation length of \a str.
 * \param [in] format   Format string as in man (3) snprintf.
 */
void                ini_snprintf (char *str, size_t size,
                                  const char *format, ...)
  __attribute__((format (printf, 3, 4)));

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif /* !INISTRING_H */
