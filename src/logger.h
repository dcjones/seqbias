/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef PEAKOLATOR_LOGGER_H
#define PEAKOLATOR_LOGGER_H

/**
 * \file
 * \brief C interface to the log writer class.
 *
 * I have a few pure C source files still, and I want everything to use the same
 * logging sysetm, thus I provide a minimal C interface.
 */


#ifdef __cplusplus
extern "C" {
#endif


void logger_abort(const char* fmt, ...);


#ifdef __cplusplus
}
#endif


#endif


