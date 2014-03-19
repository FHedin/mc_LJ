/**
 * \file logger.h
 *
 * \brief Header file for managing logging, see \b #logger.c
 *
 * \authors Florent Hedin (University of Basel, Switzerland) \n
 *          Markus Meuwly (University of Basel, Switzerland)
 *
 * \copyright Copyright (c) 2014, Florent Hedin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#ifndef LOGGER_H_INCLUDED
#define LOGGER_H_INCLUDED

#include <stdint.h>

/**
*   \enum       LOG_LEVELS
*
*   \brief      The available logging Levels.
*
*   \details This will control what is written
*   to logging files ; those files are named :
*   \li error.log
*   \li warning.log
*   \li info.log
*   \li debug.log
*
*   The higher the level, the more text is written to those text files, so be careful with long simulations
*   coupled to high levels of logging.
*
*
*   The different levels, members of the \b enum \b #LOG_LEVELS are :
*   \li \b #LOG_NOTHING : At this level there is no logging at all.
*   \li \b #LOG_ERROR : An Error is a critical event which will cause the program to abort immediately,
*       for example a missing parameter for which no default value is provided.
*   \li \b #LOG_WARNING : A Warning is a non-critical event, i.e. the event is just reported but the program
*       execution continues, for example a missing or wrong parameter for which a default value is going to be used.
*   \li \b #LOG_INFO : An Info message reports to the user a possibly useful information, an event with no consequence for the simulation.
*       For example, that the trajectory was succesfully written as planed.
*   \li \b #LOG_DEBUG : A debug message is really technical and verbose, for example dumping a variable.
*       This should be enabled for bug tracking.
*
*   The logging is progressive, enabling \b #LOG_INFO means that higher levels \b #LOG_ERROR and \b #LOG_WARNING
*       are also enabled.
*
*   The default is set to \b #LOG_WARNING. \n
*
*/
typedef enum
{
    LOG_NOTHING = 0, /*!< No log file is created : this is not recommended as no information is reported */
    LOG_ERROR = 1,   /*!< Only errors are reported to \b error.log */
    LOG_WARNING = 2, /*!< Warnings are reported to \b warning.log , and also errors to \b error.log. This is the default. */
    LOG_INFO = 3,    /*!< Info messages are reported to \b info.log , and also warnings and errors to their respective files. */
    LOG_DEBUG = 4    /*!< A large amount of debugging messages are written to \b debug.log . The previous levels are still written to their respective files. */
} LOG_LEVELS;

extern LOG_LEVELS LOG_SEVERITY;

// #define LOG_PRINT(...) log_print(__FILE__, __LINE__, __VA_ARGS__ )

void init_logfiles();
void close_logfiles();
char* get_loglevel_string();
char* get_time();

uint32_t LOG_PRINT(LOG_LEVELS mesg_severity, char *fmt, ...);
uint32_t LOG_PRINT_SHORT(LOG_LEVELS mesg_severity, char *fmt, ...);

#endif // LOGGER_H_INCLUDED
