/**
 * \file logger.c
 *
 * \brief Logging functions with a basic severity/level check.
 *
 * \authors Florent Hedin (University of Basel, Switzerland) \n
 *          Markus Meuwly (University of Basel, Switzerland)
 *
 * \copyright Copyright (c) 2011-2015, Florent Hedin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include <time.h>
#include <string.h>

#include <stdarg.h>

#include "logger.h"

// those variabes are persisting but only accessible from this file
static FILE *F_ERROR , *F_WARN , *F_INFO , *F_DEBUG ;
static time_t rawtime;

/**
 * \brief   Prepares LOGGING I/O if necessary, depending of the value of \b #LOG_SEVERITY
 */
void init_logfiles()
{
    F_ERROR = F_WARN = F_INFO = F_DEBUG = NULL;

    if(LOG_SEVERITY > LOG_NOTHING )
        F_ERROR = fopen("error.log","wt");

    if(LOG_SEVERITY > LOG_ERROR )
        F_WARN = fopen("warning.log","wt");

    if(LOG_SEVERITY > LOG_WARNING )
        F_INFO = fopen("info.log","wt");

    if(LOG_SEVERITY > LOG_INFO )
        F_DEBUG = fopen("debug.log","wt");
}

/**
 * \brief   Before exiting the program, closes properly the logging files.
 */
void close_logfiles()
{
    if(LOG_SEVERITY > LOG_NOTHING )
        fclose(F_ERROR);

    if(LOG_SEVERITY > LOG_ERROR )
        fclose(F_WARN);

    if(LOG_SEVERITY > LOG_WARNING )
        fclose(F_INFO);

    if(LOG_SEVERITY > LOG_INFO )
        fclose(F_DEBUG);
}

/**
 * \brief   Returns a string corresponding to the logging level.
 * \return An array of char containing the description of the logging level, for example "LOG_ERROR"
 */
char* get_loglevel_string()
{
    char* str_log=NULL;

    switch(LOG_SEVERITY)
    {
    case LOG_NOTHING:
        str_log = "LOG_NOTHING";
        break;

    case LOG_ERROR:
        str_log = "LOG_ERROR";
        break;

    case LOG_WARNING:
        str_log = "LOG_WARNING";
        break;

    case LOG_INFO:
        str_log = "LOG_INFO";
        break;

    case LOG_DEBUG:
        str_log = "LOG_DEBUG";
        break;

    default:
        str_log = "UNKNOWN";
        break;
    }

    return str_log;
}

/**
 * \brief Puts in a string the current date, later used for logging.
 *
 * \details This function returns the date with the following format :
 * \li DAY/MONTH/YEAR-HH:MM:SS
 * \li For example : 17/10/2013-16:26:36
 *
 * \param date A character string filled with the date and hour
 */
static void get_time_ptr(char date[32])
{
    struct tm * timeinfo;
    time_t newtime;

    time(&newtime);

    if(newtime != rawtime)
    {
        rawtime = newtime;
        timeinfo = localtime(&rawtime);
        strftime(date,32,"%d/%b/%Y-%H:%M:%S",timeinfo);
    }
}

/**
 * \brief Puts in a string the current date, and returns it
 *
 * \details This function returns the date with the following format :
 * \li DAY/MONTH/YEAR-HH:MM:SS
 * \li For example : 17/10/2013-16:26:36
 *
 * \return A character string filled with the date and hour
 */
char* get_time()
{
    static char date[32];
    struct tm * timeinfo;
    time_t newtime;

    time(&newtime);

    timeinfo = localtime(&newtime);
    strftime(date,32,"%d/%b/%Y-%H:%M:%S",timeinfo);

    return date;
}

/**
 * \brief log_print
 *
 * \param mesg_severity
 * \param fmt
 * \param ...
 */
uint32_t LOG_PRINT(LOG_LEVELS mesg_severity, char *fmt, ...)
{
    // if the severity of the current message is at least equal to the gloval level we print it
    if (mesg_severity <= LOG_SEVERITY)
    {

        va_list list;

        FILE *FP=NULL;

        static char event_date[32]="";
        char message[50]="";

        // get the date string
        get_time_ptr(event_date);

        // depending of the severity, chose the correct file for writting
        switch(mesg_severity)
        {
        case LOG_ERROR:
            FP = F_ERROR;
            sprintf(message,"[Error @ ");
            break;

        case LOG_WARNING:
            FP = F_WARN;
            sprintf(message,"[Warning @ ");
            break;

        case LOG_INFO:
            FP = F_INFO;
            sprintf(message,"[Info @ ");
            break;

        case LOG_DEBUG:
            FP = F_DEBUG;
            sprintf(message,"[Debug @ ");
            break;

        default:
            return -1;
            break;
        }

        // write date contained in event_date to message
        strncat(message,event_date,32);
        strncat(message,"]\t",2);

        // print message to file
        fprintf(FP,"%s",message);

        // prepare processing of the variable arguments list
        // fmt is the last non-optional argument
        va_start(list,fmt);

        // now forward what have to be printed to vfprintf
        vfprintf(FP,fmt,list);

        va_end(list);

        return 0;

    } // end of if (mesg_severity >= LOG_SEVERITY)

    return 1;
}

/**
 * \brief Same as \b #LOG_PRINT but without the date appended at the beginning of the line
 *
 * \param mesg_severity
 * \param fmt
 * \param ...
 */
uint32_t LOG_PRINT_SHORT(LOG_LEVELS mesg_severity, char *fmt, ...)
{
    // if the severity of the current message is at least equal to the gloval level we print it
    if (mesg_severity <= LOG_SEVERITY)
    {

        va_list list;

        FILE *FP=NULL;

        // depending of the severity, chose the correct file for writting
        switch(mesg_severity)
        {
        case LOG_ERROR:
            FP = F_ERROR;
            break;

        case LOG_WARNING:
            FP = F_WARN;
            break;

        case LOG_INFO:
            FP = F_INFO;
            break;

        case LOG_DEBUG:
            FP = F_DEBUG;
            break;

        default:
            return -1;
            break;
        }

        // prepare processing of the variable arguments list
        // fmt is the last non-optional argument
        va_start(list,fmt);

        // now forward what have to be printed to vfprintf
        vfprintf(FP,fmt,list);

        va_end(list);

        return 0;

    } // end of if (mesg_severity >= LOG_SEVERITY)

    return 1;
}
