/** 
 * \file logger.c
 *
 * \brief Logging functions with a basic severity/level check.
 *
 * \authors Florent Hedin (University of Basel, Switzerland) \n
 *          Markus Meuwly (University of Basel, Switzerland)
 * 
 * \copyright Copyright (c) 2013, Florent Hedin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *            
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <time.h>
#include <string.h>

#include "logger.h"

static FILE *F_ERROR , *F_WARN , *F_INFO , *F_DEBUG ;
static char EVENT_DATE[32];

/**
 * \brief   Prepares LOGGING I/O if necessary, depending of the value of \b #LOG_SEVERITY
 */
void init_logfiles()
{
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
 * \brief Puts in a string the current date, then used for logging.
 */
static void print_time()
{
    time_t rawtime;
    struct tm * timeinfo;

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    
    strftime(EVENT_DATE,32,"%d/%m/%Y-%T-%Z",timeinfo);
}

/**
 * \brief log_print
 * \param filename
 * \param line
 * \param fmt
 */
void log_print(char* filename, int line, char *fmt,...)
{
    va_list list;
    char *p, *r;
    int e;

    

    fprintf(fp,"%s ",print_time());
    va_start( list, fmt );

    for ( p = fmt ; *p ; ++p )
    {
        if ( *p != '%' )//If simple string
        {
            fputc( *p,fp );
        }
        else
        {
            switch ( *++p )
            {
                /* string */
            case 's':
            {
                r = va_arg( list, char * );

                fprintf(fp,"%s", r);
                continue;
            }

            /* integer */
            case 'd':
            {
                e = va_arg( list, int );

                fprintf(fp,"%d", e);
                continue;
            }

            default:
                fputc( *p, fp );
            }
        }
    }
    va_end( list );
    fprintf(fp," [%s][line: %d] ",filename,line);
    fputc( '\n', fp );
    SESSION_TRACKER++;
    fclose(fp);
}
