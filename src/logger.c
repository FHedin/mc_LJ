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

#include <time.h>
#include <string.h>

#include <stdarg.h>

#include "logger.h"

// those variabes arepersistent but only accessible from this file
static FILE *F_ERROR , *F_WARN , *F_INFO , *F_DEBUG ;
static char EVENT_DATE[32];

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
 * \brief Puts in a string the current date, later used for logging.
 */
static void get_time()
{
    time_t rawtime;
    struct tm * timeinfo;

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    
    strftime(EVENT_DATE,32,"%d/%m/%Y-%T-%Z",timeinfo);
}

/**
 * \brief log_print
 * 
 * \param filename
 * \param line
 * \param mesg_severity
 * \param fmt
 * \param ...
 * 
 * 
 */
uint32_t LOG_PRINT(LOG_LEVELS mesg_severity, char *fmt, ...)
{
    // if the severity of the current message is at least equal to the gloval level we print it
    if (mesg_severity <= LOG_SEVERITY)
    {
        
        va_list list;
        char* fmt_buff;
        
        FILE *FP=NULL;

        // get the date string 
	    get_time();
        
        // depending of the severity, chose the correct file for writting
        switch(mesg_severity)
        {
            case LOG_ERROR:
                FP = F_ERROR;
                fprintf(FP,"[Error @ ");
                break;
                
            case LOG_WARNING:
                FP = F_WARN;
                fprintf(FP,"[Warning @ ");
                break;
                
            case LOG_INFO:
                FP = F_INFO;
                fprintf(FP,"[Info @ ");
                break;
                
            case LOG_DEBUG:
                FP = F_DEBUG;
                fprintf(FP,"[Debug @ ");
                break;
                
            default:
                return -1;
                break;
        }
        
        // write date contained in EVENT_DATE
        fprintf(FP,"%s]\t",EVENT_DATE);
        
        // prepare processing of the variable arguments list
        // fmt is the last non-optional argument
        va_start( list, fmt );
        
        /* fmt is a printf-like string containing % for indicating which type of variables
         * are going to be printed : by iterating through this string we can know what are
         * the primitive types of the passed variables.
        */
        // This for means that we iterate until *fmt_buff == '\0', the end of string
        for ( fmt_buff = fmt ; *fmt_buff ; ++fmt_buff )
        {
            // if not % the character is just forwarded to the log file
            if( *fmt_buff != '%' )
            {
                fputc(*fmt_buff,FP);
            }
            // else we need to determine the type of the variable
            else
            {
                switch(*++fmt_buff)
                {
                    // a variable argument is a string
                    case 's':
                    {
                        char *tmp_s = va_arg( list, char* );
                        fprintf(FP,"%s",tmp_s);
                        continue;
                    }
                    break;
                    
                    // a variable argument is an integer
                    case 'd':
                    {
                        int32_t tmp_d = va_arg( list, int32_t );
                        fprintf(FP,"%d",tmp_d);
                        continue;
                    }
                    break;
                    
                    default:
                        fputc(*fmt_buff,FP);
                        break;
                }
            }
        } // end of for on the format string
 
        va_end( list );
//         fputc('\n',FP);
        
        return 0;

    } // end of if (mesg_severity >= LOG_SEVERITY)
    
    return 1;
}
