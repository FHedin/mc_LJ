/**
 * \file plugins_lua.h
 *
 * \brief Header file for plugins_lua.c
 *
 * \authors Florent Hedin (University of Basel, Switzerland) \n
 *          Markus Meuwly (University of Basel, Switzerland)
 *
 * \copyright Copyright (c) 2011-2015, Florent Hédin, Markus Meuwly, and the University of Basel. \n
 *            All rights reserved. \n
 *            The 3-clause BSD license is applied to this software. \n
 *            See LICENSE.txt
 *
 */

#ifndef PLUGINS_LUA_H_INCLUDED
#define PLUGINS_LUA_H_INCLUDED

#ifdef LUA_PLUGINS

typedef enum
{
    PAIR=0,
    FFI=1
}LUA_PLUGIN_TYPE;

typedef enum
{
    POTENTIAL=0,
    GRADIENT=1
}LUA_FUNCTION_TYPE;

extern LUA_PLUGIN_TYPE lua_plugin_type;

/*
 * Maximum number of Lua functions it is possible to register
 * Can be redefined when compiling 
 */
#ifndef LUA_MAX_FUNCTIONS_NUMBER
#define LUA_MAX_FUNCTIONS_NUMBER    12
#endif

/*
 * Maximum length for the name of Lua functions registred
 * Can be redefined when compiling 
 */
#ifndef LUA_FUNCTIONS_NAMELEN
#define LUA_FUNCTIONS_NAMELEN   256
#endif

void init_lua(char *plugin_file_name);
void register_lua_function(char *plugin_function_name, LUA_FUNCTION_TYPE type);
void end_lua();

double get_lua_V(ATOM at[], DATA *dat, int32_t candidate);
void get_lua_DV(ATOM at[], DATA *dat, double fx[], double fy[], double fz[]);

double get_lua_V_ffi(ATOM at[], DATA *dat, int32_t candidate);
void get_lua_DV_ffi(ATOM at[], DATA *dat, double fx[], double fy[], double fz[]);

#endif //LUA_PLUGINS

#endif // PLUGINS_LUA_H_INCLUDED
