/*
 * Copyright (c) 2013, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
 * 
 */

#ifndef PLUGINS_LUA_H_INCLUDED
#define PLUGINS_LUA_H_INCLUDED

void init_lua(char plugin_file_name[]);
// void register_lua_function(char plugin_function_name[]);
void end_lua();

/* compatible with pointer of type 
 'double (*get_ENER)' from ener.h*/
double get_lua_V(ATOM at[], DATA *dat, int32_t candidate);

#endif // PLUGINS_LUA_H_INCLUDED
