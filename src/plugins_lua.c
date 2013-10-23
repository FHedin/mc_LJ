/*
 * Copyright (c) 2013, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
 * 
 */
 
#include <stdlib.h>

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>
 
#include "global.h"
#include "plugins_lua.h"
#include "logger.h"

static lua_State* L;

void init_lua(char plugin_file_name[])
{
    int err=0;
    
    /* initialize Lua */
    L = luaL_newstate();

    /* load Lua base libraries */
    luaL_openlibs(L);
    
    /* load the script */
    err = luaL_dofile(L, plugin_file_name);
    if (err != 0)
        LOG_PRINT(LOG_ERROR,"Erro while initialising Lua : %s \n",lua_tostring(L, -1));
    else
        LOG_PRINT(LOG_INFO,"Lua successfully initialised.\n");
   
}

// void register_lua_function(char plugin_function_name[])
// {
// }

void end_lua()
{
    /* cleanup Lua */
    lua_close(L);
}

double get_lua_V(ATOM at[], DATA *dat, int32_t candidate)
{
    uint32_t i,j;
    double dx1,dy1,dz1;
    double dx2,dy2,dz2;
    
    double energy=0.0;
    
    if (candidate==-1)
    {
        for (i=0; i<(dat->natom); i++)
        {
            dx1=at[i].x;
            dy1=at[i].y;
            dz1=at[i].z;
            
            for (j=i+1; j<(dat->natom); j++)
            {
                dx2=at[j].x;
                dy2=at[j].y;
                dz2=at[j].z;
                
//                 LOG_PRINT(LOG_DEBUG,"From get_lua_V full : [i,j] = %d,%d\n",i,j);
                
                /*
                 * lua interface code goes here 
                 */
                
                /* the lua function name */
                lua_getglobal(L, "lj_v_n_m_pair");
    
                // first push arguments to the stack
                lua_pushnumber(L, dx1);
                lua_pushnumber(L, dy1);
                lua_pushnumber(L, dz1);
                
                lua_pushnumber(L, dx2);
                lua_pushnumber(L, dy2);
                lua_pushnumber(L, dz2);
                
                lua_pushnumber(L, at[i].ljp.eps);
                lua_pushnumber(L, at[j].ljp.eps);
                
                lua_pushnumber(L, at[i].ljp.sig);
                lua_pushnumber(L, at[j].ljp.sig);     
                
                lua_pushnumber(L, 4.0);
                lua_pushnumber(L, 12.0);
                lua_pushnumber(L, 6.0);
                
                /* call the function with 13 arguments, return 1 result */
                lua_call(L, 13, 1);
                
                /* get the result and pop arguments for next iteration */
                energy += (double)lua_tonumber(L, -1);
                lua_pop(L, 1);
            }
        }
    }
    else
    {
        i = (uint32_t) candidate;
        
        dx1=at[i].x;
        dy1=at[i].y;
        dz1=at[i].z;
        
        for (j=0; j<(dat->natom); j++)
        {
            if (j!=i)
            {
                dx2=at[j].x;
                dy2=at[j].y;
                dz2=at[j].z;
                
//                 LOG_PRINT(LOG_DEBUG,"From get_lua_V candidate : [i,j] = %d,%d\n",i,j);

                /*
                 * lua interface code goes here 
                 */
                
                /* the lua function name */
                lua_getglobal(L, "lj_v_n_m_pair");
                
                lua_pushnumber(L, dx1);
                lua_pushnumber(L, dy1);
                lua_pushnumber(L, dz1);
                
                lua_pushnumber(L, dx2);
                lua_pushnumber(L, dy2);
                lua_pushnumber(L, dz2);
                
                lua_pushnumber(L, at[i].ljp.eps);
                lua_pushnumber(L, at[j].ljp.eps);
                
                lua_pushnumber(L, at[i].ljp.sig);
                lua_pushnumber(L, at[j].ljp.sig);     
                
                lua_pushnumber(L, 4.0);
                lua_pushnumber(L, 12.0);
                lua_pushnumber(L, 6.0);
                
                /* call the function with 13 arguments, return 1 result */
                lua_call(L, 13, 1);
                
                /* get the result */
                energy += (double)lua_tonumber(L, -1);
                lua_pop(L, 1);
            }
        }
    } // end if-else on candidate
    
//     LOG_PRINT(LOG_DEBUG,"LUA n,m potential : %lf\n",energy);
    
    return energy;
}

