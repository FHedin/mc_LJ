/**
 * \file plugins_lua.c
 *
 * \brief File containing functions used for interacting with the user defined lua scripts
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

#ifdef LUA_PLUGINS

#include <stdlib.h>
#include <string.h>

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#include "global.h"
#include "plugins_lua.h"
#include "logger.h"

static lua_State* L;

static uint32_t is_lua_used;

static char lua_file[FILENAME_MAX];

static char lua_function[LUA_MAX_FUNCTIONS_NUMBER][LUA_FUNCTIONS_NAMELEN];

void init_lua(char *plugin_file_name)
{
    int32_t err=0;
    
    strcpy(lua_file,plugin_file_name);

    /* initialize Lua */
    L = luaL_newstate();

    /* load Lua base libraries */
    luaL_openlibs(L);

    /* load the script */
    err = luaL_dofile(L, lua_file);
    if (err != 0)
    {
        LOG_PRINT(LOG_ERROR,"Erro while initialising Lua : %s \n",lua_tostring(L, -1));
        exit(-1);
    }
    else
        LOG_PRINT(LOG_INFO,"Lua successfully initialised : file %s opened. \n",lua_file);

}

// type is 0 for energy, 1 for gradient
void register_lua_function(char *plugin_function_name, LUA_FUNCTION_TYPE type)
{
    LOG_PRINT(LOG_INFO,"Registering Lua function : %s \n",plugin_function_name);
    
    if(type == POTENTIAL)
        strcpy(lua_function[POTENTIAL],plugin_function_name);
    else if(type == GRADIENT)
        strcpy(lua_function[GRADIENT],plugin_function_name);
    else
    {
        LOG_PRINT(LOG_ERROR,"Unknown type %d for register_lua_function.\n",type);
        exit(-1);
    }
}

void end_lua()
{
    if(is_lua_used)
    {
        /* cleanup Lua */
        lua_close(L);
    }
}

/*
 * This interface calls a lua script evaluating a Lennard Jobes like potential, pair by pair
 * see plugins/lj_n_m.lua
 */
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
                lua_getglobal(L, lua_function[POTENTIAL]);

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

                /* call the function with 10 arguments, return 1 result */
                lua_call(L, 10, 1);

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
                lua_getglobal(L, lua_function[0]);

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

                /* call the function with 10 arguments, return 1 result */
                lua_call(L,10,1);

                /* get the result */
                energy += (double)lua_tonumber(L, -1);
                lua_pop(L,1);
            }
        }
    } // end if-else on candidate

    LOG_PRINT(LOG_DEBUG,"LUA pair potential : %lf\n",energy);

    return energy;
}

void get_lua_DV(ATOM at[], DATA *dat, double fx[], double fy[], double fz[])
{
    uint32_t i=0 , j=0 ;
    double dx1,dy1,dz1;
    double dx2,dy2,dz2;
    
    for (i=0 ; i < dat->natom ; i++ )
    {
        fx[i] = 0.0 ;
        fy[i] = 0.0 ;
        fz[i] = 0.0 ;
        
        dx1=at[i].x;
        dy1=at[i].y;
        dz1=at[i].z;
        
        for (j=0 ; j < dat->natom ; j++ )
        {
            if (i==j) continue ;
            
            dx2=at[j].x;
            dy2=at[j].y;
            dz2=at[j].z;
            
            lua_getglobal(L, lua_function[GRADIENT]);

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

            /* call the function with 10 arguments, return 3 results */
            lua_call(L,10,3);
            
            fx[i] += (double)lua_tonumber(L, -1);
            fy[i] += (double)lua_tonumber(L, -2);
            fz[i] += (double)lua_tonumber(L, -3);
            lua_pop(L,3);
        }
    }
    
    LOG_PRINT(LOG_DEBUG,"LUA pair gradient done.\n");
}

/*
 * This interface calls a lua script evaluating a Lennard Jobes like potential
 * This version directly sends the ATOM structure, the Lua side requires LUAJIT/FFI
 * see plugins/lj_n_m_ffi.lua
 * 
 * This is the recommended way of calling lua script
 */
double get_lua_V_ffi(ATOM at[], DATA *dat, int32_t candidate)
{
    double energy=0.0;

    lua_getglobal(L, lua_function[POTENTIAL]);

    lua_pushinteger(L, dat->natom);
    lua_pushlightuserdata(L, (void*) at);
    lua_pushinteger(L, candidate);
    lua_call(L,3,1);

    /* get the result */
    energy += (double)lua_tonumber(L, -1);
    lua_pop(L, 1);

    LOG_PRINT(LOG_DEBUG,"LUA FFI potential : %lf\n",energy);

    return energy;
}

/*
 * This is the recommended way of calling lua script 
 */
void get_lua_DV_ffi(ATOM at[], DATA *dat, double fx[], double fy[], double fz[])
{
    lua_getglobal(L, lua_function[GRADIENT]);
    
//     function lj_dv_n_m_ffi(natom, at_list, fx, fy, fz)
    lua_pushinteger(L, dat->natom);
    lua_pushlightuserdata(L, (void*) at);
    lua_pushlightuserdata(L, (void*) fx);
    lua_pushlightuserdata(L, (void*) fy);
    lua_pushlightuserdata(L, (void*) fz);
    
    lua_call(L,5,0);
    
    LOG_PRINT(LOG_DEBUG,"LUA FFI gradient done.\n");
}

#endif //LUA_PLUGINS
