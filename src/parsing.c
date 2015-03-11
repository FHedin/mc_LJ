/**
 * \file parsing.c
 *
 * \brief File containing functions for parsing the input file
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "global.h"
#include "io.h"
#include "ener.h"
#include "tools.h"
#include "logger.h"
#include "plugins_lua.h"

///the array of LJ-params size
static uint32_t lj_size = 0 ;

/**
 * @brief his function parses the input file, fills many fields of the DATA and SPDAT structures,
 * and allocates the ATOM list.
 * 
 * @note strcasecmp(...) is used so the input file is case insensitive
 * 
 * @param fname Path to the input file to open
 * @param dat Common data
 * @param spdat Common data for Spatial Averaging
 * @param at The atom list
 */
void parse_from_file(char fname[], DATA *dat, SPDAT *spdat, ATOM **at)
{
    LJPARAMS *ljpars=NULL;

    char buff1[FILENAME_MAX]="", *buff2=NULL, *buff3=NULL ;

    FILE *ifile=NULL;
    ifile=fopen(fname,"r");

    if (ifile==NULL)
    {
        LOG_PRINT(LOG_ERROR,"Error while opening the file '%s'\n",fname);
        exit(-1);
    }

    ///iterate over each line of the text file
    while(fgets(buff1,FILENAME_MAX,ifile)!=NULL)
    {
        // skip comment line, but print it to LOG_INFO , may be useful for debugging a bad input file
        if (buff1[0]=='#')
        {
            LOG_PRINT(LOG_INFO,"Skipping line %s",buff1);
            continue;
        }

        buff2=strtok(buff1," \n\t");

        while (buff2 != NULL)
        {
            buff3=strtok(NULL," \n\t");

            ///to know which MC method we use
            if (!strcasecmp(buff2,"METHOD"))
            {
                ///no extra parameter for Metropolis currently
                if (!strcasecmp(buff3,"METROP"))
                    sprintf(dat->method,"%s",buff3);
                ///for spatial averaging extra parameters are required
                ///see literature for more details
                else if (!strcasecmp(buff3,"SPAV"))
                {
                    char *weps=NULL , *meps=NULL , *neps=NULL;

                    ///weps is the width of the gaussian distribution
                    weps=strtok(NULL," \n\t");
                    weps=strtok(NULL," \n\t");
                    spdat->weps = atof(weps);

                    ///meps is the number of sets
                    meps=strtok(NULL," \n\t");
                    meps=strtok(NULL," \n\t");
                    spdat->meps = (uint32_t) atoi(meps);

                    ///neps the number of replicated structures
                    neps=strtok(NULL," \n\t");
                    neps=strtok(NULL," \n\t");
                    spdat->neps = (uint32_t) atoi(neps);

                    sprintf(dat->method,"%s",buff3);
                }
                else
                {
                    LOG_PRINT(LOG_WARNING,"%s %s is unknown. Should be METROP or SPAV.\n",buff2,buff3);
                }
            }
            ///get type of potential we plan to use
            else if (!strcasecmp(buff2,"POTENTIAL"))
            {
                ///available : hard coded LJ potential, hard coded Aziz potential, and user defined potential read from LUA script
                if (strcasecmp(buff3,"LJ") && strcasecmp(buff3,"AZIZ") && strcasecmp(buff3,"PLUGIN"))
                {
                    LOG_PRINT(LOG_WARNING,"%s %s is unknown. Should be LJ or AZIZ or PLUGIN.\n",buff2,buff3);
                }
                else
                {
                    ///the user of pointers to functions avoids the use of if(...) in energy functions so code is faster
                    if (!strcasecmp(buff3,"AZIZ"))
                        get_ENER = &(get_AZIZ_V);
                    else if (!strcasecmp(buff3,"LJ"))
                    {
                        get_ENER = &(get_LJ_V);
                        get_DV = &(get_LJ_DV);
                    }
#ifdef LUA_PLUGINS
                    /**
                     * for the lua plugin potentials the two mode are 
                     *  - PAIR where the user will write code returning the v and dv contribution for 2 atoms only
                     *  - FFI where the user will have access from the LUA code to the data structures defined in global.h ; only available if your
                     *      LUA installation come with the FFI package, included by default with luajit (see online)
                     * 
                     */
                    else if (!strcasecmp(buff3,"PLUGIN"))
                    {
                        char *buff4=NULL , *buff5=NULL , *buff6=NULL , *buff7=NULL;
//                         PLUGIN_TYPE plug_type;
                        buff4=strtok(NULL," \n\t");
                        buff5=strtok(NULL," \n\t");
                        buff6=strtok(NULL," \n\t");
                        buff7=strtok(NULL," \n\t");
                        
                        if (!strcasecmp(buff4,"PAIR"))
                        {
                            lua_plugin_type = PAIR;
                            get_ENER = &(get_lua_V);
                            get_DV = &(get_lua_DV);
                        }
                        else if (!strcasecmp(buff4,"FFI"))
                        {
                            lua_plugin_type = FFI;
                            get_ENER = &(get_lua_V_ffi);
                            get_DV = &(get_lua_DV_ffi);
                        }
                        else
                        {
                            LOG_PRINT(LOG_ERROR,"Plugin type %s is unknown. Must be PAIR or FFI.\n",buff4);
                            exit(-1);
                        }
                        
                        /// open lua file and register function name
                        init_lua(buff5);
                        register_lua_function(buff6,POTENTIAL);
                        register_lua_function(buff7,GRADIENT);
                    }
#endif
                }
            }
            /// units to use during simulation
            /// reduced units recommended, CHARMM units are experimental
            else if (!strcasecmp(buff2,"UNITS"))
            {
                if (!strcasecmp(buff3,"REDUCED"))
                    charmm_units=0;
                else if (!strcasecmp(buff3,"CHARMM"))
                    charmm_units=1;
                else
                {
                    LOG_PRINT(LOG_WARNING,"%s %s is unknown. Must be CHARMM or REDUCED.\n",buff2,buff3);
                }
            }
            /// section where saving of energy, coordinates and trajectory is handled
            else if (!strcasecmp(buff2,"SAVE"))
            {
                ///energy saving
                if (!strcasecmp(buff3,"ENER"))
                {
                    char *title=NULL , *each=NULL;
                    title = strtok(NULL," \n\t\'");
                    sprintf(io.etitle,"%s",title);
                    each = strtok(NULL," \n\t");
                    each = strtok(NULL," \n\t");
                    io.esave = (uint32_t) atoi(each);
                }
                ///coordinates saving
                else if (!strcasecmp(buff3,"COOR"))
                {
                    char *what=NULL , *type=NULL , *title=NULL ;
                    what = strtok(NULL," \n\t");
                    ///for saving initial coordinates i.e. before simulation starts
                    if (!strcasecmp(what,"FIRST"))
                    {
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.crdtitle_first,"%s",title);
                    }
                    ///for saving last coordinates i.e. after simulation ends
                    else if (!strcasecmp(what,"LAST"))
                    {
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.crdtitle_last,"%s",title);
                    }
                    else if (!strcasecmp(what,"TRAJ"))
                    {
                        char *each=NULL;
                        ///only dcd type for the moment so useless variable
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.trajtitle,"%s",title);
                        each = strtok(NULL," \n\t"); //junk
                        each = strtok(NULL," \n\t");
                        io.trsave = (uint32_t) atoi(each);
                    }
                }
            }
            /// define number of atoms 
            else if (!strcasecmp(buff2,"NATOMS"))
            {
                dat->natom = (uint32_t) atoi(buff3);
                *at = malloc(dat->natom*sizeof(ATOM));
                build_cluster(*at,dat,0,dat->natom,-1);	///< initialise the cluster with atoms at infinity initially
            }
            /// define temperature
            else if (!strcasecmp(buff2,"TEMP"))
                dat->T = atof(buff3);
            /// define number of steps as 64 bits integer
            else if (!strcasecmp(buff2,"NSTEPS"))
                dat->nsteps = (uint64_t) strtoull(buff3,NULL,0);    //atol(buff3);
            /// defines the dmax and if we want it fixed or automatically updated
            else if (!strcasecmp(buff2,"DMAX"))
            {
                char *mode=NULL , *each=NULL , *target=NULL;
                dat->d_max=atof(buff3);
                mode=strtok(NULL," \n\t");
                dat->d_max_when=0;
                dat->d_max_tgt=0.0;

                /// if we want dmax to be automatically updated. see global.h and tools.c for details
                if (!strcasecmp(mode,"UPDATE"))
                {
                    each=strtok(NULL," \n\t");
                    target=strtok(NULL," \n\t"); //junk
                    target=strtok(NULL," \n\t");
                    dat->d_max_when = (uint32_t) atoi(each);
                    dat->d_max_tgt=atof(target);
                }
            }
            /// define a list of LJ parameters
            else if (!strcasecmp(buff2,"LJPARAMS"))
            {
                char *type=buff3 , *epsi=NULL , *sigma=NULL ;

                /// the array grows each time a LJPARAMS section is detected
                ljpars=(LJPARAMS*)realloc(ljpars,(lj_size+1)*sizeof(LJPARAMS));

                sprintf(ljpars[lj_size].sym,"%s",type);

                epsi=strtok(NULL," \n\t");
                epsi=strtok(NULL," \n\t");
                ljpars[lj_size].eps=atof(epsi);

                sigma=strtok(NULL," \n\t");
                sigma=strtok(NULL," \n\t");
                ljpars[lj_size].sig=atof(sigma);

                lj_size++;
            }
            /// for building atom list manually
            else if (!strcasecmp(buff2,"ATOM"))
            {
                uint32_t i=0,j=0,k=0,l=0;
                char *from=buff3, *to=NULL, *type=NULL, *coor=NULL;
                
                to=strtok(NULL," \n\t");
                to=strtok(NULL," \n\t");
                type=strtok(NULL," \n\t");
                coor=strtok(NULL," \n\t");
                coor=strtok(NULL," \n\t");
                
                j = (uint32_t) atoi(from) - 1;
                if (!strcasecmp(to,"END"))
                    k = dat->natom;
                else
                    k = (uint32_t) atoi(to);
                
                if (k<=0 || k>dat->natom)
                    k=dat->natom;

                LOG_PRINT(LOG_INFO,"Building an atomic list from index %d to %d and of type  %s.\n",j,k-1,type);
                
                for(i=j; i<k; i++)
                {
                    sprintf((*at)[i].sym,"%s",type);
                    for(l=0; l<lj_size; l++)
                    {
                        if (!strcasecmp(ljpars[l].sym,(*at)[i].sym))
                        {
                            (*at)[i].ljp.eps=ljpars[l].eps;
                            (*at)[i].ljp.sig=ljpars[l].sig;
                            break;
                        }
                    }
                }
                ///randomly distribute atoms
                if(!strcasecmp(coor,"RANDOM"))
                {
                    build_cluster(*at,dat,j,k,1);
                }
                ///put all atoms at origin
                else if (!strcasecmp(coor,"ZERO"))
                {
                    build_cluster(*at,dat,j,k,0);
                }
                ///read starting cnfiguration from file
                else if (!strcasecmp(coor,"FILE"))
                {
                    // initial structure read from an xyz file
                    char *fi=NULL;
                    fi=strtok(NULL," \n\t'");

                    FILE *start=NULL;
                    start = fopen(fi,"r");
                    if (start==NULL)
                    {
                        LOG_PRINT(LOG_ERROR,"Error while opening initial structure file %s\n",fi);
                        exit(-1);
                    }

                    read_xyz(*at,dat,start);

                    fclose(start);
                }
                else
                {
                    ///random is default
                    build_cluster(*at,dat,j,k,1);
                }
            }

            buff2=strtok(NULL," \n\t");
        }
    }

    fclose(ifile);
    free(ljpars);
}
