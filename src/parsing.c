/*
 * Copyright (c) 2014, Florent Hedin, Markus Meuwly, and the University of Basel
 * All rights reserved.
 *
 * The 3-clause BSD license is applied to this software.
 * see LICENSE.txt
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

static uint32_t lj_size = 0 ;

void parse_from_file(char fname[], DATA *dat, SPDAT *spdat, ATOM **at)
{
    LJPARAMS *ljpars=NULL;

    char buff1[FILENAME_MAX]="", *buff2=NULL, *buff3=NULL ;

    FILE *ifile=NULL;
    ifile=fopen(fname,"r");

    if (ifile==NULL)
    {
//         fprintf(stdout,"Error while opening the file '%s'\n",fname);
        LOG_PRINT(LOG_ERROR,"Error while opening the file '%s'\n",fname);
        exit(-1);
    }

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

            if (!strcasecmp(buff2,"METHOD"))
            {
                if (!strcasecmp(buff3,"METROP"))
                    sprintf(dat->method,"%s",buff3);
                else if (!strcasecmp(buff3,"SPAV"))
                {
                    char *weps=NULL , *meps=NULL , *neps=NULL;

                    weps=strtok(NULL," \n\t");
                    weps=strtok(NULL," \n\t");
                    spdat->weps = atof(weps);

                    meps=strtok(NULL," \n\t");
                    meps=strtok(NULL," \n\t");
                    spdat->meps = (uint32_t) atoi(meps);

                    neps=strtok(NULL," \n\t");
                    neps=strtok(NULL," \n\t");
                    spdat->neps = (uint32_t) atoi(neps);

                    sprintf(dat->method,"%s",buff3);
                }
                else
                {
                    LOG_PRINT(LOG_WARNING,"%s %s is unknown. Should be METROP or SPAV.\n",buff2,buff3);
//                     fprintf(stdout,"Warning : %s %s is unknown. Should be METROP or SPAV.\n",buff2,buff3);
                }
            }
            else if (!strcasecmp(buff2,"POTENTIAL"))
            {
                if (strcasecmp(buff3,"LJ") && strcasecmp(buff3,"AZIZ") && strcasecmp(buff3,"PLUGIN"))
                {
                    LOG_PRINT(LOG_WARNING,"%s %s is unknown. Should be LJ or AZIZ or PLUGIN.\n",buff2,buff3);
//                     fprintf(stdout,"Warning : %s %s is unknown. Should be LJ or AZIZ or PLUGIN.\n",buff2,buff3);
                }
                else
                {
                    if (!strcasecmp(buff3,"AZIZ"))
                        get_ENER = &(get_AZIZ_V);
                    else if (!strcasecmp(buff3,"LJ"))
                    {
                        get_ENER = &(get_LJ_V);
                        get_DV = &(get_LJ_DV);
                    }
#ifdef LUA_PLUGINS
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
                        
                        // open lua file and register function name
                        init_lua(buff5);
                        register_lua_function(buff6,POTENTIAL);
                        register_lua_function(buff7,GRADIENT);
                    }
#endif
                }
            }
            else if (!strcasecmp(buff2,"UNITS"))
            {
                if (!strcasecmp(buff3,"REDUCED"))
                    charmm_units=0;
                else if (!strcasecmp(buff3,"CHARMM"))
                    charmm_units=1;
                else
                {
                    LOG_PRINT(LOG_WARNING,"%s %s is unknown. Must be CHARMM or REDUCED.\n",buff2,buff3);
//                     fprintf(stdout,"Warning : %s %s is unknown. Must be CHARMM or REDUCED.\n",buff2,buff3);
                }
            }
            else if (!strcasecmp(buff2,"SAVE"))
            {
//                 fprintf(stderr,"Entering save\n");
//                 fflush(stderr);
                if (!strcasecmp(buff3,"ENER"))
                {
                    char *title=NULL , *each=NULL;
                    title = strtok(NULL," \n\t\'");
                    sprintf(io.etitle,"%s",title);
                    each = strtok(NULL," \n\t"); //junk
                    each = strtok(NULL," \n\t");
                    io.esave = (uint32_t) atoi(each);
//                     fprintf(stderr,"Save ener file\n");
//                     fflush(stderr);
                }
                else if (!strcasecmp(buff3,"COOR"))
                {
                    char *what=NULL , *type=NULL , *title=NULL ;
                    what = strtok(NULL," \n\t");
                    if (!strcasecmp(what,"FIRST"))
                    {
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.crdtitle_first,"%s",title);
//                         fprintf(stderr,"Save coor first\n");
//                         fflush(stderr);
                    }
                    else if (!strcasecmp(what,"LAST"))
                    {
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.crdtitle_last,"%s",title);
//                         fprintf(stderr,"Save coor last\n");
//                         fflush(stderr);
                    }
                    else if (!strcasecmp(what,"TRAJ"))
                    {
                        char *each=NULL;
                        type = strtok(NULL," \n\t");
                        title = strtok(NULL," \n\t\'");
                        sprintf(io.trajtitle,"%s",title);
                        each = strtok(NULL," \n\t"); //junk
                        each = strtok(NULL," \n\t");
                        io.trsave = (uint32_t) atoi(each);
//                         fprintf(stderr,"Save coor traj\n");
//                         fflush(stderr);
                    }
                }
            }
            else if (!strcasecmp(buff2,"NATOMS"))
            {
                dat->natom = (uint32_t) atoi(buff3);
                *at = malloc(dat->natom*sizeof(ATOM));
                build_cluster(*at,dat,0,dat->natom,-1);	//initialise
            }
            else if (!strcasecmp(buff2,"TEMP"))
                dat->T = atof(buff3);
            else if (!strcasecmp(buff2,"NSTEPS"))
                dat->nsteps = (uint64_t) strtoull(buff3,NULL,0);    //atol(buff3);
            else if (!strcasecmp(buff2,"DMAX"))
            {
                char *mode=NULL , *each=NULL , *target=NULL;
                dat->d_max=atof(buff3);
                mode=strtok(NULL," \n\t");
                dat->d_max_when=0;
                dat->d_max_tgt=0.0;

                if (!strcasecmp(mode,"UPDATE"))
                {
                    each=strtok(NULL," \n\t");
                    target=strtok(NULL," \n\t"); //junk
                    target=strtok(NULL," \n\t");
                    dat->d_max_when = (uint32_t) atoi(each);
                    dat->d_max_tgt=atof(target);
                }
            }
            else if (!strcasecmp(buff2,"LJPARAMS"))
            {
                char *type=buff3 , *epsi=NULL , *sigma=NULL ;

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
            else if (!strcasecmp(buff2,"ATOM"))
            {
                uint32_t i=0,j=0,k=0,l=0;
                char *from=buff3 , *to=NULL , *type=NULL , 	*coor=NULL;

//                 fprintf(stdout,"Building an atom list.\n");
                
                to=strtok(NULL," \n\t");
                to=strtok(NULL," \n\t");
                type=strtok(NULL," \n\t");
                coor=strtok(NULL," \n\t");
                coor=strtok(NULL," \n\t");

//                 fprintf(stdout,"Building an atom list done.\n");
//                 fflush(stdout);
                
                j = (uint32_t) atoi(from) - 1;
                k = (uint32_t) atoi(to);

//                 fprintf(stdout,"Building an atom list : from to : done.\n");
//                 fflush(stdout);
                
                if (k<=0 || k>dat->natom)
                    k=dat->natom;

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
                
//                 fprintf(stdout,"LJ params arrays filled with params\n");
//                 fflush(stdout);

                if(!strcasecmp(coor,"RANDOM"))
                {
                    build_cluster(*at,dat,j,k,1);
//                     fprintf(stdout,"Initial cluster build by random\n");
//                     fflush(stdout);
                }
                else if (!strcasecmp(coor,"ZERO"))
                {
                    build_cluster(*at,dat,j,k,0);
//                     fprintf(stdout,"Initial cluster build by zero\n");
//                     fflush(stdout);
                }
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
//                         fprintf(stdout,"Error while opening initial structure file %s\n",fi);
                        exit(-1);
                    }

                    read_xyz(*at,dat,start);
                    
//                     steepd_ini(*at,dat);

                    fclose(start);
                    
//                     fprintf(stdout,"Initial cluster build by file\n");
//                     fflush(stdout);
                }
                else
                {
                    //random is default
                    build_cluster(*at,dat,j,k,1);
//                     fprintf(stdout,"Initial cluster build by random\n");
//                     fflush(stdout);
                }
            }

            buff2=strtok(NULL," \n\t");
        }
    }

    fclose(ifile);
    free(ljpars);
}
