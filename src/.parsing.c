#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>

#include "global.h"
#include "io.h"
#include "ener.h"
#include "tools.h"

static int lj_size = 0 ;

ATOM* parse_from_file(char fname[], DATA *dat, SPDAT *spdat)
{
	int i=0;
	
	ATOM *at=NULL;
	LJPARAMS *ljpars=NULL;
	
	char buff1[1024]="", *buff2=NULL, *buff3=NULL ;
	
	FILE *ifile=NULL;
	ifile=fopen(fname,"r");
	
	if (ifile==NULL)
	{
		printf("Error while opening the file '%s'\n",fname);
		exit(-1);
	}
	
	while(fgets(buff1,1024,ifile)!=NULL)
	{
		buff2=strtok(buff1," \n\t");
		
		while (buff2 != NULL)
		{
			buff3=strtok(NULL," \n\t");
// 			printf("%s %s\n",buff2,buff3);
			
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
					spdat->meps = atoi(meps);
					
					neps=strtok(NULL," \n\t");
					neps=strtok(NULL," \n\t");
					spdat->neps = atoi(neps);
					
					sprintf(dat->method,"%s",buff3);
				}
				else
					fprintf(stdout,"Warning : %s %s is unknown. Should be METROP or SPAV.\n",buff2,buff3);
			}
			else if (!strcasecmp(buff2,"POTENTIAL"))
			{
				if (strcasecmp(buff3,"LJ") && strcasecmp(buff3,"AZIZ"))
					fprintf(stdout,"Warning : %s %s is unknown. Should be LJ or AZIZ.\n",buff2,buff3);
				else
				{
					if (!strcasecmp(buff3,"AZIZ"))
						get_ENER = &(get_AZIZ_V);
					else if (!strcasecmp(buff3,"LJ"))
					{
						get_ENER = &(get_LJ_V);
						get_DV = &(get_LJ_DV);
					}
				}
			}
			else if (!strcasecmp(buff2,"UNITS"))
			{
				if (!strcasecmp(buff3,"REDUCED"))
					charmm_units=0;
				else if (!strcasecmp(buff3,"CHARMM"))
					charmm_units=1;
				else
					fprintf(stdout,"Warning : %s %s is unknown. Must be CHARMM or REDUCED.\n",buff2,buff3);
			}
			else if (!strcasecmp(buff2,"SAVE"))
            {
				if (!strcasecmp(buff3,"ENER"))
				{
					char *title=NULL , *each=NULL;
					title = strtok(NULL," \n\t\'");
					sprintf(io.etitle,"%s",title);
					each = strtok(NULL," \n\t"); //junk
					each = strtok(NULL," \n\t");
					io.esave = atoi(each);
				}
				else if (!strcasecmp(buff3,"COOR"))
				{
					char *what=NULL , *type=NULL , *title=NULL ;
					what = strtok(NULL," \n\t");
					if (!strcasecmp(what,"FIRST"))
					{
						type = strtok(NULL," \n\t");
						title = strtok(NULL," \n\t\'");
						sprintf(io.crdtitle,"%s",title);
					}
					else if (!strcasecmp(what,"TRAJ"))
					{
						char *each=NULL;
						type = strtok(NULL," \n\t");
						title = strtok(NULL," \n\t\'");
						sprintf(io.trajtitle,"%s",title);
						each = strtok(NULL," \n\t"); //junk
						each = strtok(NULL," \n\t");
						io.trsave = atoi(each);
					}
				}
            }
			else if (!strcasecmp(buff2,"NATOMS"))
			{
				dat->natom=atoi(buff3);
				at = calloc(dat->natom,sizeof *at);
				build_cluster(at,dat,0,dat->natom,-1);	//initialise
			}
			else if (!strcasecmp(buff2,"TEMP"))
				dat->T=atof(buff3);
			else if (!strcasecmp(buff2,"NSTEPS"))
				dat->nsteps=atoi(buff3);
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
					dat->d_max_when=atoi(each);
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
				int i=0,j=0,k=0,l=0;
				char *from=buff3 , *to=NULL , *type=NULL , 	*coor=NULL;
				
				to=strtok(NULL," \n\t");
				to=strtok(NULL," \n\t");
				type=strtok(NULL," \n\t");
				coor=strtok(NULL," \n\t");
				coor=strtok(NULL," \n\t");
				
				j=atoi(from)-1;
				k=atoi(to);
				
				if (k<=0 || k>dat->natom)
					k=dat->natom;
				
				for(i=j;i<k;i++)
				{
					sprintf(at[i].sym,"%s",type);
					for(l=0;l<lj_size;l++)
					{
						if (!strcasecmp(ljpars[l].sym,at[i].sym))
						{
							at[i].ljp.eps=ljpars[l].eps;
							at[i].ljp.sig=ljpars[l].sig;
							break;
						}
					}
				}
				/*
				if(!strcasecmp(coor,"RANDOM"))
					build_cluster(at,dat,j,k,1);
				else if (!strcasecmp(coor,"ZERO"))
					build_cluster(at,dat,j,k,0);
				else
				{
					//random is default
					build_cluster(at,dat,j,k,1);
				}
				*/
				
			}
			buff2=strtok(NULL," \n\t");
		}	
	}
	
	fclose(ifile);
	free(ljpars);
    
    build_cluster(at,dat,0,dat->natom,1);
	
	return at;
}
