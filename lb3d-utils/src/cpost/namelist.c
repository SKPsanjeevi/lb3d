/* Code to parse a reduced subset of the Fortran namelist file format.
 *
 * This is ick on a stick.
 *
 * The grammar used is:
 *	namelist := ( w* group )*
 *	group:= w* & w* identifier (w* pair)* w* /
 *	pair := identifier w* = w* value
 *	value:= quoted-v | identifier
 *	quoted-v:= ` string `
 *	identifier:= (continuous finite-length series of non-whitespace characters,
 *			not including '=', '\'' )
 *	string := ( continuous finite-length series of non-whitespace characters)
 *	w* := (whitespace character, including CR)
 *
 * Basically, should cope with anything except arrays.
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#include "namelist.h"

/* Take a char pointer, advance it until either it points at non-whitespace,
 * or it points at end-of-string.
 */

char *skip_whitespace (char *p)
{
	while ( (*p!=0) && (
			   (*p==' ')
			|| (*p=='\t')
			|| (*p=='\r')
			|| (*p=='\n')
			))
			p++;
	return p;
}

char *skip_to_whitespace (char *p)
{
	while (	   (*p!=0)
		&& (*p!=' ')
		&& (*p!='\t')
		&& (*p!='\r')
		&& (*p!='\n')
		)
			p++;
	return p;
}

/* Take a char pointer, advance it until it points at either end-of-string, 
 * or whitespace, or an equals.
 */

char *nmlist_scan_identifier(char *p)
{
	while (    (*p!=0) 
		&& (*p!=' ')
		&& (*p!='\r')
		&& (*p!='\n')
		&& (*p!='\t')
		&& (*p!='=')
		&& (*p!='/')
		)
		p++;
	return p;
}

char *nmlist_copy_identifier(char *p,char **dest)
{
	char *dp;
	dp=*dest;
	while (    (*p!=0) 
		&& (*p!=' ')
		&& (*p!='\r')
		&& (*p!='\n')
		&& (*p!='\t')
		&& (*p!='=')
		&& (*p!='/')
		&& (*p!='\'')
		)
		*dp++=*p++;
	*dest=dp;
	return p;
}

char *nmlist_scan_value(char *buffer)
{
	char *p;
	p=buffer;
	if ('\''==*p) {
		p++;
		while ( (*p!=0) && (*p!='\'') ) {
			p++;
		}
		if (*p!='\'') {
			fprintf(stderr,"Unterminated quote\n");
			return NULL;
		}
		p++;
		return p;
	} else {
		return nmlist_scan_identifier(p);
	}
}

char *nmlist_copy_value(char *buffer, char **dest)
{
	char *p;
	char *dp;
	dp=*dest;
	p=buffer;
	if ('\''==*p) {
		p++;
		while ( (*p!=0) && (*p!='\'') ) {
			*dp++ = *p++;
		}
		if (*p!='\'') {
			fprintf(stderr,"Unterminated quote even after parsing!\n");
			return NULL;
		}
		p++;
		*dest=dp;
		return p;
	} else {
		return nmlist_copy_identifier(p,dest);
	}
}

/* Scan the pair in the given buffer.
 * If it parses, place the number of bytes of storage required in
 * pairsize, and return a pointer to the character after the pair.
 * Othersize, return NULL.
 */

char *nmlist_scan_pair(unsigned int *pairsize, char *buffer)
{
	unsigned int size;
	char *p,*q;

	/* 	identifier	*/

	p=nmlist_scan_identifier(buffer);
	if (buffer==p) {
		fprintf(stderr,"Corrupted pair - bad key\n");
		return NULL;
	}

	size = p-buffer + 1;

	/*	w*	*/

	p=skip_whitespace(p);

	/*	=	*/

	if ('='!=*p) {
		fprintf(stderr,"Corrupted pair - missing '='\n");
		fflush(stderr);
		printf("p=<%s>\n",p);
		return NULL;
	}
	p++;

	/*	w*	*/

	p=skip_whitespace(p);

	/*	value	*/

	if (NULL==(q=nmlist_scan_value(p))) {
		fprintf(stderr,"Corrupted pair - bad list\n");
		return NULL;
	}

	if ('\''==*p) {
		*pairsize = size + (q-p) -1; /* Don't allocate space for quotes. */
	} else {
		*pairsize = size + (q-p) +1;
	}

	return q;

}

char *nmlist_copy_pair(char **dest, nmlist_pair *pair, unsigned int *pairsize, char *buffer)
{
	char *p,*q,*dp;


	/* 	identifier	*/

	pair->key=*dest; /* Fill in the key pointer */
	p=nmlist_copy_identifier(buffer,dest);

	dp=*dest;
	*dp++=0;		/* Null-terminate the key */
	*dest=dp;

	/*	w*	*/

	p=skip_whitespace(p);

	/*	=	*/

	if ('='!=*p) {
		fprintf(stderr,"Corrupted pair after parse! - missing '='\n");
		return NULL;
	}
	p++;

	/*	w*	*/

	p=skip_whitespace(p);

	/*	value	*/

	pair->value = *dest;	/* Fill in the value pointer */
	if (NULL==(q=nmlist_copy_value(p,dest))) {
		fprintf(stderr,"Corrupted pair after parse! - bad list\n");
		return NULL;
	}
	dp=*dest;
	*dp++=0;			/* Null-terminate the value */
	*dest=dp;

	return q;

}

/* Scan through the supplied buffer,
 * returning either the number of bytes of storage which must
 * be allocated to hold all required information,
 * or <0 to indicate an error.
 * If successful, fills in group->npairs as well.
 */

int nmlist_scan_group (nmlist_group *group, char *buffer)
{
	unsigned int size=0;
	unsigned int pairsize,npairs=0;
	char *p,*q;

	/*	w*	*/

	p=skip_whitespace(buffer);

	/*	&	*/

	if ('&' != *p++) {
		return -1;
	}

	/*	w*	*/

	p=skip_whitespace(p);

	/*	identifier	*/

	q=nmlist_scan_identifier(p);

	if (q==p) {
		fprintf(stderr,"Group truncated before identifier\n");
		return -1;
	}

	size+=(q-p) +1;	/* Leave space for NULL */
	p=q;

	/*	w*	*/

	p=skip_whitespace(p);

	/*	pair*	*/

	while ( (*p!=0) && (*p!='/') ) {
		if (NULL==(q=nmlist_scan_pair(&pairsize,p))) {
			fprintf(stderr,"Error scanning pair\n");
			return -1;
		}
		npairs++;

		size += pairsize;
	/*	w*	*/
		p=skip_whitespace(q);
	}

	if ('/'!=*p) {
		fprintf(stderr,"Group missing trailing '/'\n");
		return -1;
	}

	group->npairs=npairs;
	return size;


}

char *nmlist_parse_group(nmlist_group *group, char *buffer)
{
	int siz;
	unsigned int pairsize,i;
	char *p,*q;
	char *text,*tptr;

	if ( (siz=nmlist_scan_group(group,buffer)) <0) {
		fprintf(stderr,"Error : group did not parse.\n");
		return NULL;
	}

	/************* Now assume that the group parses... ****************/
	/* ( This is ick. ) */

	if (NULL==(text=malloc((size_t)siz))) {
		perror("malloc()");
		printf("size=%d\n",siz);
		return NULL;
	}
	tptr=text;

	if (NULL==(group->pairs=malloc(sizeof(nmlist_pair)*(group->npairs)))) {
		perror("malloc()");
		free(text);
		return NULL;
	}

	/*	w*	*/

	p=skip_whitespace(buffer);

	/*	&	*/

	if ('&' != *p++) {
		fprintf(stderr,"Group suddenly does not parse!\n");
		free(text);
		free(group->pairs);
		return NULL;
	}

	/*	w*	*/

	p=skip_whitespace(p);

	/*	identifier	*/

	p=nmlist_copy_identifier(p,&tptr);
	*tptr++=0; /* Null-terminate group name */
	group->name=text;

	/*	w*	*/

	p=skip_whitespace(p);

	/*	pair*	*/

	for (i=0;i<(group->npairs) ; i++) {
		if (NULL==
			(q=nmlist_copy_pair(&tptr,&(group->pairs[i]),&pairsize,p))
			) {
			fprintf(stderr,"Error scanning pair AFTER PARSE!\n");
			free(text);
			free(group->pairs);
			return NULL;
		}

	/*	w*	*/
		p=skip_whitespace(q);
	}

	if ('/'!=*p) {
		fprintf(stderr,"Group missing trailing '/' AFTER PARSE!\n");
		free(text);
		free(group->pairs);
		return NULL;
	}
	p++;

	if (siz!=(tptr-text)) {
		printf("%cDANGEROUS ERROR! Memory accounting error of %ld bytes!\n",
			7,siz-(tptr-text));
		free(text);
		free(group->pairs);
		return NULL;
	}
	return p;

}

void nmlist_dump_group(nmlist_group *group)
{
	int i;
	printf("\t& %s\n",group->name);
	for (i=0;i<group->npairs;i++) {
		printf("\t\t%s = %s\n",
			group->pairs[i].key,
			group->pairs[i].value);
	}
	printf("/\n");

}

void nmlist_dealloc_group(nmlist_group *group)
{
	free(group->name);
	free(group->pairs);
}


/*
 * Open the specified file, and parse it as a Fortran namelist according
 * to the subset grammar specified above.
 *
 * Returns an array of struct nmlist_group_s; the last one has NULL name.
 * Each such struct has a pointer, `pairs', which points to an array of
 * char pointers, each element of which points to a null-term `key' string,
 * which is immediately followed by the null-term `value' string.
 *
 * Returns -1 on error.
 */ 

int nmlist_parse_file(nmlist *nl, char *filename)
{
	FILE *file;
	struct stat statbuf;
	off_t filesize;
	char *buffer; /* Contains the file */
	nmlist_group	*groups=NULL,*tmp;
	unsigned int	ngroups=0,i;
	char *p;

	if (0!=stat(filename,&statbuf)) {
		perror("stat()");
		return -1;
	}

	filesize = statbuf.st_size ;

	if (NULL==(buffer=malloc(filesize+1))) {
		perror("malloc()");
		return -1;
	}

	if (NULL==(file=fopen(filename,"r"))) {
		perror("fopen()");
		free(buffer);
		return -1;
	}

	if (-1==(fread(buffer,1,filesize,file))) {
		perror("fread()");
		free(buffer);
		return -1;
	}

	buffer[filesize]=0; /* Ensure null-terminated. */

	fclose(file);

	p=skip_whitespace(buffer);
	while (0!=*p) {
		ngroups++;
		/* Extend the groups list */
		if (NULL==(tmp=realloc(groups,sizeof(nmlist_group)*ngroups))) {
			perror("realloc()");
			for (i=0;i<(ngroups-1);i++) {
				nmlist_dealloc_group(&(groups[i]));
			}
			return -1;
		}
		groups=tmp;
		/* groups contains enough space for another element */

		p=nmlist_parse_group(&(groups[ngroups-1]),p);
		if (NULL==p) {
			fprintf(stderr,"Error parsing namelist.\n");
			for (i=0;i<(ngroups-1);i++) {
				nmlist_dealloc_group(&(groups[i]));
			}
			return -1;
		}
		p=skip_whitespace(p);
	}

	/* File successfully parsed! */

	free(buffer);
	nl->ngroups=ngroups;
	nl->groups=groups;
	return 0;
	
}

void nmlist_dealloc(nmlist *nl)
{
	int i;
	for (i=0;i<nl->ngroups;i++) {
		nmlist_dealloc_group(&(nl->groups[i]));
	}
	free(nl->groups);
}

char *nmlist_lookup(nmlist *nl, char *groupname, char *keyname)
{
	int i,j;
	nmlist_group *group;
	nmlist_pair *pair;

	i=nl->ngroups-1;

	while (i>=0) {
		group=&(nl->groups[i]);
		if (0==strcmp(groupname,group->name)) {
			/* Found a group with the right name */
			j=group->npairs-1;
			while (j>=0) {
				pair=&(group->pairs[j]);
				if (0==strcmp(keyname,pair->key)) {
					return pair->value;
				}
				j--;
			}
			/* Key not found in group - look for another group
			 * with the same name.
			 */
		}
				
		i--;
	}

	return NULL;
}

/* Look up the given key, and
 * write its value into the given double.
 * Returns zero on success, else error.
 */

int nmlist_lookup_double(nmlist *nl, char *groupname, char *keyname,double *x)
{
	char *value;
	if (NULL==(value=nmlist_lookup(nl,groupname,keyname))) {
		return -1;
	}

	*x=atof(value);
	return 0;
}

int nmlist_lookup_int(nmlist *nl, char *groupname, char *keyname,int *x)
{
	char *value;
	if (NULL==(value=nmlist_lookup(nl,groupname,keyname))) {
		return -1;
	}

	*x=atoi(value);
	return 0;
}

/*
 * Look up the given key: return 0 if it's ".false.", 1 if it's ".true.",
 * and -1 on error.
 */

int nmlist_lookup_bool(nmlist *nl, char *groupname, char *keyname)
{
	char *value;
	if (NULL==(value=nmlist_lookup(nl,groupname,keyname))) {
		return -1;
	}
	if (0==strcasecmp(".true.",value)) {
		return 1;
	}
	if (0==strcasecmp(".false.",value)) {
		return 0;
	}
	return -1;
}



nmlist_group *nmlist_get_group(nmlist *nl, char *groupname)
{
	int i;
	nmlist_group *group;

	i=nl->ngroups-1;

	while (i>=0) {
		group=&(nl->groups[i]);
		if (0==strcmp(groupname,group->name)) {
			return group;
		}
	}

	return NULL;
}

