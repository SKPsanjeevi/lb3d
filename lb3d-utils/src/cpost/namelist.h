#ifndef _NAMELIST_H
#define _NAMELIST_H

struct nmlist_pair_s {
	char *key;
	char *value;
};

typedef struct nmlist_pair_s nmlist_pair;

struct nmlist_group_s {
	char *name;	/* Group name */
	unsigned int	npairs; /* Number of pairs */
	nmlist_pair *pairs; /* Array of struct pairs */
};

typedef struct nmlist_group_s nmlist_group;

struct nmlist_s {
	unsigned int	ngroups;
	nmlist_group	*groups;
};

typedef struct nmlist_s nmlist;

char *skip_whitespace (char *);
char *skip_to_whitespace (char *);
char *nmlist_parse_group(nmlist_group *, char *);
void nmlist_dump_group(nmlist_group *);
void nmlist_dealloc_group(nmlist_group *);
int nmlist_parse_file(nmlist *nl, char *);
void nmlist_dealloc(nmlist *);
char *nmlist_lookup(nmlist *, char *, char *);
nmlist_group *nmlist_get_group(nmlist *, char *);
int nmlist_lookup_double(nmlist *, char *, char *,double *);
int nmlist_lookup_int(nmlist *, char *, char *,int *);
int nmlist_lookup_bool(nmlist *, char *, char *);
#endif /* _NAMELIST_H */
