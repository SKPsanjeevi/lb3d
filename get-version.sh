#!/bin/bash

function usage()
{
echo "Syntax is $0 [-m|-h|-c|-t]" >&2
echo  >&2
echo "  -m : Generate Makefile section."  >&2
echo "  -h : Generate code for Fortran header file."  >&2
echo "  -c : Generate code for C header file."  >&2
echo "  -t : Generate code for TeX \newcommand."  >&2
echo  >&2
exit -1
}

# If no flag is set, show usage and abort.
if [[ -z "$1" ]]
then
    usage
else

# Otherwise, proceed with finding versioning information.

# Run git diff to check for local changes. Return code 0 indicates no changes, 1 indicates changes are present.
# Exit code 129 is an error when not inside a Git repo. Use that code to avoid errors for the describe and log commands.
git diff --quiet 2> /dev/null
DIFFCODE=$?

if [[ "$DIFFCODE" -ne 129 ]]
then
    DESC=`git describe --always --tag` > /dev/null 2>&1
    BRANCH=`git rev-parse --abbrev-ref HEAD | tr '\n' ' '`
    if [[ "$DIFFCODE" -ne 0 ]]
    then
        LOCALCHANGES=" (with local changes)"
    fi
else
    DESC="exported"
    BRANCH="exported"
fi

SVNREV=`svnversion 2>/dev/null`

case "$1" in

    "-m")
# Section for the Makefile.
# Do NOT expand variables in this case - hence the double quotes around "EOF".
read -d '' MAKESECTION <<"EOF"

VTMP="lbe_version.tmp"
VH="lbe_version.h"

version:
	@./../get-version.sh -h > $(VTMP)
	@echo "#define LBE_FLAGS \\"$(DEFFLAGS) $(MAKEFFLAGS)\\"" >> $(VTMP)
	@echo "#define LBE_PLATFORM \\"$(PLATFORM)\\"" >> $(VTMP)
	@diff $(VTMP) $(VH) > /dev/null 2>&1 ; if [ "$$?" -eq 0 ] ; then rm -rf $(VTMP) ; else echo "Creating new $(VH)"; mv $(VTMP) $(VH); fi
EOF
	echo "$MAKESECTION"
	;;

    "-h")
# Code for Fortran header file
read -d '' FORTRANHEADER <<EOF
! Putting revision and flag information in some pragma for later use.

#define GIT_DESC \"$DESC\"
#define GIT_DIFFCODE $DIFFCODE
#define GIT_LOCALCHANGES \"$LOCALCHANGES\"

#define GIT_BRANCH \"$BRANCH\"

#define SVN_REVISION \"$SVNREV\"

EOF
	echo "$FORTRANHEADER"
	;;

    "-c")
# Code for C header file
read -d '' CHEADER <<EOF
// Putting revision and flag information in some pragma for later use.

#define GIT_DESC \"$DESC\"
#define GIT_DIFFCODE $DIFFCODE
#define GIT_LOCALCHANGES \"$LOCALCHANGES\"

#define GIT_BRANCH \"$BRANCH\"

#define SVN_REVISION \"$SVNREV\"

EOF
	echo "$CHEADER"
	;;

    "-t")
	# TeX newcommand
	echo "\newcommand{\revision}{exported}" | sed "s/exported/$DESC/"
	;;

    *)
	echo "Unknown option '$1'." >&2
	usage
esac

fi

