#!/bin/bash

## Pass in a tar
tar_file=$1

## Clean up all
rm -f ../c/*.c ../c/*~
rm -f ../netlib/*.f ../netlib/*~

## 
if [[ "$tar_file" == *.tgz || "$tar_file" == *.tar.gz ]]
then
    (cd ../netlib; tar -xzvf $tar_file; rm -f ._*.f; cat support.list | xargs rm -f; cat install.list | xargs rm -f)
elif [[ "$tar_file" == *.tar ]]
then
    ## xerbla.f and any .f using f90 intrinsic (maxloc) are not imported
    (cd ../netlib; tar -xvf $tar_file; rm -f ._*.f; cat support.list | xargs rm -f; cat install.list | xargs rm -f)
else
    echo "Please pass a netlibs tar file as arg 1"
    exit
fi

## Polish fortran files.
##  Syntax error:  CHARACTER(1) --> CHARACTER
##  F90 Support : CYCLE --> CALL F90_CYCLE (temporary)
##  F90 Support : EXIT --> CALL F90_EXIT (temporary)
##  etc.
./preprosess.sh


## Run f2c 
f2c -A -R -a ../netlib/*.f ../netlib/*.F -d../c/
rm -f ._* ../netlib/*.f ../netlib/*.F

## Remove ftnlen 
##  ftnlen arguments should be carefully removed to preserve correct fortran-level interfaces
##  most f2c'ed LAPACK routines do not need ftnlen arguments e.g., CHARACTER(1) do not need.
##  the script will replace the arguments and function definitions.
./remove_ftnlen.sh

## Change temporary tweaks
##  f90_cycle__() --> continue
##  f90_exit__() --> break
##  etc.
./postprocess.sh

## Remove comma expressions
##  f2c generates comma expressions for max and min macros. Such an expression should be properly removed. 
./remove_comma_expr.sh
# ./remove_comma_abs.sh

## Generate a prototypes header
## ./headergen.sh


## Apply artistic style
./astyle.sh


## Disable xblas-related code so shared libraries link properly.
./disable_xblas.sh


echo " "
echo " Note that string-related f2c functions should not be used here (e.g., s_copy, s_cmp)."
echo " It is a total mess when an actual string is mixed in C and FORTRAN.                  "
echo " s_cat is an exception and it is modified from the original definition in f2clib.     "
echo " " 
echo "                       Do not link against -lf2c.                                     "
