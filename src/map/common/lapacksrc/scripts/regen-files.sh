#!/bin/bash

## Pass in a tar
tar_file=$1

## Clean up all
rm -f ../c/*.c ../c/*~
rm -f ../netlib/*.[fF] ../netlib/*.[fF]90 ../netlib/*~
rm -rf ../netlib/lapack-*/

## 
if [[ "$tar_file" == *.tgz || "$tar_file" == *.tar.gz ]]
then
    (cd ../netlib; tar -xzf $tar_file ;
    cp lapack-*/SRC/*.[fF] .;
    cp lapack-*/SRC/*.[fF]90 .;
    rm -f ._*.f;
    cat support.list | xargs rm -f;
    cat install.list | xargs rm -f)
    full_tar=1
elif [[ "$tar_file" == *.tar ]]
then
    (cd ../netlib;
    tar -xf $tar_file;
    rm -f ._*.f;
    cat support.list | xargs rm -f;
    cat install.list | xargs rm -f)
    full_tar=0
else
    echo "Please pass a netlibs tar file as arg 1"
    exit
fi

## Move to the lapack dir
mkdir ../fortran/
mv ../netlib/*.[fF] ../fortran/;
mv ../netlib/*.[fF]90 ../fortran/;

## Updating other areas in libflame
if [[ "$full_tar" == 1 ]]
then
    ## Get the name of the directory
    new_test_dir=$(echo $tar_file | sed "s/\.tgz//" | sed "s/\.tar.gz//" | sed "s/\..\/netlib\///")

    ## Call the script that makes a new netlibs-test
    (cd ../../../../../netlib-test/;
    ./create_new_testdir.sh ../src/map/lapack2flamec/f2c/netlib/$new_test_dir $new_test_dir)

    ## Regen the f2c in src/flablas/f2c/
    (cd ../../../../flablas/f2c/;
    rm -f *.[fF] *.[fF]90 *.c;
    cp ../../map/lapack2flamec/f2c/netlib/$new_test_dir/BLAS/SRC/* .;
    ./regen-files.sh)
else
    echo " "
    echo " Skipping other installs because they are not in the tar "
    echo " "
fi

## Disable xblas-related code so shared libraries link properly.
./disable_xblas.sh

## Cleanup
rm -f ._* ../netlib/*.[fF] ../netlib/*.[fF]90
rm -rf ../netlib/lapack-*/

echo " "
echo " Note that string-related f2c functions should not be used here (e.g., s_copy, s_cmp)."
echo " It is a total mess when an actual string is mixed in C and FORTRAN.                  "
echo " s_cat is an exception and it is modified from the original definition in f2clib.     "
echo " " 
echo "                       Do not link against -lf2c.                                     "
