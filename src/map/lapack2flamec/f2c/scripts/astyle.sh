#!/bin/bash


main()
{
    # local variables
    script_name=${0##*/}
    src_dir="../c"

    echo " "
    echo " "
    echo " "${script_name}
    echo " "
    echo " Objective:"
    echo "   The script reformats  files."
    echo " "
    
    ast=$(which astyle)
    
    if [ -f "$ast" ]
    then
        astyle --style=allman ${src_dir}/*.c ${src_dir}/*.h    
        rm -f ${src_dir}/*.orig
    else
        echo "$ast is not found."
    fi
}

main "$@"
