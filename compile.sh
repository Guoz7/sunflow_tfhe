#!/bin/bash
gdbdebug=1 
flag="-L/usr/local/lib"
if [ "$gdbdebug" = "1" ]; then 
 flag+=" -g"
fi

if [ "$1" = "static" ]; then
    flag+=" -static"
fi
gcc $flag  tfhe_sunflow.c  
./a.out

####################
