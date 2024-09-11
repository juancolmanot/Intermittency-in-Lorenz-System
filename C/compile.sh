#!/bin/bash

help_panel(){
    echo -e "COMPILE"
    echo -e "\tcompile: Es una herramienta que permite compilar código en C."
    echo -e "SINOPSIS"
    echo -e "\tcompile [programa a compilar] [argumentos]"
    echo -e "OPCIONES"
    echo -e "\tPrograma a compilar"
    echo -e "\t\tNombre del programa a ejecutar (debe estar previsto en el Makefile)"
}

compile_script(){
    make clean > /dev/null
    echo "[+] Make command info"
    n_headers=$(grep -c '.h"' src/$1.c)
    if [ n_headers == 0 ]; then
        make MAIN=$1
    else
        local_headers=""
        common_headers=""
        for header in $(grep '.h"' src/$1.c | awk -F'/' '{print $NF}' | awk -F'"' '{print $1}'); do
            if [ -f "../../common/include/$header" ]; then
                common_headers+="${header%.h} "
            else
                local_headers+="${header%.h} "
            fi
        done
        main_source=$(basename src/$1.c .c)
        make MAIN=$1 LOCAL_HEADERS="$local_headers" COMMON_HEADERS="$common_headers"
    fi
}

if [ -z $1 ]; then
    help_panel
    exit 1
elif [ $1 == "-h" ]; then
    help_panel
    exit 1
else
    compile_script $1 "$@"
fi
exit 0
