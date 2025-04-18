#!/bin/bash

#command
usage() {
    echo "Usage: $0 --phecode <PHECODE> --pop <POP>"
    exit 1
}

#command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --phecode)
            PHECODE=$2
            shift 2
            ;;
        --pop)
            POP=$2
            shift 2
            ;;
        *)
            echo "unknown flag: $1"
            usage
            ;;
    esac
done

#check for required arguments
if [[ -z "$PHECODE" || -z "$POP" ]]; then
    usage
fi

#download hail table
python "/myrepo/pull_data.py" --phecode "$PHECODE" --pop "$POP"

#format hail tables
Rscript "/myrepo/table_format.R" --phecode "$PHECODE" --pop "$POP"

#GWAS qqman
#python "/myrepo/qqman.py" --phecode "$PHECODE" --pop "$POP"

#locuszoomR
#Rscript "/myrepo/locuszoom.R" --phecode "$PHECODE" --pop "$POP"