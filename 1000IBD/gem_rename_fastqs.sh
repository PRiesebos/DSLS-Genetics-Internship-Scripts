#!/usr/bin/env bash
###################################################################
#Script Name	  : gem_rename_fastqs.sh
#Description	  : 
#Args           :
#Author       	: Roy Oelen
#example        : ./gem_rename_fastqs.sh \
# \

###################################################################

# the parameters supplied
INPUT_DIR=$1
SEARCH='fq'
REPLACE='fastq'
FILEREGEX='.+'${SEARCH}.+

# define function to check files
RECURSIVE_RENAME () {
  # get passed variable
  PASSED_DIR=$1
  # list all files that have the search
  dirlist=(${PASSED_DIR}*)
  # loop the files and directories
  for e in "${dirlist[@]}"
  do
  # if is is a directory, go deeper
  if [ -d "$e" ];
  then
    RECURSIVE_RENAME ${e}'/'
  fi
  # if it is a file try a rename
  if [ -f "$e" ];
  then
    # now check if it matches with our search
    if [[ ${e##*/} =~ ${FILEREGEX} ]];
    then
      # get just the basename
      base_e="$(basename -- $e)"
      # do a search and replace
      e_replaced=$(echo $base_e | sed "s/$SEARCH/$REPLACE/g")
      # now actually move the file
      mv ${PASSED_DIR}'/'${base_e} ${PASSED_DIR}'/'${e_replaced} 
    fi
  fi
  done
}
# call function
RECURSIVE_RENAME ${INPUT_DIR}

