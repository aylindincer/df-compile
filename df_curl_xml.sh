#!/bin/sh

#==============================================================================
# Script: df_curl_xml.sh

# Author: Aylin Dincer
# Date updated: 09/15/2020

# Purpose: Pull MR / PET visit info and processed imaging data for the semi-annual
# data release.
#==============================================================================

#==============================================================================
# Usage:
# ./df_curl_xml.sh <alias> <secret> <input_dir> <output_dir>
#==============================================================================
#

#==============================================================================
# Required Inputs
# <alias> : alias token
# <secret> : secret tocken 
# <input_dir> : the directory where the xml files are located
# <output_dir> : the directory where you would like the pulled data to be saved.
#==============================================================================

#==============================================================================
# Ouput
# Each xml file performs a data search on the XNAT platform. The filtered imaging
# data is then saved as a comma separated value filetype. There should be a total
# of 11 csv files with the following nomenclature: [YYYYMMDD_xml_filename.csv]
#==============================================================================


# Define required arguments
ALIAS=$1
SECRET=$2
INPATH=$3
OUTPATH=$4

# Set date
now=$(date +"%Y%m%d")

# Define function to run curl call for each xml script
run_xml (){
	echo -e "----------------"
	echo -e "Running $1.xml: \n"
	curl -k -u ${ALIAS}:${SECRET} -X POST "https://xnat.org/d" --data-binary "@${INPATH}/$1.xml" > ${OUTPATH}/${now}_$1.csv

}


run_xml fsmr
run_xml fspet
run_xml manpup
run_xml mr
run_xml pet
run_xml petmr
run_xml pup
run_xml radreadmr
run_xml radreadpet
run_xml wmhmr
run_xml wmhpet
