# Table of Contents
- [Data Release Overview](#data-release-overview)
- [XNAT Tokens](#xnat-tokens)
- [General instructions to run scripts](#general-instructions-to-run-scripts)
- [Extract data via curl call](#extract-data-via-curl-call)
- [Running the DF compile script](#running-the-df-compile-script)

<br>
<br>

# Data Release Overview
This repository contains scripts to compile the semi-annual data release.

# XNAT Tokens
XNAT tokens are a secure way of authenticating your request so that XNAT can determine whether you are allowed to perform the action that you are requesting. The tokens expire in a short period of time. Do not share your token with anyone else.

1. Log in to XNAT
2. Click on your username in the very top right navigation bar to go to your User Account page. 
3. Click on the "Manage Alias Tokens" tab. 
4. Click "Create Alias Token".
5. Click on the generated "Alias" link or the "view" button under the "Actions" column to view the full token
6. Copy the "alias" and "secret" strings and paste them into the following commands below when specified.

More detailed documentation can be found here: https://wiki.xnat.org/documentation/how-to-use-xnat/generating-an-alias-token-for-scripted-authentication


# General instructions to run scripts

1. Download repository

2. Make sure your working directory has the following files/folders:
```
data/
demographics/
required_docs/
xml_files/
xml_pull/
df_compile.R
curl_xml.sh

```

3. Check to make sure the `xml_files/` directory has 11 '.xml' files with the following names:
```
fsmr.xml
fspet.xml
manpup.xml
mr.xml
pet.xml
petmr.xml
pup.xml
radreadmr.xml
radreadpet.xml
wmhmr.xml
wmhpet.xml

```

4. On the terminal, change directories to where the above scripts and folders using ( `cd` ).


# Extract data via curl call
**df_curl_xml.sh**

This script pulls data specifically identified in xml files.

<br>

**General Usage:**
```
./curl_xml.sh <alias> <secret> <input_dir> <output_dir>
```

<br>

**Required inputs:**

`<alias>`: Obtain alias token using the instructions under "XNAT token"
  
`<secret>`: Obtain secret token using the instructions under "XNAT token"

`<input_dir>` : the directory where the xml files are located. Input directory must contain the following files: 


` <output_dir> `: the directory where you would like the pulled data to be saved. (xml_pull/)

 <br>
 
**Example Usage**

```
./curl_xml.sh ALIAS_TOKEN SECRET_TOKEN xml_files/ xml_pull/
```

where ALIAS_TOKEN is your specified alias token, SECRET_TOKEN is your specificed secret token

<br>

**Script output**

This script organizes the files into folders like this:

```
xml_pull/YYYYMMDD_*.csv
```
where YYYYMMDD is the year, month, and day the script was run and * is the data that is pulled from each xml.

For example, 20201119_mr.csv was pulled on November 19, 2020 and it's pulling data from the mr.xml file


# Running the DF compile script
**df_compile.R**

This script cleans and organizes all imaging data, reformat for distribution, and outputs the demographic tables needed for the data ditionary and creates a CSV used for adding data to the XNAT database.

1. Make sure you ran the df_curl_xml.sh script and 11 .csv files were created with the same date.

2. In the R script, update the following variables under the section "User Input"

    * `xml_date` variable set to the targeted xml file dates with the following format: YYYYMMDD

    * `vDF` variable to the current data freeze release with the following format: DF#

    * Update any changes to the project list under the `projects` variable

3.  Run the entire R script.

4. There will be a set of variables remaining that are used for troubleshooting. Check each session and fix if needed.
    * `dup_` variables  = 0 if there are no duplicate processing and are > 0 if duplicates are found.

    * `misproc_` variables = 0 if all sessions that weren't processed have a MR/PET Status, and are > 0 if session wasn't processed and also has no status.
  
    * `reqstat_` variable is a list of PET sessions that initally required other info before processing. A check will need to be done for each session to determine if that is still the case.

**Script output**

This script creates the following files:

```
data/YYYYMMDD_DF#_3TFS.csv
data/YYYYMMDD_DF#_15TFS.csv
data/YYYYMMDD_DF#_AV45.csv
data/YYYYMMDD_DF#_AV45_man.csv
data/YYYYMMDD_DF#_AV1451.csv
data/YYYYMMDD_DF#_AV1451_man.csv
data/YYYYMMDD_DF#_FDG.csv
data/YYYYMMDD_DF#_FDG_man.csv
data/YYYYMMDD_DF#_mastercnda.csv
data/YYYYMMDD_DF#_PIB.csv
data/YYYYMMDD_DF#_PIB_man.csv
data/YYYYMMDD_DF#_RADREAD.csv
data/YYYYMMDD_DF#_WMH.csv

demographics/YYYYMMDD_DF#_3TFS_dem.csv
demographics/YYYYMMDD_DF#_15TFS_dem.csv
demographics/YYYYMMDD_DF#_AV45_dem.csv
demographics/YYYYMMDD_DF#_AV45_man_dem.csv
demographics/YYYYMMDD_DF#_AV1451_dem.csv
demographics/YYYYMMDD_DF#_AV1451_man_dem.csv
demographics/YYYYMMDD_DF#_FDG_dem.csv
demographics/YYYYMMDD_DF#_FDG_man_dem.csv
demographics/YYYYMMDD_DF#_PIB_dem.csv
demographics/YYYYMMDD_DF#_PIB_man_dem.csv
demographics/YYYYMMDD_DF#_RADREAD_dem.csv
demographics/YYYYMMDD_DF#_WMH_dem.csv
```

`3TFS` : 3T MRI sessions with FreeSurfer data

`15TFS` : 1.5T MRI sessions with FreeSurfer data

`AV45` : AV45-PET sessions with FreeSurfer-derived PUP data

`AV45_man` : AV45-PET sessions with manual traced PUP data

`AV1451` : AV1451-PET sessions with FreeSurfer-derived PUP data

`AV1451_man` : AV1451-PET sessions with manual traced PUP data

`FDG` : FDG-PET sesisons with FreeSurfer-derived PUP data

`FDG_man` : FDG-PET sesisons with manual traced PUP data

`PIB` : PiB-PET sessions with FreeSurfer-derived PUP data

`PIB_man`: PiB-PET sessions with manual traced PUP data

`RADREAD`: Radiolgoical reads for the MR sessions

`WMH` : White matter hyperintensity volume

`dem` : demographc data for each modality