#==================
# Script: df_compile.R
# Author: Aylin Dincer
# Date updated: 11/19/2020

# Purpose: Clean and organize all imaging data and reformat for distribution.

#Instruction:
# 1. First run the df_curl_xml.sh script before moving forward.
#
# 2. Make sure you have the following directory structure and files in your working dir:
#     WORKING_DIR/data
#     WORKING_DIR/demographics
#     WORKING_DIR/xml_pull
#     WORKING_DIR/required_docs/fs_columns.csv
#     WORKING_DIR/required_docs/mr_fieldstrength_missing.csv
#     WORKING_DIR//required_docsnon.csv
#
# 2. If already ran the above bash script, do the following:
#     - update `xml_date` variable to the targeted xml dates in filename
#     - update `vDF` variable to the current data freeze release with the following format: DF#
#     - update any change to the project list
#
# 3. Run R script
#
# 4. Troubleshooting / Catching missed data
#       - `dup_` variables  = 0 if there are no duplicate processing and are > 0 if duplicates are found.
#       - `misproc_` variables = 0 if all sessions that weren't processed have a MR/PET Status, and are > 0 if session wasn't processed and also has no status
#       - `reqstat_` variables 
#==================


# Setup Workspace
#==================
# Clear variables
rm(list = ls())

# Load in libraries
library(stringr)
library(dplyr)
#==================


# USER INPUT
#==================
# Date on XML
xml_date <- "20201119"

# ADRC Data Freeze version
vDF <- "DF17"

# CNDA projects that should be included:
projects <- c("project1",
              "project2",
              "project3",
              "project4",
              "project5",
              "project6",
              "project7",
              "project8",
              "project9",
              "project10",
              "project11",
              "project12",
              "project13",
              "project14",
              "project15",
              "project16")
#==================






#Do Not Edit Beyond This Point Unless Needed.
#=================================================





#Functions
#==================
demographicsfs <- function(mydata, field_strength) {
  passedfs <- c("Passed",
                "Passed with edits")
  
  mrfs_usable <- filter(mydata, FS_QC_Status %in% passedfs)
  
  
  dem_table <- data.frame(matrix(ncol = 4, nrow = 1))
  colnames(dem_table) <- c("n", "percent_male", "long_n", "total_session")
  
  
  dem_table$n[1] <- length(unique(mrfs_usable$MAP))
  dem_table$percent_male[1] <- round((sum(mrfs_usable$Sex == "M")/length(mrfs_usable$MRI_Accession))*100,1)
  nlong <- mrfs_usable %>%
    group_by(MAP) %>% 
    dplyr::filter( n() > 1 )
  dem_table$long_n[1] <- length(unique(nlong$MAP))
  dem_table$total_session[1] <- length(unique(mrfs_usable$MRI_Accession))
  write.csv(dem_table, file = paste0("demographics/",xml_date,"_",vDF, "_",field_strength,"FS_dem.csv"), na = "", row.names = FALSE)
  
}


demographicspup <- function(mydata, tracer, proctype) {
  passedpup <- c("Passed",
                "Passed with edits")
  
  petpup_usable <- filter(mydata, PUP_QC_Status %in% passedpup)
  
  
  dem_table <- data.frame(matrix(ncol = 4, nrow = 1))
  colnames(dem_table) <- c("n", "percent_male", "long_n", "total_session")
  
  
  dem_table$n[1] <- length(unique(petpup_usable$MAP))
  dem_table$percent_male[1] <- round((sum(petpup_usable$Sex == "M")/length(petpup_usable$PET_Accession))*100,1)
  nlong <- petpup_usable %>%
    group_by(MAP) %>% 
    dplyr::filter( n() > 1 )
  dem_table$long_n[1] <- length(unique(nlong$MAP))
  dem_table$total_session[1] <- length(unique(petpup_usable$PET_Accession))
  write.csv(dem_table, file = paste0("demographics/",xml_date,"_",vDF, "_",tracer,"_",proctype,"_dem.csv"), na = "", row.names = FALSE)
  
}
#==================













# Read CSVs
#===================================
# XML data
#==================
mr <- read.csv(paste0("xml_pull/",xml_date,"_mr.csv"), 
                 header = TRUE, stringsAsFactors = FALSE)
petmr <- read.csv(paste0("xml_pull/",xml_date,"_petmr.csv"), 
               header = TRUE, stringsAsFactors = FALSE)
fsmr <- read.csv(paste0("xml_pull/",xml_date,"_fsmr.csv"), 
                  header = TRUE, stringsAsFactors = FALSE)
fspet <- read.csv(paste0("xml_pull/",xml_date,"_fspet.csv"), 
                  header = TRUE, stringsAsFactors = FALSE)
pet <- read.csv(paste0("xml_pull/",xml_date,"_pet.csv"), 
                header = TRUE, stringsAsFactors = FALSE)
pup <- read.csv(paste0("xml_pull/",xml_date,"_pup.csv"), 
                header = TRUE, stringsAsFactors = FALSE)
manpup <- read.csv(paste0("xml_pull/",xml_date,"_manpup.csv"), 
                header = TRUE, stringsAsFactors = FALSE)
wmhmr <- read.csv(paste0("xml_pull/",xml_date,"_wmhmr.csv"), 
                   header = TRUE, stringsAsFactors = FALSE)
wmhpet <- read.csv(paste0("xml_pull/",xml_date,"_wmhpet.csv"), 
                  header = TRUE, stringsAsFactors = FALSE)
radreadmr <- read.csv(paste0("xml_pull/",xml_date,"_radreadmr.csv"), 
                  header = TRUE, stringsAsFactors = FALSE)
radreadpet <- read.csv(paste0("xml_pull/",xml_date,"_radreadpet.csv"), 
                   header = TRUE, stringsAsFactors = FALSE)
#==================

# Read other required csvs
#==================
# CSV to easily change FS variable to biostats variant
fs_names <- read.csv("required_docs/fs_columns.csv",
                     header = TRUE, stringsAsFactors = FALSE)

# CSV to fill in missing field strength values
fstrength_miss <- read.csv("required_docs/mr_fieldstrength_missing.csv", 
                           header = TRUE, stringsAsFactors = FALSE)

# CSV to remove non-K participants
non_k <- read.csv("required_docs/non_k.csv",
                     header = TRUE, stringsAsFactors = FALSE)
non_k <- non_k[c("MAP", "STUDY")]
#==================
#===================================
















#MR SECTION
#===================================

# Clean data frame of only needed columns

# If any additional variables need to be included, do it here.
#===================================
#MR
#==================
mr <- mr[,c("label",
            "session_id",
            "xnat_subjectdata_subject_label",
            "xnat_subjectdata_subjectid",
            "xnat_subjectdata_map",
            'age',
            'xnat_subjectdata_gender_text',
            "date",
            "project",
            "xnat_subjectdata_projects",
            "xnat_subjectdata_xnat_subjectdata_field_map_grant",
            "scanner_csv",
            "field_strength",
            "xnat_mrsessiondata_field_map_mrstatus",
            "xnat_mrsessiondata_field_map_mrnotes")]

colnames(mr) <- c("MR_Session", 
                  "MRI_Accession",
                  "Subject",
                  "Subject_Accession",
                  "MAP",
                  "Age",
                  "Sex",
                  "MR_Date",
                  "Project",
                  "Projects",
                  "Grant",
                  "Scanner",
                  "Field_Strength",
                  "MR_status",
                  "MR_notes")
#==================

# If any additional variables need to be included, do it here.
#===================================
#PETMR
#==================
petmr$Field_Strength <- "3"
petmr <- petmr[,c("label",
            "session_id",
            "xnat_subjectdata_subject_label",
            "xnat_subjectdata_subjectid",
            "xnat_subjectdata_map",
            'gen_age',
            'xnat_subjectdata_gender_text',
            "date",
            "project",
            "xnat_subjectdata_projects",
            "xnat_subjectdata_xnat_subjectdata_field_map_grant",
            "scanner",
            "Field_Strength",
            "xnat_petsessiondata_field_map_mrstatus",
            "xnat_petsessiondata_field_map_mrnotes")]

colnames(petmr) <- c("MR_Session", 
                  "MRI_Accession",
                  "Subject",
                  "Subject_Accession",
                  "MAP",
                  "Age",
                  "Sex",
                  "MR_Date",
                  "Project",
                  "Projects",
                  "Grant",
                  "Scanner",
                  "Field_Strength",
                  "MR_status",
                  "MR_notes")
#==================

# If any additional variables need to be included, do it here.
#===================================
#FSMR
#==================
fsmr <- fsmr[ , !(names(fsmr) %in% c("xnat_mrsessiondata_label",
                                     "xnat_mrsessiondata_insert_user",
                                     "xnat_mrsessiondata_insert_date",
                                     "xnat_mrsessiondata_projects",
                                     "xnat_mrsessiondata_project",
                                     "fs_qc_method",
                                     "fs_qc_date",
                                     "fs_qc_validatedby",
                                     "quarantine_status"))]

colnames(fsmr) = fs_names$biostat_cols[match(colnames(fsmr), fs_names$xml_cols)]
#==================

# If any additional variables need to be included, do it here.
#===================================
#FSPET
#==================
fspet <- fspet[ , !(names(fspet) %in% c("xnat_petsessiondata_label",
                                     "xnat_petsessiondata_insert_user",
                                     "xnat_petsessiondata_insert_date",
                                     "xnat_petsessiondata_projects",
                                     "xnat_petsessiondata_project",
                                     "fs_qc_method",
                                     "fs_qc_date",
                                     "fs_qc_validatedby",
                                     "quarantine_status"))]

colnames(fspet) = fs_names$biostat_cols[match(colnames(fspet), fs_names$xml_cols)]
#==================
#===================================

# Remove and clean up MR data
#===================================
# Fill in missing field strength
#==================
mr <- left_join(mr, fstrength_miss, by = 'MR_Session')
mr <-  mutate(mr, Field_Strength = coalesce(Field_Strength.x,Field_Strength.y))
mr$Field_Strength <- round(mr$Field_Strength, 1)
mr <- mr[ , !(names(mr) %in% c("Field_Strength.x","Field_Strength.y"))]
#==================

# Keep only PETMR data
#==================
petmr_scanner <- c("CCIR mMR",
                   "CCIRWP-PETMR",
                   "CCIR PET/MR",
                   "WP-PETMRA",
                   "WP-PETMRA",
                   "SIEMENS Biograph_mMR CCIRWP-PETMR")
petmr <- filter(petmr, Scanner %in% petmr_scanner)

#==================

# Combine MR and PETMR session
#==================
mr_all <- rbind(mr, petmr)
#==================

# Keep only the ADRC CNDA projects
# This is a precaution (since CNDA can do weird things)
#==================
mr_all <- filter(mr_all, Project %in% projects)
#==================

# Remove MR sessions (duplicates, split sessions, tests)
#==================
rm_mrstatus <- c("Duplicate Session",
                 "Split Session - MPRAGE in Another Session",
                 "Phantom or Test Session")

mr_all <- filter(mr_all, !MR_status %in% rm_mrstatus)
#==================
#===================================

# Remove non-K participants
#===================================
# Remove non-MAP
#==================
mr_all <- mr_all[!(is.na(mr_all$MAP) | mr_all$MAP==""), ]

#==================

# Identify non_k2
#==================
non_k2 <- non_k[!(is.na(non_k$MAP) | non_k$MAP==""), ]
non_k2$MAP <- as.character(non_k2$MAP)

mr_all <- left_join(mr_all, non_k2, by = "MAP")
#==================
#===================================

# Remove and clean up Fs data
#===================================
# Seperate 1.5T and 3T MR Data
#==================
mr_15T <- filter(mr_all, Field_Strength == '1.5')
mr_3T <- filter(mr_all, Field_Strength == '3')
#==================

# Combine fspet and fsmr
#==================
fs <- rbind(fspet, fsmr)
#==================

# Remove unusable FS values
#==================
rm_fs <- c("Failed",
           "",
           "Quarantined",
           "Failed-needs reprocessing")
fs[fs$FS_QC_Status %in% rm_fs, (grepl('MR_', colnames(fs)) |  grepl('surfarea', colnames(fs)))] <- NA
#==================

#Keep only FS with a final QC status
#==================
keep_fs <- c("Passed",
                 "Failed",
                 "Passed with edits")

fs <- filter(fs, FS_QC_Status %in% keep_fs)
#==================

# Keep only FS5.3 and 5.0
#==================
fs_ver <- c("freesurfer-Linux-centos6_x86_64-stable-pub-v5.3.0-HCP-patch",
            "freesurfer-Linux-centos5_x86_64-stable-pub-v5.0.0")
fs <- filter(fs, FS_Version %in% fs_ver)
#==================

# Check if multiple FS within a session
# fs_check should equal 0
#==================
dup_fs <- fs %>% 
  group_by(MRI_Accession) %>% 
  dplyr::filter( n() > 1 )
#==================
#===================================

#Combine MR data with FS data
#===================================
# Merge data
#==================
mrfs15T <- left_join(mr_15T,fs,by = "MRI_Accession")
mrfs3T <- left_join(mr_3T,fs,by = "MRI_Accession")
#==================

# Reorder by FS QC Status
#==================
mrfs15T <- arrange(mrfs15T, FS_QC_Status)
mrfs3T <- arrange(mrfs3T, FS_QC_Status)
#==================

# Write FS output
#==================
write.csv(mrfs15T, file = paste0("data/",xml_date, "_",vDF, "_15TFS.csv"), na = "", row.names = FALSE)
write.csv(mrfs3T, file = paste0("data/",xml_date, "_",vDF, "_3TFS.csv"), na = "", row.names = FALSE)
#==================
#===================================

# MR Demographics
#===================================
demographicsfs(mrfs15T, '15T')
demographicsfs(mrfs3T, '3T')
#===================================















#PET SECTION
#===================================

# Clean data frame of only needed columns
#===================================

# If any additional variables need to be included, do it here.
#===================================
#PET
#==================
pet <- pet[,c("label",
            "session_id",
            "xnat_subjectdata_subject_label",
            "xnat_subjectdata_subjectid",
            "xnat_subjectdata_map",
            'gen_age',
            'xnat_subjectdata_gender_text',
            "date",
            "project",
            "xnat_subjectdata_projects",
            "tracer_name",
            "xnat_subjectdata_xnat_subjectdata_field_map_grant",
            "scanner_name",
            "xnat_petsessiondata_field_map_petstatus",
            "xnat_petsessiondata_field_map_petnotes",
            "xnat_petsessiondata_field_map_tbichronicblows",
            "xnat_petsessiondata_field_map_tbisingleincident")]

colnames(pet) <- c("PET_Session", 
                  "PET_Accession",
                  "Subject",
                  "Subject_Accession",
                  "MAP",
                  "Age",
                  "Sex",
                  "PET_Date",
                  "Project",
                  "Projects",
                  "Tracer",
                  "Grant",
                  "Scanner",
                  "PET_status",
                  "PET_notes",
                  "TBIChronicBlows",
                  "TBISingleIncident")
#==================

# If any additional variables need to be included, do it here.
#===================================
#PUP
#==================
pup <- pup[ , !(names(pup) %in% c("xnat_petsessiondata_label",
                                  "xnat_petsessiondata_insert_user",
                                  "xnat_petsessiondata_insert_date",
                                  "xnat_petsessiondata_projects",
                                  "xnat_petsessiondata_project",
                                  "label",
                                  "insert_user",
                                  "insert_date",
                                  "projects",
                                  "project",
                                  "mrid",
                                  "pet_fspet_fs_qcstat",
                                  "pet_fspet_fs_qcnote",
                                  "pet_tc_qc_date",
                                  "pet_tc_qc_method",
                                  "pet_tc_qc_validatedby",
                                  "quarantine_status"))]

colnames(pup)[colnames(pup)=="xnat_petsessiondata_session_id"] <- "PET_Accession"
colnames(pup)[colnames(pup)=="expt_id"] <- "PUP_ID"
colnames(pup)[colnames(pup)=="date"] <- "PUP_Date"
colnames(pup)[colnames(pup)=="proctype"] <- "PUP_Type"
colnames(pup)[colnames(pup)=="model"] <- "PUP_Model"
colnames(pup)[colnames(pup)=="pup_puptimecoursedata_field_map_pupstandardprotocol"] <- "PUP_Standard_Protocol"
colnames(pup)[colnames(pup)=="pup_puptimecoursedata_field_map_pupprotocolreason"] <- "PUP_Protocol_Reason"
colnames(pup)[colnames(pup)=="pet_fspet_mr_label"] <- "Processed_with_MR_Session"
colnames(pup)[colnames(pup)=="pet_fspet_mr_date"] <- "Processed_with_MR_Date"
colnames(pup)[colnames(pup)=="fsid"] <- "Processed_with_FS"
colnames(pup)[colnames(pup)=="pet_tc_qc_status"] <- "PUP_QC_Status"
colnames(pup)[colnames(pup)=="pet_tc_qc_note"] <- "PUP_QC_Notes"


manpup <- manpup[ , !(names(manpup) %in% c("xnat_petsessiondata_label",
                                  "xnat_petsessiondata_insert_user",
                                  "xnat_petsessiondata_insert_date",
                                  "xnat_petsessiondata_projects",
                                  "xnat_petsessiondata_project",
                                  "label",
                                  "insert_user",
                                  "insert_date",
                                  "projects",
                                  "project",
                                  "mrid",
                                  "pet_fspet_fs_qcstat",
                                  "pet_fspet_fs_qcnote",
                                  "pet_tc_qc_date",
                                  "pet_tc_qc_method",
                                  "pet_tc_qc_validatedby",
                                  "quarantine_status"))]

colnames(manpup)[colnames(manpup)=="xnat_petsessiondata_session_id"] <- "PET_Accession"
colnames(manpup)[colnames(manpup)=="expt_id"] <- "PUP_ID"
colnames(manpup)[colnames(manpup)=="date"] <- "PUP_Date"
colnames(manpup)[colnames(manpup)=="proctype"] <- "PUP_Type"
colnames(manpup)[colnames(manpup)=="model"] <- "PUP_Model"
colnames(manpup)[colnames(manpup)=="pup_puptimecoursedata_field_map_pupstandardprotocol"] <- "PUP_Standard_Protocol"
colnames(manpup)[colnames(manpup)=="pup_puptimecoursedata_field_map_pupprotocolreason"] <- "PUP_Protocol_Reason"
colnames(manpup)[colnames(manpup)=="pet_fspet_mr_label"] <- "Processed_with_MR_Session"
colnames(manpup)[colnames(manpup)=="pet_fspet_mr_date"] <- "Processed_with_MR_Date"
colnames(manpup)[colnames(manpup)=="fsid"] <- "Processed_with_FS"
colnames(manpup)[colnames(manpup)=="pet_tc_qc_status"] <- "PUP_QC_Status"
colnames(manpup)[colnames(manpup)=="pet_tc_qc_note"] <- "PUP_QC_Notes"

#==================
#===================================


# Remove and clean up MR data
#===================================
# Keep only the necessary projects
# This is a precaution
#==================
pet <- filter(pet, Project %in% projects)
#==================

# Remove MR sessions (duplicates, split sessions, tests)
#==================
rm_petstatus <- c("Duplicate Session",
                  "Split Session",
                  "Phantom or Test Session")

pet <- filter(pet, !PET_status %in% rm_petstatus)
#==================

# Remove non-K participants
#==================
# Remove non-MAP
#==================
pet <- pet[!(is.na(pet$MAP) | pet$MAP==""), ]
pet$MAP <- as.character(pet$MAP)
#==================

# Identify non_k
#==================
pet <- left_join(pet, non_k, by = "MAP")
#==================

#Keep PiB, AV45, AV1451, FDG Tracers
#==================
pet <-  filter(pet, Tracer %in% c("AV45","FDG","PIB","AV-1451"))
#==================
#===================================

# Remove and clean up PUP data
#===================================

#Remove PUPS run as a nonstandard protocol
#==================
pup <- filter(pup, !PUP_Standard_Protocol == 'False')
manpup <- filter(manpup, !PUP_Standard_Protocol == 'False')
#==================

# Remove unusable PUP values
#==================
rm_pup <- c("Failed",
           "",
           "Quarantined",
           "Failed-needs reprocessing")
pup[pup$PUP_QC_Status %in% rm_pup, grepl('pet_f', colnames(pup))] <- NA
manpup[manpup$PUP_QC_Status %in% rm_pup, grepl('pet_m', colnames(manpup))] <- NA
#==================

#Keep only PUPs with a final QC status
#==================
keep_pup <- c("Passed",
             "Failed",
             "Quarantined",
             "Passed with edits")

pup <- filter(pup, PUP_QC_Status %in% keep_pup)
manpup <- filter(manpup, PUP_QC_Status %in% keep_pup)
#==================

#Separate Manual and FreeSurfer proctypes
#==================
pup_man <-filter(manpup, PUP_Type == 'Manual')
pup_fs <- filter(pup, PUP_Type == 'FreeSurfer')
#==================

# Check if multiple PUPs within a session
# pup_check should equal 0
#==================
dup_pup_man<- pup_man %>% 
  group_by(PET_Accession) %>% 
  dplyr::filter( n() > 1 )


dup_pup_fs<- pup_fs %>% 
  group_by(PET_Accession) %>% 
  dplyr::filter( n() > 1 )
#==================
#===================================

#Separate into tracer
#===================================
# Merge PET with PUP data
#==================
petpup_fs <- left_join(pet,pup_fs,by = "PET_Accession")
petpup_man <- inner_join(pet,pup_man,by = "PET_Accession")
#==================

# Filter by tracer for FS proctype
#==================
PIB_fs <- filter(petpup_fs, Tracer == 'PIB')
AV45_fs <- filter(petpup_fs, Tracer == 'AV45')
AV1451_fs <- filter(petpup_fs, Tracer == 'AV-1451')
FDG_fs <- filter(petpup_fs, Tracer == 'FDG')
#==================

# Filter by tracer for manual proctype
#==================
PIB_man <- filter(petpup_man, Tracer == 'PIB')
AV45_man <- filter(petpup_man, Tracer == 'AV45')
AV1451_man <- filter(petpup_man, Tracer == 'AV-1451')
FDG_man <- filter(petpup_man, Tracer == 'FDG')
#==================
#===================================

# Add/Remove summary PET measures
#===================================
#Braak Regions and Tauopathy measure
#==================
AV1451_fs <- AV1451_fs %>% dplyr::mutate(Braak1_2 = AV1451_fs$pet_fsuvr_rsf_tot_ctx_entorhinal,
                                         Braak3_4 = rowMeans(AV1451_fs[,c("pet_fsuvr_rsf_tot_amygdala", 
                                                                          "pet_fsuvr_rsf_tot_accumbens",
                                                                          "pet_fsuvr_rsf_tot_hippocampus",
                                                                          "pet_fsuvr_rsf_tot_ctx_insula",
                                                                          "pet_fsuvr_rsf_tot_ctx_medorbfrn",
                                                                          "pet_fsuvr_rsf_tot_ctx_latorbfrn",
                                                                          "pet_fsuvr_rsf_tot_ctx_parsorbls",
                                                                          "pet_fsuvr_rsf_tot_ctx_parahpcmpl",
                                                                          "pet_fsuvr_rsf_tot_ctx_rosantcng",
                                                                          "pet_fsuvr_rsf_tot_ctx_postcng",
                                                                          "pet_fsuvr_rsf_tot_ctx_caudantcng",
                                                                          "pet_fsuvr_rsf_tot_ctx_isthmuscng")],na.rm = TRUE),
                                         Braak5_6 = rowMeans(AV1451_fs[,c("pet_fsuvr_rsf_tot_ctx_sstsbank",
                                                                          "pet_fsuvr_rsf_tot_ctx_midtmp",
                                                                          "pet_fsuvr_rsf_tot_ctx_caudmidfrn",
                                                                          "pet_fsuvr_rsf_tot_ctx_fusiform",
                                                                          "pet_fsuvr_rsf_tot_ctx_inferprtl",
                                                                          "pet_fsuvr_rsf_tot_ctx_infertmp",
                                                                          "pet_fsuvr_rsf_tot_ctx_latocc",
                                                                          "pet_fsuvr_rsf_tot_ctx_supramrgnl",
                                                                          "pet_fsuvr_rsf_tot_ctx_precuneus")], na.rm = TRUE),
                                         Tauopathy = rowMeans(AV1451_fs[,c("pet_fsuvr_rsf_tot_ctx_entorhinal",
                                                                           "pet_fsuvr_rsf_tot_amygdala",
                                                                           "pet_fsuvr_rsf_tot_ctx_infertmp",
                                                                           "pet_fsuvr_rsf_tot_ctx_latocc")], na.rm = TRUE))
#==================

#Calculate Centiloid
#==================
PIB_fs <- PIB_fs %>% dplyr::mutate(Centiloid_PiB_fBP_TOT_CORTMEAN = (126.7*PIB_fs$pet_fbp_tot_cortmean)-4.4,
                                   Centiloid_PiB_fSUVR_TOT_CORTMEAN = (111.8*PIB_fs$pet_fsuvr_tot_cortmean)-119.3,
                                   Centiloid_PiB_fBP_rsf_TOT_CORTMEAN = (53.9*PIB_fs$pet_fbp_rsf_tot_cortmean)-4.7,
                                   Centiloid_PiB_fSUVR_rsf_TOT_CORTMEAN = (45.0*PIB_fs$pet_fsuvr_rsf_tot_cortmean)-47.5)

PIB_man <- PIB_man %>% dplyr::mutate(Centiloid_PiB_mBP_TOT_CORTMEAN = (126.7*PIB_man$pet_mbp_tot_cortmean)-4.4,
                                   Centiloid_PiB_mSUVR_TOT_CORTMEAN = (111.8*PIB_man$pet_msuvr_tot_cortmean)-119.3)

AV45_fs <- AV45_fs %>% dplyr::mutate(Centiloid_AV45_fSUVR_TOT_CORTMEAN = (163.6*AV45_fs$pet_fsuvr_tot_cortmean)-181.0,
                                   Centiloid_AV45_fSUVR_rsf_TOT_CORTMEAN = (53.6*AV45_fs$pet_fsuvr_rsf_tot_cortmean)-43.2)

AV45_man <- AV45_man %>% dplyr::mutate(Centiloid_AV45_mSUVR_TOT_CORTMEAN = (163.6*AV45_man$pet_msuvr_tot_cortmean)-181.0)
#==================

#Remove Amyloid Summary measure from other non-amyloid tracers
#==================
AV1451_fs <- select(AV1451_fs, -grep('cortmean', colnames(AV1451_fs)))
FDG_fs <- select(FDG_fs, -grep('cortmean', colnames(FDG_fs)))
AV1451_man <- select(AV1451_man, -grep('cortmean', colnames(AV1451_man)))
FDG_man <- select(FDG_man, -grep('cortmean', colnames(FDG_man)))
#==================


# Additional Formating
#===================================
# Reorder by PUP QC Status
#==================
PIB_fs <- arrange(PIB_fs, PUP_QC_Status)
AV45_fs <- arrange(AV45_fs, PUP_QC_Status)
AV1451_fs <- arrange(AV1451_fs, PUP_QC_Status)
FDG_fs <- arrange(FDG_fs, PUP_QC_Status)
#==================

# Rename prefix on pup varaibles for BIOSTATS
#==================
PIB_fs <- rename_with(PIB_fs, ~ gsub("pet_f", "pib_f", .x, fixed = TRUE))
AV45_fs <- rename_with(AV45_fs, ~ gsub("pet_f", "av45_f", .x, fixed = TRUE))
AV1451_fs <- rename_with(AV1451_fs, ~ gsub("pet_f", "tau_f", .x, fixed = TRUE))
FDG_fs <- rename_with(FDG_fs, ~ gsub("pet_f", "fdg_f", .x, fixed = TRUE))

PIB_man <- rename_with(PIB_man, ~ gsub("pet_m", "pib_m", .x, fixed = TRUE))
AV45_man <- rename_with(AV45_man, ~ gsub("pet_m", "av45_m", .x, fixed = TRUE))
AV1451_man <- rename_with(AV1451_man, ~ gsub("pet_m", "tau_m", .x, fixed = TRUE))
FDG_man <- rename_with(FDG_man, ~ gsub("pet_m", "fdg_m", .x, fixed = TRUE))
#==================
#===================================

# PET Demographics
#===================================
demographicspup(PIB_fs, 'PIB', 'fs')
demographicspup(AV45_fs, 'AV45', 'fs')
demographicspup(AV1451_fs, 'AV1451', 'fs')
demographicspup(FDG_fs, 'FDG', 'fs')

demographicspup(PIB_man, 'PIB', 'man')
demographicspup(AV45_man, 'AV45', 'man')
demographicspup(AV1451_man, 'AV1451', 'man')
demographicspup(FDG_man, 'FDG', 'man')
#===================================

# Write PUP output
#==================
write.csv(PIB_fs, file = paste0("data/",xml_date, "_",vDF,"_PIB.csv"), na = "", row.names = FALSE)
write.csv(AV45_fs, file = paste0("data/",xml_date,"_",vDF, "_AV45.csv"), na = "", row.names = FALSE)
write.csv(AV1451_fs, file = paste0("data/",xml_date, "_",vDF,"_AV1451.csv"), na = "", row.names = FALSE)
write.csv(FDG_fs, file = paste0("data/",xml_date,"_",vDF, "_FDG.csv"), na = "", row.names = FALSE)
write.csv(PIB_man, file = paste0("data/",xml_date, "_",vDF,"_PIB_man.csv"), na = "", row.names = FALSE)
write.csv(AV45_man, file = paste0("data/",xml_date, "_",vDF,"_AV45_man.csv"), na = "", row.names = FALSE)
write.csv(AV1451_man, file = paste0("data/",xml_date, "_",vDF,"_AV1451_man.csv"), na = "", row.names = FALSE)
write.csv(FDG_man, file = paste0("data/",xml_date, "_",vDF,"_FDG_man.csv"), na = "", row.names = FALSE)
#==================
#===================================













# WMH SECTION
#===================================

# Clean data frame of only needed columns
#==================
# wmh mr
wmhmr <- wmhmr[ , !(names(wmhmr) %in% c("xnat_mrsessiondata_label",
                                        "xnat_mrsessiondata_insert_user",
                                        "xnat_mrsessiondata_insert_date",
                                        "xnat_mrsessiondata_projects",
                                        "xnat_mrsessiondata_project",
                                        "quarantine_status"))]

colnames(wmhmr) <- c("MRI_Accession",
                     "WMH_ID",
                     "WMH_volume")

# wmh pet
wmhpet <- wmhpet[ , !(names(wmhpet) %in% c("xnat_petsessiondata_label",
                                  "xnat_petsessiondata_insert_user",
                                  "xnat_petsessiondata_insert_date",
                                  "xnat_petsessiondata_projects",
                                  "xnat_petsessiondata_project",
                                  "label",
                                  "insert_user",
                                  "insert_date",
                                  "projects",
                                  "project",
                                  "quarantine_status"))]

colnames(wmhpet) <- c("MRI_Accession",
                      "WMH_ID",
                      "WMH_volume")

#==================

# Combine WMH processing from pet and mr sessions
#==================
wmh <- rbind(wmhmr,wmhpet)
#==================

# Check if multiple WMH within a session
# pup_check should equal 0
#==================
dup_wmh<- wmh %>% 
  group_by(MRI_Accession) %>% 
  dplyr::filter( n() > 1 )
#==================

# Merge WMH to session data
#==================
mrwmh <- inner_join(mr_all,wmh,by = "MRI_Accession")
mrwmh <- arrange(mrwmh, MR_Date)
#==================

# Write WMH output
#==================
write.csv(mrwmh, file = paste0("data/",xml_date, "_",vDF,"_WMH.csv"), na = "", row.names = FALSE)
#==================
#===================================














# RAD READS SECTION
#===================================
# Clean data frame of only needed columns
#==================
# radread mr
radreadmr <- radreadmr[ , !(names(radreadmr) %in% c("modality",
                                                    "reader",
                                                    "history",
                                        "quarantine_status"))]

colnames(radreadmr)[colnames(radreadmr) == "expt_id"] <- "RADREAD_ID"
colnames(radreadmr)[colnames(radreadmr) == "date"] <- "RADREAD_Date"
colnames(radreadmr)[colnames(radreadmr) == "session_id"] <- "MRI_Accession"

# radread pet
radreadpet <- radreadpet[ , !(names(radreadpet) %in% c("modality",
                                                    "reader",
                                                    "history",
                                                    "quarantine_status"))]
colnames(radreadpet)[colnames(radreadpet) == "expt_id"] <- "RADREAD_ID"
colnames(radreadpet)[colnames(radreadpet) == "date"] <- "RADREAD_Date"
colnames(radreadpet)[colnames(radreadpet) == "session_id"] <- "MRI_Accession"

#==================

# Combine radread processing from pet and mr sessions
#==================
radread <- rbind(radreadmr,radreadpet)
#==================

# Keep most recent rad read
#==================
radread <- radread %>% 
  group_by(MRI_Accession) %>% 
  arrange(RADREAD_Date, .by_group = TRUE) %>% 
  dplyr::filter(row_number()==n())
#==================

# check_radread should equal 0
#==================
dup_radread <- radread %>% 
  group_by(MRI_Accession) %>% 
  dplyr::filter( n() > 1 )
#==================

# Merge radread to session data
#==================
mrradread <- inner_join(mr_all,radread,by = "MRI_Accession")
mrradread <- arrange(mrradread, MR_Date)
#==================

# Write radread output
#==================
write.csv(mrradread, file = paste0("data/",xml_date, "_",vDF,"_RADREAD.csv"), na = "", row.names = FALSE)
#==================

#===================================







# MASTER PROJECT SECTION
#===================================
# Select only the necessary columns
#==================
cnda_mrfs15T <- mrfs15T[, c("Project", "Subject_Accession", "MRI_Accession","FS_ID","FS_QC_Status", "MAP", "MR_Session","STUDY")]
cnda_mrfs3T <- mrfs3T[, c("Project", "Subject_Accession", "MRI_Accession", "FS_ID", "FS_QC_Status", "MAP", "MR_Session","STUDY")]
cnda_petpup_fs <- petpup_fs[, c("Project", "Subject_Accession","PET_Accession", "PUP_ID", "PUP_QC_Status", "MAP", "PET_Session","STUDY")]
cnda_mrwmh <- mrwmh[, c("Project", "Subject_Accession","MRI_Accession", "WMH_ID", "MAP", "MR_Session","STUDY")]
cnda_mrradread <- mrradread[, c("Project", "Subject_Accession","MRI_Accession", "RADREAD_ID", "MAP", "MR_Session","STUDY")]
#==================

# Blank ASSESSOR ID if failed QC
#==================
cnda_mrfs15T$FS_ID <- ifelse(cnda_mrfs15T$FS_QC_Status %in% c("Failed", "Quarantined", NA), NA, cnda_mrfs15T$FS_ID )
cnda_mrfs3T$FS_ID <- ifelse(cnda_mrfs3T$FS_QC_Status %in% c("Failed", "Quarantined", NA), NA, cnda_mrfs3T$FS_ID )
cnda_petpup_fs$PUP_ID <- ifelse(cnda_petpup_fs$PUP_QC_Status %in% c("Failed", "Quarantined", NA), NA, cnda_petpup_fs$PUP_ID )
#==================

# Remove selected participants
#==================
remove_study <- c("non-k1", "non-k2")
cnda_mrfs15T <- cnda_mrfs15T[!cnda_mrfs15T$STUDY %in% remove_study ,]
cnda_mrfs3T <- cnda_mrfs3T[!cnda_mrfs3T$STUDY %in% remove_study ,]
cnda_petpup_fs <- cnda_petpup_fs[!cnda_petpup_fs$STUDY %in% remove_study ,]
cnda_mrwmh <- cnda_mrwmh[!cnda_mrwmh$STUDY %in% remove_study ,]
cnda_mrradread <- cnda_mrradread[!cnda_mrradread$STUDY %in% remove_study ,]
#==================

# Remove QC Status from data frames
#==================
cnda_mrfs15T <- cnda_mrfs15T[ ,  c("Project", "Subject_Accession","MRI_Accession","FS_ID","MAP","MR_Session")]
cnda_mrfs3T <- cnda_mrfs3T[ ,  c("Project", "Subject_Accession","MRI_Accession","FS_ID","MAP","MR_Session")]
cnda_petpup_fs <- cnda_petpup_fs[ ,  c("Project", "Subject_Accession","PET_Accession","PUP_ID","MAP","PET_Session")]
cnda_mrwmh <- cnda_mrwmh[ ,  c("Project", "Subject_Accession","MRI_Accession","WMH_ID","MAP","MR_Session")]
cnda_mrradread <- cnda_mrradread[ ,  c("Project", "Subject_Accession","MRI_Accession","RADREAD_ID","MAP","MR_Session")]
#==================

# Rename columns
#==================
master_cnda_cols <- c("PROJECT_ID", "SUBJECT_ID","EXPERIMENT_ID","ASSESSOR_ID", "SUBJECT_LABEL", "EXPERIMENT_LABEL")

colnames(cnda_mrfs15T) <- master_cnda_cols
colnames(cnda_mrfs3T) <- master_cnda_cols
colnames(cnda_petpup_fs) <- master_cnda_cols
colnames(cnda_mrwmh) <- master_cnda_cols
colnames(cnda_mrradread) <- master_cnda_cols
#==================

# Merge
#==================
cnda_final <- rbind(cnda_mrfs15T, cnda_mrfs3T, cnda_petpup_fs, cnda_mrwmh, cnda_mrradread)
cnda_final$NEW_PROJECT <- "DF"
write.csv(cnda_final, file = paste0("data/",xml_date, "_",vDF,"_mastercnda.csv"), na = "", row.names = FALSE)
#==================
#===============================================

# Troubleshooting steps
#===============================================

# Determine any session that weren't processed and 
# do not have a PET/MR status
#==================
missproc_mrfs15T <- filter(mrfs15T, (is.na(MR_status) | MR_status == "") & (is.na(FS_QC_Status) | FS_QC_Status == ""))
missproc_mrfs3T <- filter(mrfs3T, (is.na(MR_status) | MR_status == "") & (is.na(FS_QC_Status) | FS_QC_Status == ""))
missproc_petpup_fs <- filter(petpup_fs, (is.na(PET_status) | PET_status == "") & (is.na(PUP_QC_Status) | PUP_QC_Status == ""))
#==================

# Check if any PET sessions that initally required other info can now be processed.
#==================
reqstat_petpup_fs <- filter(petpup_fs, grepl('Require', PET_status))
#==================



# Clean up variables 
#===============================================
rm(AV1451_fs,AV1451_man,AV45_fs,AV45_man,cnda_final,cnda_mrfs15T,cnda_mrfs3T,cnda_mrradread,
   cnda_mrwmh,cnda_petpup_fs,FDG_fs,FDG_man,fs,fs_names,fsmr,fspet,
   fstrength_miss,manpup,mr,mr_15T,mr_3T,mr_all,mrfs15T,
   mrfs3T,mrradread,mrwmh,non_k,pet,petmr,petpup_fs,petpup_man,
   PIB_fs,PIB_man,pup, pup_fs, pup_man, radread, radreadmr, radreadpet,
   wmh, wmhmr, wmhpet, fs_ver,projects, keep_fs, keep_pup, master_cnda_cols,
   petmr_scanner, remove_study, rm_fs, rm_mrstatus, rm_petstatus, rm_pup)
#===============================================
