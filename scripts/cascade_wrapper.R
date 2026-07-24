#===================================== After clonning the repository==========================================================================================================================================

 

#      cascade_app/

#      └─── cascade_script_v2.R

#      └─── data/

#      |         └─── in/

#      |               └─── sample_a_pos.mzML

#      |                        └─── sample_a_neg.mzML

#      |                        └─── sample_a_features_pos.csv

#      |                        └─── sample_a_features_neg.csv

#      |                        └─── sample_a_annotations_pos.tsv

#      |                        └─── sample_a_annotations_neg.tsv

 

#==============================================================================================================================================================================================================

 

#===================================== REQUIREMENTS ===========================================================================================================================================================

 

#   A file (.mzML) containing DDA MS data with an additional detector (PDA, ELSD, CAD)

#   In case you don’t know how to obtain it, see: wiki/How-to-create-a-compliant-mzML-file

#   A file (.csv) containing features, as obtained by mzmine comprehensive export

#   A file (.tsv) containing annotations, as obtained by TIMA

# ===========================================================================================================================================================================================================

 

#===================================== INSTALLATION ==========================================================================================================================================================

 

install.packages(

  "cascade",

  repos = c(

    "https://adafede.r-universe.dev",

    "https://bioc.r-universe.dev",

    "https://cloud.r-project.org"

  )

)

# ===========================================================================================================================================================================================================

 

#===================================== CASCADE STARTS HERE===================================================================================================================================================

#load packagea

library(cascade)

library(dplyr)  

 

# Check documentation for the package

#  ?cascade

 

#check all the functions available

# help(package = "cascade") #shows all the functions available in the package

 

# =================================== PATHS ==================================================================================================================================================================

 

# Define the main folder of your project and file names

data_path <- "C:/Users/legregam/OneDrive - Université de Genève/Documents/Cascade/calendula140min"              #"path/to/your/data"  # Define the path to your data
filename_pos <- "20240628_LQ_Calendula_1_140min_Pos.mzML"          #your indexed mzml in pos mode 
filename_neg <- "20240628_LQ_Calendula_1_140min_Neg.mzML"          #your indexed mzml in neg mode 
feature_table_pos <- "Calendula_pos_quant_full.csv"                #your MZmine features table in POS mod               #your MZmine features table in pos mode 
feature_table_neg <- "calendula_neg_quant_full.csv"                #your MZmine features table in neg mode
annotations_results_pos <- "calendula_pos_tima_results.tsv"      #your TimaR annotatios results table in pos mode
annotations_results_neg <- "calendula_neg_tima_results.tsv"      #your TimaR annotatios results table in neg mode

 

# ===========================================================================================================================================================================================================

# Do not modify this section

data_in_path       <- file.path(data_path, "in/")

file_negative      <- file.path(data_in_path, filename_neg)

file_positive      <- file.path(data_in_path, filename_pos)

features_pos       <- file.path(data_in_path, feature_table_pos)

features_neg       <- file.path(data_in_path, feature_table_neg)

annotations_pos    <- file.path(data_in_path, annotations_results_pos)

annotations_neg    <- file.path(data_in_path, annotations_results_neg)

 

#create automatically a path to export data

export_path <- file.path(data_path, "results")  # Save results inside the 'cascade' folder

if (!dir.exists(export_path)) {

  dir.create(export_path)  # Create the directory if it doesn't exist

}

# ================================= Validate the file paths =================================================================================================================================================

if (!file.exists(file_negative)) {

  stop("Negative mzML file does not exist in the specified data_path.")

}

if (!file.exists(file_positive)) {

  stop("Positive mzML file does not exist in the specified data_path.")

}

if (!file.exists(features_pos)) {

  stop("Features table in pos file does not exist in the specified data_path.")

}

if (!file.exists(features_neg)) {

  stop("Features table in neg file does not exist in the specified data_path.")

}

if (!file.exists(annotations_pos)) {

  stop("Annotation results in pos file does not exist in the specified data_path.")

}

if (!file.exists(annotations_neg)) {

  stop("Annotation results in neg file does not exist in the specified data_path.")

}

# ===========================================================================================================================================================================================================

 

# ================================= VARIABLES ===============================================================================================================================================================

 

# 0.1. define the shifts between the CAD and PDA traces in reference to the MS trace

cad_shift <- 0.2

pda_shift <- 0.3

ms_shift  <- 0.0

 

# 0.1. define the rt window you want to consider

time_min <- 0.5

time_max <- 120

 

 

# ===========================================================================================================================================================================================================

 

# ================================= ALIGNMENT ===============================================================================================================================================================

 

# 1. Apply the check_chromatograms_alignment function to verify the data

 

result_plot <- check_chromatograms_alignment(

  file_negative = file_negative,

  file_positive = file_positive,

  time_min = time_min,             # Adjust as needed, chromatographic run time

  time_max = time_max,           # Adjust as needed, chromatographic run time

  cad_shift = cad_shift,           # CAD shift

  pda_shift = pda_shift,           # PDA shift

  fourier_components = 0.01,  # Fourier components

  frequency = 1,              # frequency

  resample = 1,               # Resampling factor

  chromatograms = c("bpi_pos", "cad_pos", "pda_pos", "bpi_neg"),  # Chromatograms to plot

  type = "baselined",         # Choose "baselined" or "improved"

  normalize_intensity = TRUE, # Normalize intensity?

  normalize_time = FALSE      # Normalize time?

)

 

# The function will return a plot

print("Alignment check completed. Resulting plot:")

print(result_plot)

 

#NOTE: modify the variable accordingly to achieve the best alignment according to your experiment data

# ===========================================================================================================================================================================================================

 

# ================================= INTEGRATIONS ============================================================================================================================================================

 

# 2. Apply the checK_peaks_integrations

df <- data.table::fread(features_pos)
names(df) |> grep("datafile", x = _, value = TRUE)

#2.1. cad trace

 

result_peaks_cad <- check_peaks_integration(

  file = file_positive,        # Positive mzML file

  features = features_pos,     # Features table

  detector = "cad",            # Specify detector (e.g., "cad, bpi or pda")

  chromatogram = "baselined",  # Use "baselined" or "improved" chromatograms

  min_area = 0.001,            # Minimum area threshold for integration

  min_intensity = 5E4,         # Minimum intensity threshold

  shift = cad_shift,                # Shift (CAD shift value)

  show_example = FALSE,        # Use your data, not example data

  fourier_components = 0.01,   # Fourier components for smoothing

  time_min = time_min,              # Minimum retention time

  time_max = time_max,            # Maximum retention time

  frequency = 1,               # Frequency

  resample = 1                 # Resampling factor

)

# Display the result

result_peaks_cad

 

#2.2. pda trace

 

result_peaks_pda <- check_peaks_integration(

  file = file_positive,        # Positive mzML file

  features = features_pos,     # Features table

  detector = "pda",            # Specify detector (e.g., "cad")

  chromatogram = "baselined",  # Use "baselined" or "improved" chromatograms

  min_area = 0.005,            # Minimum area threshold for integration

  min_intensity = 1e5,         # Minimum intensity threshold

  shift = pda_shift,           # Shift (CAD shift value)

  show_example = FALSE,        # Use your data, not example data

  fourier_components = 0.01,   # Fourier components for smoothing

  time_min = time_min,         # Minimum retention time

  time_max = time_max,         # Maximum retention time

  frequency = 1,               # Frequency

  resample = 1                 # Resampling factor

)

# Display the result

result_peaks_pda