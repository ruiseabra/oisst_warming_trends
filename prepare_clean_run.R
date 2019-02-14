# DELETE FOLDERS AND FILES FROM PREVIOUS RUNS OF THE CODE
# TO ALLOW FOR A NEW CLEAN RUN

# ATTENTION!! This script will delete all outputs from previous runs
# The code is wrapped in "if (FALSE)" statements to avoid accidental deletion of important data
# Run the code inside the if statements manually

# delete outputs from "oisst_download.R"
if (FALSE) {
  unlink("oisst_repository/", recursive = TRUE)
  unlink("oisst_download_outputs/", recursive = TRUE)
} 

# delete outputs from "oisst_warming_trends_analysis.R"
if (FALSE) {
  # "oisst_warming_trends_inputs/" should just have the gshhs coastline data, and thus should not be deleted
  unlink("oisst_warming_trends_outputs/", recursive = TRUE)
  unlink("final_stats/", recursive = TRUE)
  unlink("final_plots/", recursive = TRUE)
} 