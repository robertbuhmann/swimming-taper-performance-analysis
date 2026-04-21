# Swimming Taper Performance Analysis

Reproducible analysis code and supporting data for a study examining changes in heart rate variability, sleep, and mental fatigue across a pre-competition taper, and their associations with swimming performance.

## Repository contents

### Analysis scripts
- `primary_aim_analysis.R`  
  Primary athlete-level analysis examining associations between taper-period monitoring variables and performance change.

- `Taper_period_analysis.R`  
  Mixed-model analysis examining differences in monitoring variables between the pre-taper and taper periods.

### Data files

#### Raw data
- `taper_period_analysis.csv`  
  Raw dataset used for the taper-period analyses.

#### Processed data
- `swimming_training_load_imputed_2026-03-31.rds`  
  Processed dataset containing imputed mental fatigue data and derived variables used in the analyses.

- `m_primary_data_used.RDS`  
  Analysis-ready dataset used for the primary athlete-level model.

- `m_delta_data_used.RDS`  
  Analysis-ready dataset used for the delta model.

## Minimum reproducible workflow

### Primary aim
Run:

```r
source("primary_aim_analysis.R")
```

This script uses the processed data files to fit the athlete-level models examining associations between taper-period variables and competition performance.

### Secondary aim
Run:

```r
source("Taper_period_analysis.R")
```

This script uses `taper_period_analysis.csv` to examine differences between the pre-taper and taper periods for heart rate variability, sleep quality, sleep duration, and mental fatigue.


## File descriptions

### `taper_period_analysis.csv`
Contains the raw repeated-measures monitoring dataset used for the taper-period analyses.

### `swimming_training_load_imputed_2026-03-31.rds`
Contains the processed monitoring dataset after mental fatigue imputation and data preparation.

### `m_primary_data_used.RDS`
Contains the complete-case athlete-level data used in the primary performance model.

### `m_delta_data_used.RDS`
Contains the complete-case athlete-level data used in the delta model.

## Reproducibility note

This repository is intended to provide the minimum files required to reproduce the analyses reported in the manuscript. The included processed datasets are analysis-ready files derived from the broader study dataset and are provided to allow direct reproduction of the reported models.

## Citation

If using these materials, please cite the associated manuscript and Zenodo record.
