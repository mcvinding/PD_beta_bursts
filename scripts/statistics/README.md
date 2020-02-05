## Scripts for statistical inference

### Import scripts

`import2r.R`: import output of beta burst analysis (<script>) from Matlab to R and arrange data.

`import2r_extended.R`: Import output of extended analysis across threshold (`burstanalysis_extended.m`) from Matlab to R.

`import_UPDRS.R`: import MDS-UPDRS-III scores from database to R and arrange data with beta burst summary data (output from *import2r.R*).

### Analysis scripts

`fooof_stats.R`: 

`relpow_analysis.R`: import relative beta power from Matlab to R. Make summaries and do statistical comparison. Export data for plotting.

`roc_analyusis.R`: calculate area under the ROC curves.

`stats_extended.R`: Analysis of burst rate and calculate AU-ROC across thresholds.

`stats_main.R`: Analysis and summaries of burst rate, burst duration, intet-burst interval, and peak amplitude.

`UPDRS_stats.R`: Regression analysis of MDS-UPDRS-III scores and beta burst rate.