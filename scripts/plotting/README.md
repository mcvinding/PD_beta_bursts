## Scripts for making plots and figures

`ROC_plots.R`: Plot ROC curves (Fig. 5 and Fig. S2).

`UPDRSmod_plots.R`: Plot MSD-UPDRS-III  scores against burst/min and prediction from GLM analysis with K-fold cross-validation (Fig. 7).

`eventplot.m`: plot ROI time-series data timelocked to peak of bursts (Fig. 4A).

`example_figures.m`: Plot part of the ROI time-series.

`plot_PSD.m` Plot PSD (Fig. 2A).

`plot_extended.R`: Plot n events across threshold in the extended analysis (`burstanalysis_extended.m`), output of the statistical analysis across thresholds (`statistics/stats_extended.R`) and area under the ROC curve across thresholds (Fig. 6).

`plot_fooof.R`: make dotplots of summary measures of the FOOOF (see `scripts/fooof_analysis.py`): beta peak power (Fig 2C), 1/f slope coefficients (Fig. 2D), and 1/f intervept coefficients(Fig. 2E).

`plot_relpow.R`: make dotplots of relative the beta power (Fig 2B).

`plots_main.R`: make dotplots of burst rate (Fig 3), pooled density plots of event duration (Fig. 4B), and pooled density plots of inter-event interval (Fig 4C).

`threshold_plot.m`: plot Pearson's correlation vs. threshold used to determine the threshold (Fig. 1B).
