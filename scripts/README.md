## Analysis scripts

`PDbb_get_betatrace.py`: Extrate ROI time-series, band-pass filter, and do Hilbert transform.

`PDbb_makeMNE.py`: do MEG source reconstruction with dSPM.

`PDbb_srcMaps.py`: plot and inspect dSPM source reconstructions.

`burstanalysis_extended.m`: Run the analysis in `burstanalysis_stc.m` across all thresholds.

`burstanalysis_stc.m`: find bursts in beta trace and make summaries.

`burstsummary.m`: Arrange output of `burstanalysis_stc.m` for I/O.

`fooof_analysis.py`: 1/f regression analysis of PSD.

`import_filter_ica_bb.py`: Import raw data and do ICA to remove artefacts.

`poweranalysis.m`: calculate PSD.
