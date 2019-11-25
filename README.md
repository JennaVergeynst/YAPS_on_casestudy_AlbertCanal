# Case study in the Albert Canal using YAPS

This repository provides examples for positioning acoustic telemetry detections with YAPS (Baktoft, Gjelland, Økland & Thygesen (2017): Positioning of aquatic animals based on time-of-arrival and random walk models using YAPS (Yet Another Positioning Solver). DOI:10.1038/s41598-017-14278-z), on a case study in the Albert Canal, Belgium. Different wrapping functions around YAPS and a script to run an example track are provided. The detection data and resulting YAPS tracks are available at http://doi.org/10.5281/zenodo.3544748 (Vergeynst, Jenna, Pauwels, Ine, Baeyens, Raf, Mouton, Ans, & Coeck, Johan. (2019). Acoustic positioning telemetry in a case study on the Albert Canal (Belgium): VPS and YAPS positions [Data set]. Zenodo.).

To run the code, installation of the YAPS package is required: https://github.com/baktoft/yaps.

## Preparatory steps

The raw detection data of this dataset were synchronised by use of following code: https://github.com/JennaVergeynst/time_synchronization.

To create the Time Of Arrival (TOA) matrix required by the YAPS algorithm, following code can be used: https://github.com/JennaVergeynst/prepare_toa_for_yaps.

Finally, the open-source Acoustic Positioning Telemetry website provides an overview of available python code and R packages: https://apostel-team.github.io/APosTel-team.github.io/.

The following sections give information on the parameters used to prepare the TOA data for each example. Explanation on the parameters can be found [here](https://github.com/JennaVergeynst/prepare_toa_for_yaps/blob/master/Prepare_toa_data.py).

## Fish positioning
Run the file fish_positions/run_yaps_on_fish_track.R to create YAPS tracks for fish with ID 100. To create the TOA matrix (TOA_data_100.csv), following parameters were used as input for the TOA preparation function:
- burst interval [17, 35] ([min_burst, max_burst])
- max_time=1
- pas_tol=30
- min_track_length=5

## Fixed tag positioning
Run the file fixed_positions/run_yaps_on_fixed_tag.R to create YAPS positions for synchronisation tag S15 (ID 65043), collocated with receiver ST15. Input parameters for the TOA preparation function:
- burst interval [540, 660]
- pas_tol=60
- min_track_length=10

## Test tracks
Run the file test_positions/run_yaps_on_test_track.R to create YAPS tracks for 3 test tracks:
- test 255: burst interval [17,35]
- test 53429: burst interval [15,30]
- test 16200: burst interval [10,18]
- pas_tol = 5
- max_time = 0.1
- min_track_length = 10