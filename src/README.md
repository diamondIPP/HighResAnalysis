# High Resolution Analysis Code

## Run logs
There are three beam tests with high resolution data:
 - [CERN 09/2018](https://docs.google.com/spreadsheets/d/1KoDi9OLU0SiqtLvTGgglQGHm51xn5J0ErAs89h6a-Rc/edit#gid=0)
 - [CERN 10/2018](https://docs.google.com/spreadsheets/d/1t-MXNW0eN9tkGZSakfPdmnd_wcq4cX14Nw0bQ2ma_OQ/edit#gid=0)
 - [DESY 12/2019](https://docs.google.com/spreadsheets/d/1vtwJnPLbk0M1UztpSX9SZNsPYAyMCO0TnYQzD6jQWoo/edit#gid=0)

### Conversion to json files
 - requires _client_secret.json_ file which is stored on mutter:/home/reichmann/software/HighResAnalysis/config 
 - convert online runlogs to a runlog.json file with [src/spreadsheet.py](spreadsheet.py):
    ```shell
    make_runlog <YYYYMM>
    ```
 - runlogs are stored in the data directory which may be set in the config

## Configuration
 - [config](../config) directory contains two files per default:
   - [default.ini](../config/default.ini)
   - [dia_info.json](../config/dia_info.json)
 - create main.ini based on default.ini and adjust to your preferences
   - [data]/server names the host where the raw data is stored
 - move _client_secret.json_ file to [config](../config) directory as described above

## Raw file conversion
 - handled by [converter.py](converter.py)
 - different raw files for CERN and DESY
### CERN
 - subclasses [converter.py](converter.py) in [cern/converter.py](../cern/converter.py)
 - requires three raw files:
   - RawDUT: ROOT file from the DUT: ljutel_<cms-run>.root
   - RawTEL: binary file from the telescope: acq<tel-run>.bin
   - RawREF: ROOT file from the reference plane (for alignment): anchor<tel-run>.root
 - step -2: convert the RawTEL file to ROOT (with judith)
 - step -1: convert _adc_ to _vcal_ in the RawDUT file (see [ADC Calibration](#adc-calibration)))
 - step  0: merge the two ROOT files for TEL and DUT
### DESY
 - only one raw file: run<nr>_<date>.raw
 - step 0: convert binary file to ROOT file
### Conversion
 - ROOT files after step 0 are identical for CERN & DESY
 - steps:
 1. noisescan             (proteus)
 2. telescope alignment   (proteus)
 3. track reconstruction  (proteus)
 4. root -> hdf5          (python)
 - steps 1 & 2 are skipped if the alignment data exists in the [proteus directory](../proteus)

## Event Alignment
 - DESY: EUDAQ data always aligned
 - CERN: 
   - Telescope sometimes merges two triggers -> event misalignment between RawDUT and RawTEL
   - [event_alignment.py](../cern/event_alignment.py) identifies the misaligned events
   - realignment applied in [adc.py](../cern/adc.py) (step -1)

## ADC Calibration
 - pxar calibration files stored in [calibration dir](../calibration)
 - named: phCal0.dat, phCal185.dat, ...
 - phCal0.dat is then valid for runs 0-184 