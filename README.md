# Analysis Software for the High Resolution Tests at CERN & DESY

## Requirements
 - python>=3.6
 - (py)ROOT>=6.22
 - [proteus](https://github.com/diamondIPP/proteus)
 - [judith](https://github.com/diamondIPP/judith) (only for CERN data)
 - [eudaq2](https://github.com/diamondIPP/eudaq-2) (only for DESY data)

## Installation
 - downloading code
   ```shell 
   git clone --recurse-submodules https://github.com/diamondIPP/HighResAnalysis.git 
   cd HighResAnalysis
   ```
 - python3, pip, virtualenv (requires sudo)
   ```shell 
    make prepare-dev 
   ```
 - virtual environment and python packages (included in previous step)
    ```shell
    make venv
   ```
 - required python packages are listed in [requirements](requirements.txt)
 - activate virtual environment and install aliases
    ```shell
   source .bash_aliases
   ```

## Running
 ```shell
analyse <run_number> <dut_number=0>
```
 - for more information run  ``` analyse -h ```

## Further reading
More information how the code is structured may be found [here](src/README.md).