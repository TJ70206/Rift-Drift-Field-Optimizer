# RDFO for CEC2017

This folder is used to run a minimal `RDFO` experiment on `CEC2017`.

## Main Files
- `runCEC2017.m`: main script. Set `functionId` manually before running.
- `RDFO.m`: RDFO algorithm.
- `Get_Functions_cec2017.m`: CEC2017 function interface.
- `cec17_func.mexw64`: compiled CEC2017 evaluator.
- `input_data17/`: CEC2017 data files.

## Usage
1. Open `runCEC2017.m`.
2. Set `functionId`.
3. Run the script in MATLAB.

## Results
Results are saved to:
- `results/Fx_D30/result.mat`
