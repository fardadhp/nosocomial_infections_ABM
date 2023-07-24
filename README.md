Run from terminal:
$ monte_carlo_simulation.py [--nsamples] [--iterations] [--T0] [--burnIn] [--simLength] [--ifCalibration]

$ monte_carlo_simulation.py --help \
<pre>
Help
optional arguments: 
  -h, --help            show this help message and exit 
  -n, --nsamples        integer, number of set of samples (draws from distributions) (default: 100) 
  -i, --iterations      integer, number of simulations per set of samples (default: 30) 
  -t, --T0              string, starting date of simulation, format: 'YYYY-MM-DD' (default: None) 
  -b, --burnIn          integer, burn-in period (default: 30) 
  -l, --simLength       integer, simulation length (default: 360) 
  -c, --ifCalibration   True/False (default: True) 
</pre>
Further details: \
NSAMPLES: Set to a multyiple of your number of cores to maximize multiprocessing capabilities

ITERATIONS: A minimum of 30 is suggested for confidence interval calculations

T0: If not known, ignore or use None. The program will use the date of the first event in the event_queue file as T0.

IFCALIBRATION: If set to True, the sampling distributions will be read from *calibration/mc_parameters_ranges.csv*, otherwise, from *data/mc_parameters_ranges.csv*. This is in case a different set of prior distributions should be used for calibration vs. uncertainty analysis. The results of the simulation will be saved to *calibration* or *monte carlo* folders, respectively.
