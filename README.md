# Random DNA Chemistry (Python version)

Random DNA strand displacement circuits with Gillespie algorithm implementation, written in Python 3.7.3

## Steps to obtain Random DNA Strand Displacement RC simulation results

* Change the params in params_dsd.py (change the hold time tau value and the experiment index) For example, if you change the value for t_hold to be 0.1 and the exp to be "exp1", then your result files will be in redesign/visualDSD/tau
* dsd_test.py: It will prompt you on what to do (load the file to visual DSD, then download the result and put in the directory, etc.) I had to do it manually but if we want to do it automatically we have to find a way to interface with visual DSD
* dsd_combine.py: This will combine all the results from Visual DSD into a big CSV file
* compare_long_term.py: this will consume the big CSV file for that experiment and output result for short term memory task
* compare_short_term.py: this will consume the big CSV file for that experiment and output result for short term memory task
* Go back to step 1 and change to the next experiment and repeat step 1 to 5. Repeat until a desire number of experiment (I did 10 experiments)
* After completing the number of experiment for t_hold = 0.1, then change t_hold to 0.2 and repeat step 1 to 6.

Overview: 
* Obtain results for various tau (input hold time) values: 0.1, 0.2, 0.3, 0.4, 0.5, 0.5 seconds
* For each tau value: do 10 experiments ('exp1', 'exp2', 'exp3', ..., 'exp9', 'exp10')

Go to `redesign/` directory
For each tau value (Change `time_params['t_hold']` in `params_dsd.py`):
    For each experiment (Change `exp` string in `params_dsd.py`):
        Follow steps in `dsd_test.py` to generate simulation results for different periods on that experiment
        Run `dsd_combine.py` to combine all the CSV files into one big CSV file containing all the simulation results
        Run `compare_short_term.py` to compute the NRMSE for Short-term Memory Task
        Run `compare_long_term.py` to compute the NRMSE for Long-term Memory Task
    Run `show_task_results.py` to compute the average and variance NRMSE of both tasks for this tau value
    
