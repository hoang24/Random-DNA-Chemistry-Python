import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import pandas as pd
from simulate_perturbed_chem import plot_concentration


directory = 'exp1'

df_dsdResult = pd.read_csv(f'visualDSD/{directory}/dsdResult.csv', engine='python')
dsdResult = df_dsdResult.to_dict('list')
dsdTime = np.array(dsdResult['Time '])
dsdConcentration = {key: dsdResult[key] for key in dsdResult if key != 'Time '}

df_pyResult = pd.read_csv(f'visualDSD/{directory}/pyResult.csv', engine='python')
pyResult = df_pyResult.to_dict('list')
pyTime = np.array(pyResult['Time'])
pyConcentration = {key: pyResult[key] for key in pyResult if key != 'Time'}

# plot_concentration(time_lookup=dsdTime, concentration_lookup=dsdConcentration, plot_name='visualDSD/dsdResult')
# plot_concentration(time_lookup=pyTime, concentration_lookup=pyConcentration, plot_name='visualDSD/pyResult')

def exp_decay_decrease(time, initial, rate_const):
    '''
        Exponential decay function (decreasing form): y = Ae^(-kx)
    '''
    species_count = initial * np.exp(-rate_const*time)
    return species_count

def exp_decay_increase(time, final, rate_const):
    '''
        Exponential decay function (increasing form): y = A - Ae^(-kx)
    '''
    species_count = final * (1 - np.exp(-rate_const*time))
    return species_count

for dsdKey, pyKey in zip(dsdConcentration.keys(), pyConcentration.keys()):
    # Convert the concentration data to numpy array
    dsdValue = np.array(dsdConcentration[dsdKey])
    pyValue = np.array(pyConcentration[pyKey])

    # Curve fitting the data from Visual DSD to return the parameters (initial/final and rate_const)
    paramsD, param_covD = curve_fit(f=exp_decay_decrease, xdata=dsdTime, ydata=dsdValue)
    paramsI, param_covI = curve_fit(f=exp_decay_increase, xdata=dsdTime, ydata=dsdValue)

    # Calculate the one standard deviation errors on the parameters use
    perrD = np.sqrt(np.diag(param_covD))
    perrI = np.sqrt(np.diag(param_covI))
    # if perrD[0] < perrI[0] and perrD[1] < perrI[1]:
    #     params = paramsD
    # elif perrI[0] < perrD[0] and perrI[1] < perrD[1]:
    #     params = paramsI
    # else:
    #     raise Exception('Neither exp_decay_decrease or exp_decay_increase is a good fit.')

    # print("Coefficients:")
    # print(params)
    # print('One standard deviation errors:')
    # print(f'decreasing: {perrD}; increasing: {perrI}')

    # Create a fit curve of the Visual DSD data on the DSD timeline and Python timeline
    fit_curveD = paramsD[0] * np.exp(-paramsD[1] * pyTime)
    fit_curveI = paramsI[0] * (1 - np.exp(-paramsI[1] * pyTime))

    # Generate plots
    plt.plot(dsdTime, dsdValue, ':', color='red', label=f'DSD Result: {dsdKey}')
    plt.plot(pyTime, pyValue, ':', color='blue', label=f'Python Result: {pyKey}')
    # plt.plot(pyTime, fit_curveD, '--', color='black', label=f'Fit Curve using Exponential Decay (Decreasing Form)')
    # plt.plot(pyTime, fit_curveI, '--', color='purple', label=f'Fit Curve using Exponential Decay (Increasing Form)')
    # plt.yscale('log')
    plt.title('Species Count of DSD and Python Results and Fit Curve')
    plt.xlabel('Time (s)')
    plt.ylabel('Species Count (species)')
    plt.legend()
    plt.show()

    # Comparison
    # diff_py_dsd_D = fit_curveD - pyValue
    # diff_py_dsd_I = fit_curveI - pyValue
    # plt.plot(pyTime, diff_py_dsd_D, '--', color='black', label='Error using Exponential Decay (Decreasing Form) Fit Curve')
    # plt.plot(pyTime, diff_py_dsd_I, '--', color='purple', label='Error using Exponential Decay (Increasing Form) Fit Curve')
    # plt.title('Difference of Python Result to DSD Result')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Species Count (species)')
    # plt.legend()
    # plt.show()
