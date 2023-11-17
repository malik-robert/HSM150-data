## Data from the experiments on the two-motor HSM150 system
- In the folder **Data**, there are two subfolders: **Control** and **Identification**.

### Identification
- This folder contains data from system identification.
- The data are labeled as `idf_sign_input_number.mat`
  - `sign` ... the input sign: negative-`neg` or positive-`pos`
  - `input` ... the magnitude of the input torque for the given experiment
  - `number` ... the number of the experiment

- e.g. `idf_neg_input_0.04.mat` contains the system response for the step of the input torque `-0.01 N m` at the operating point `-0.03 N m`. The resultant input torque is `-0.04 N m`.

### Control
- This folder contains data from the experiments of the system control.
- The data are labeled as `MIMO_Feedback_ReferenceSignal.mat`
  - `ReferenceSignal` ... the type of the reference signal, either `Constant` - step responses or `Harmonic` - responses to the harmonic reference signal

### Data Structure
- Each `.mat` contains a cluster of data, including the catalog values of the HSM150 system and the system parameters used in the experiments. The most important is the variable `data`, which contains the raw measurements.
- The variable `data` has the following structure:

#### Identification
```
| Time | Position Reference | Position | Velocity | Input Torque | Torque Measurement | Motor Electric Current |
| s    | rad                | rad      | rad/s    | Nm           | Nm                 | A                      |
```

#### Control
```
| Time | Position Reference | Position | Velocity Reference | Velocity | State Error | Velocity Error | Input Torque | Torque Measurement | Load Torque | Computed Torque Method | Motor Electric Current |
| s    | rad                | rad      | rad/s              | rad/s    | rad         | rad/s          | Nm           | Nm                 | Nm          | Nm                     | A                      |
```
