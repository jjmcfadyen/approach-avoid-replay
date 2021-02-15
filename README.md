# Codebook

This is the codebook for an MEG study on approach-avoidance and replay. This work is currently in progress.

The scripts are split into three sections:
1. Experimental paradigm
2. Data preprocessing
3. Data analysis

## Experimental paradigm
The experiment was delivered via web browser (to allow for possible online data collection, although all data was ultimately collected in the lab).

### Creating experiment structures

Two structures, A and B, were created for the experiment. These were counterbalanced across participants.

The rules for an experiment structure were executed using the ***generateStructure.m*** script.

All four structures (A and B, behavioural and MEG) can be found as .js files in the paradigm folder.

### Running the experiment using jsPsych

The experiment was hosted on Google Firebase, using all scripts within the ***paradigm*** directory. The experiment was launched by the `index.html` file within the ***index*** subdirectory using Google Chrome. From there, the experiment was executed using scripts written in Javascript, supported by the jsPsych package.

Notes:

- Firebase files with security information are not included in this repository
- The experiment needed to be launched from a private or "incognito" browser tab, so that unique user IDs were generated. If multiple participants completed the experiment from a normal browser tab, only one user ID would be generated and thus the data would be overwritten.
- A hidden section at the top of the consent form (shown at initial launch by `index.html`) allows the experimenter to specify the subject number, the experiment structure to use (A or B), and whether the session was the first (the shorter behavioural session) or second (the longer MEG session, including an initial functional localiser task).

## Data preprocessing

### Behavioural data

Data was downloaded from the online Firestore database via the `exportData.py` script, which can be found under the ***preprocessing*** directory.

For the behavioural session, this produced **one .JSON file** per participant with their responses.

For the MEG session, this produced **one .JSON file for the main task** and **one .CSV file for the functional localiser.**

Note that, in the event that data had not been saved correctly to Firestore, data was also manually saved upon experiment completion by using the `jsPsych.data.get().csv()` command in the browser window.

Behavioural data was then cleaned via MATLAB in `analyseBehav.m`. Principal components analysis was conducted on the questionnaire scores using in `questionnaires_PCA.R`.

### MEG data

**Required MATLAB toolboxes**:
- SPM12
- Fieldtrip
- OSL
- FSL (not included due to large size)

The behavioural data required for this stage are the [SUBJECT ID]_fl.csv table (for the functional localiser) produced by the exportData python script, and the [SUBJECT ID]_task.csv table (for the main task) produced by the analyse_behav.m MATLAB script. The MEG data begins in raw CTF format: e.g., MG06139_RiskyRe_20200124_1.ds.

`preprocessMEG_1.m`: In the first stage, data is read in using Fieldtrip and manually cropped so that the recording ends shortly after the photodiode triggers stopped. The photodiode onsets are then identified and compared against the behavioural file. This is done semi-automatically, but note that several manual adjustments were required to remove incorrect onset/offset times.

`preprocessMEG_2.m`: In this script, the raw CTF data are imported using OSL and cropped. The trigger markers are then inserted as events. For the functional localiser, this included the image onset time with its image name. For the main task, this included the planning onset, transition screen onset, response onset, image onset (for approach trials only), and outcome screen onsets for each trial. The data are then filtered with a 0.5 Hz highpass filter, including a mains filter, and downsampled to 100 Hz. The data are then filtered with a 0.5 Hz highpass filter, including a mains filter, and downsampled to 100 Hz.
Notes:
* For subject 506599, the sensor positions appear to be missing from the first task session (MG06229_RiskyRe_20200303_05.ds). This means that bad channels couldn't be interpolated.
* This was done on a Windows machine, meaning the OSL toolboxes was edited to that it didn't load FSL (which needed a unix-based system).

`preprocessMEG_3.m`: This script conducts the ICA and epochs the data. This is the only component that requires FSL, and so needs to be done on a unix-based system. 

The final epoched .mat and .dat files are converted into Fieldtrip format and an associated information file is saved. These were then sent back to the work PC for further analyses.

Example files for one subject:

- 088674_FL_1_osl.mat
- 088674_FL_1_osl_info.mat
- 088674_FL_2_osl.mat
- 088674_FL_2_osl_info.mat *(... and so on for all 4 blocks)*
- 088674_Task_1_decision.mat
- 088674_Task_1_decision_info.mat *(... and so on for all 11 blocks)*
- 088674_Task_1_transition.mat
- 088674_Task_1_transition_info.mat *(... and so on for all 11 blocks)*
- 088674_Task_1_image.mat
- 088674_Task_1_image_info.mat *(... and so on for all 11 blocks)*
- 088674_Task_1_outcome.mat
- 088674_Task_1_outcome_info.mat *(... and so on for all 11 blocks)*
- 088674_Task_1_ICA.mat
- 088674_Task_1_ICA.dat *(... and so on, for FL and Task across all blocks)*

## Data analysis

### Creating classifiers
Classifiers were created by conducting L1 regularised logistic regression on the visually-evoked fields from the functional localiser. This was executed by a function, `q_buildClassifiers.m`, for each subject and each training time point on a cluster.

The accuracy of the classifiers was then assessed using leave-one-out cross-validation in `assess_classifiers.m` within the ***analysis*** directory.

Four training times (120, 130, 140, and 150 ms) produced classifiers with the top 80% cross-validation accuracy. Rather than select the best-performing classifier overall (130 ms) and use this to generate sequenceness for all participants during the planning window, we selected the classifier (out of the top 4) that produced the strongest overall sequenceness during planning for each participant (thus, allowing for a degree of between-subject variability). This was achieved in the `assess_optimisation.m` script in the ***analysis*** directory. The sequenceness was estimated by `ss_build.m` within the ***utils*** directory.

### Generating sequenceness

The `assess_sequenceness.m` script used the optimal classifiers from the previous step to generate sequenceness for each transition in each path during the planning period. The resultant sequenceness strength for each interval (or "lag") was then matched to the behavioural information of each trial in a large dataset (in long format) to be used in subsequent linear mixed effects modelling.

### Linear mixed effects modelling

The `lme.R` script under ***analysis*** reads in the final long-form dataset with all behavioural information and the sequenceness measures for each interval. Behavioural and behaviour-replay relationships are assessed using linear mixed effects modelling.
