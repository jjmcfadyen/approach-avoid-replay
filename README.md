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

* Firebase files with security information are not included in this repository
* The experiment needed to be launched from a private or "incognito" browser tab, so that unique user IDs were generated. If multiple participants completed the experiment from a normal browser tab, only one user ID would be generated and thus the data would be overwritten.
* A hidden section at the top of the consent form (shown at initial launch by `index.html`) allows the experimenter to specify the subject number, the experiment structure to use (A or B), and whether the session was the first (the shorter behavioural session) or second (the longer MEG session, including an initial functional localiser task).

## Data preprocessing

### Behavioural data

Data was downloaded from the online Firestore database via the `exportData.py` script, which can be found under the ***preprocessing*** directory.

For the behavioural session, this produced **one .JSON file** per participant with their responses.

For the MEG session, this produced **one .JSON file for the main task** and **one .CSV file for the functional localiser.**

Note that, in the event that data had not been saved correctly to Firestore, data was also manually saved upon experiment completion by using the `jsPsych.data.get().csv()` command in the browser window.

Behavioural data was then cleaned via MATLAB in `analyseBehav.m`.

### MEG data
