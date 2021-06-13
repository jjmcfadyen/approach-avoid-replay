import firebase_admin
from firebase_admin import credentials
from firebase_admin import firestore
import pandas as pd
import numpy as np
import json
import os
from IPython.core.debugger import set_trace

# --------------------------------------------------------------------------------------------------------
# Initialise Firebase
cred = credentials.Certificate("key.json")
firebase_admin.initialize_app(cred)

db = firestore.client()
# --------------------------------------------------------------------------------------------------------
# Define parameters

megFile = True

if megFile:
    sections = ['0_functional_localiser','1_exploration_data', '2_exploration_test', '3_value_learning', '4_value_test', '5_negator_learning','6_test']
    analyseQuestionnaires = False
else:
    sections = ['1_exploration_data', '2_exploration_test', '3_value_learning', '4_value_test', '5_negator_learning','6_test']
    analyseQuestionnaires = True
questionnaires = ['intolerance_of_uncertainty', 'risk', 'worry']
subjects = ['zdkMWBU4kFg96dOVJ0O2sp661Rr2']
collectionName = 'pilot_v2.4'
worry_categories = ['recreational', 'financial', 'social', 'healthsafety', 'ethical']

maxBlocks = 11 # maximum number of blocks you'd expect to find in any given section
maxTrials = 121 # maximum number of trials you'd expect to find in any given block in a section
# --------------------------------------------------------------------------------------------------------
# Read in subject data and save

for subject in subjects:

    subj = db.collection("iterations").document(collectionName).collection("subjects").document(subject)

    # Get subject info and rebuild timeline
    subjectID = int(subj.get().to_dict()["subjectID"])

    if analyseQuestionnaires:

        ################################################################################################
                                             # Questionnaries #
        ################################################################################################

        qTimeline = [];
        subjData = {'questionnaires': {}}
        for qName in questionnaires:

            qDoc = subj.collection("questionnaires").document(qName).get().to_dict()['trial_data']

            data = {
                 "responses": qDoc['answer_list'],
                 "timeTaken": qDoc['rt'],
                 "totalScore": qDoc['total_score'],
                 "timeStamp": qDoc['time_elapsed']
            }

            if ('categories' in qDoc):
                data = {
                 "responses": qDoc['answer_list'],
                 "categories": qDoc['categories'],
                 "timeTaken": qDoc['rt'],
                 "timeStamp": qDoc['time_elapsed'],
                 "totalScore": {'totalScore': qDoc['total_score'],
                                              'ethicalScore': qDoc['ethical_score'],
                                              'financialScore': qDoc['financial_score'],
                                              'healthScore': qDoc['healthsafety_score'],
                                              'recreationalScore': qDoc['recreational_score'],
                                              'socialScore': qDoc['social_score']}
            }

            qTimeline.append(qDoc['time_elapsed'])

            subjData['questionnaires'][qName] = data

        subjData['timeline'] = qTimeline
    else:
        subjData = {};

    ################################################################################################
                                            # Sections #
    ################################################################################################

    for sect in sections:

        thisSection = subj.collection(sect)
        print('Importing data from section ', sect, '...')

        # Create empty structures to fill in for loop
        sectionData = {'block': [],
                       'trial': [],
                       'timeStamp': []}

        if sect == '1_exploration_data' or sect == '3_value_learning':
            sectionData['choice'] = []
            sectionData['rt'] = []
            sectionData['stimuli'] = {
                               'imgs': [],
                               'isi': [],
                               'path': [],
                               'state': [],
                               'timeStamp': []
                               }

        elif sect == '2_exploration_test':
            for q in ['q1','q2']:
                sectionData[q] = { # i.e. which image was behind [door 1] / [door 2] ?
                                   'rt': [],
                                   'acc': [],
                                   'choice': [],
                                   'probes': []
                                   }
        elif sect == '4_value_test':
            sectionData['probedPath'] = [] # which path is being asked about
            sectionData['timeStamp'] = []
            for q in ['q1','q2']:
                sectionData[q] = {
                                   'rt': [],
                                   'acc': [],
                                   'choice': [], # which of the probeValues was chosen
                                   'probeImg': [], # what was the value of THIS image? (q1) OR what was the cumulative sum in this image ? (q2)
                                   'probeState': [], # image's state position
                                   'probeValues': [] # the options for which values to choose from
                                   }
        elif sect == '5_negator_learning' or sect == '6_test':
            sectionData['timeStamp'] = []
            sectionData['acc'] = []
#             sectionData['block_score'] = []
            sectionData['negators'] = []
            sectionData['probability'] = []
            sectionData['choice'] = []
            sectionData['rt'] = []
            sectionData['transition'] = []
            sectionData['stimuli'] = {
                'img': [],
                'isi': [],
                'timeStamp': [],
                'value': [],
                'path': [],
                'state': []
            }
        elif sect == '0_functional_localiser':
            sectionData['timeStamp'] = []
            sectionData['img'] = []
            sectionData['isi'] = []
            sectionData['rt'] = []
            sectionData['state'] = []
            sectionData['path'] = []
            sectionData['word1'] = []
            sectionData['word2'] = []
            sectionData['acc'] = []

        if sect == '3_value_learning':
            sectionData['stimuli'] = {
                               'imgs': [],
                               'isi': [],
                               'path': [],
                               'state': [],
                               'timeStamp': [],
                               'value': []
                               }

        for block in list(range(1,maxBlocks+1)):

            if block >= 10:
                blockstr = str(block)
            else:
                blockstr = '0' + str(block)

            for trial in list(range(1,maxTrials+1)):

                if trial < 10:
                    strTrial = '0' + str(trial)
                else:
                    strTrial = str(trial)

                try:

                    # Get relevant document (usually the 'choice' part)
                    if subject == 'j7y9DFh6gCPkGUaKeMKo0fT044j1':
                        if sect == '6_test':
                            doc = thisSection.document('block0' + str(block-1) + '_trial' + strTrial + '_choice').get().to_dict()['trial_data']
                        elif sect == '0_functional_localiser' or sect == '1_exploration_data' or sect == '3_value_learning' or sect == '5_negator_learning':
                            doc = thisSection.document('block0' + str(block) + '_trial' + strTrial + '_choice').get().to_dict()['trial_data']
                        elif sect == '2_exploration_test' or sect == '4_value_test':
                            doc = thisSection.document('block0' + str(block) + '_trial' + strTrial).get().to_dict()['trial_data']
                    else:
                        if sect == '0_functional_localiser' or sect == '1_exploration_data' or sect == '3_value_learning' or sect == '5_negator_learning' or sect == '6_test':
                            doc = thisSection.document('block' + blockstr + '_trial' + strTrial + '_choice').get().to_dict()['trial_data']
                        elif sect == '2_exploration_test' or sect == '4_value_test':
                            doc = thisSection.document('block' + blockstr + '_trial' + strTrial).get().to_dict()['trial_data']
                        if subject == 'j7y9DFh6gCPkGUaKeMKo0fT044j1':
                            if sect == '6_test':
                                doc = thisSection.document('block0' + str(block-1) + '_trial' + strTrial + '_choice').get().to_dict()['trial_data']

                    sectionData['block'].append(block)
                    sectionData['trial'].append(trial)

                    # Extract choice & stimuli (i.e. animation) information
                    if sect == '0_functional_localiser':
                        sectionData['timeStamp'].append(doc['time_elapsed'])
                        sectionData['img'].append(doc['img'])
                        sectionData['isi'].append(doc['isi'])
                        sectionData['rt'].append(doc['rt'])
                        sectionData['state'].append(doc['state'])
                        sectionData['path'].append(doc['path'])
                        sectionData['word1'].append(doc['words'][0])
                        sectionData['word2'].append(doc['words'][1])
                        sectionData['acc'].append(doc['acc'])

                    if sect == '1_exploration_data' or sect == '3_value_learning':
                        sectionData['choice'].append(doc['choice'])
                        sectionData['rt'].append(doc['rt'])
                        sectionData['timeStamp'].append(doc['time_elapsed'])

                        imgs = []
                        isis = []
                        paths = []
                        states = []
                        timeStamps = []
                        values = []

                        for stim in range(0,3):
                            if sect == '3_value_learning' and any("SUPPLY ROOM" in item for item in imgs):
                                pass
                            else:
                                stimDoc = thisSection.document('block' + blockstr + '_trial' + strTrial + '_stim' + str(stim)).get().to_dict()['trial_data']
                                imgs.append(stimDoc['img'])
                                isis.append(stimDoc['isi'])
                                paths.append(stimDoc['path'])
                                states.append(stimDoc['state'])
                                timeStamps.append(stimDoc['time_elapsed'])
                                if sect == '3_value_learning':
                                    values.append(stimDoc['value'])

                        sectionData['stimuli']['imgs'].append(imgs)
                        sectionData['stimuli']['isi'].append(isis)
                        sectionData['stimuli']['path'].append(paths)
                        sectionData['stimuli']['state'].append(states)
                        sectionData['stimuli']['timeStamp'].append(timeStamps)
                        if sect == '3_value_learning':
                            sectionData['stimuli']['value'].append(values)

                    elif sect == '2_exploration_test':
                        sectionData['timeStamp'].append(doc['time_elapsed'])
                        for q in ['q1','q2']:
                            sectionData[q]['rt'].append(doc[q + '_rt'])
                            sectionData[q]['acc'].append(doc[q + '_acc'])
                            sectionData[q]['choice'].append(doc[q + '_choice'])
                            sectionData[q]['probes'].append(doc[q + '_probe_imgs'])
                    elif sect == '4_value_test':
                        sectionData['probedPath'].append(doc['path'])
                        sectionData['timeStamp'].append(doc['time_elapsed'])
                        for q in ['q1','q2']:
                            sectionData[q]['acc'].append(doc[q + '_acc'])
                            sectionData[q]['rt'].append(doc[q + '_rt'])
                            sectionData[q]['choice'].append(doc[q + '_choice'])
                            sectionData[q]['probeImg'].append(doc[q + '_probe_img'])
                            sectionData[q]['probeState'].append(doc[q + '_probe_state'])
                            sectionData[q]['probeValues'].append(doc[q + '_probe_values'])
                    elif sect == '5_negator_learning' or sect == '6_test':
                        sectionData['timeStamp'].append(doc['time_elapsed'])

                        sectionData['acc'].append(doc['acc'])
#                         sectionData['block_score'].append(doc['block_score'])
                        sectionData['negators'].append(doc['negators'])
                        sectionData['probability'].append(doc['probability'])
                        sectionData['choice'].append(doc['resp_name'])
                        sectionData['rt'].append(doc['rt'])
                        sectionData['transition'].append(doc['transition'])

                        imgs = []
                        isis = []
                        paths = []
                        states = []
                        values = []
                        timeStamps = []

                        for stim in list(range(0,4)):
                            try:
                                stimDoc = thisSection.document('block' + blockstr + '_trial' + strTrial + '_stim' + str(stim)).get().to_dict()['trial_data']
                                imgs.append(stimDoc['img'])
                                isis.append(stimDoc['isi'])
                                paths.append(stimDoc['path'])
                                states.append(stimDoc['state'])
                                values.append(stimDoc['value'])
                                timeStamps.append(stimDoc['time_elapsed'])
                            except:
                                pass

                        sectionData['stimuli']['img'].append(imgs)
                        sectionData['stimuli']['isi'].append(isis)
                        sectionData['stimuli']['path'].append(paths)
                        sectionData['stimuli']['state'].append(states)
                        sectionData['stimuli']['value'].append(values)
                        sectionData['stimuli']['timeStamp'].append(timeStamps)

                    if sect == '6_test':
                        print(sect,'block',block,',trial',trial)

                except:
                    pass

        subjData[sect] = sectionData


    # Add in some metadata
    subjData['subjectInfo'] = {
        'subjectID': subj.get().to_dict()['subjectID'],
        'date': subj.get().to_dict()['date'],
        'time': subj.get().to_dict()['time']
    }
# --------------------------------------------------------------------------------------------------------
# Save data
for subject in subjects:

    subj = db.collection("iterations").document(collectionName).collection("subjects").document(subject)

    # Get subject info and rebuild timeline
    subjectID = int(subj.get().to_dict()["subjectID"])

    if megFile == False:
        with open(os.path.join('output',str(subjectID) + '_s1_data.json'), 'w') as fp:
            json.dump(subjData, fp)
    else:
        with open(os.path.join('output',str(subjectID) + '_s2_data.json'), 'w') as fp:
            json.dump(subjData, fp)
        # Save just the functional localiser as a CSV for MATLAB MEG analysis
        FL = pd.DataFrame({
            'Block': subjData['0_functional_localiser']['block'],
            'Trial': subjData['0_functional_localiser']['trial'],
            'Image': subjData['0_functional_localiser']['img'],
            'ISI': subjData['0_functional_localiser']['isi'],
            'RT': subjData['0_functional_localiser']['rt'],
            'Acc': subjData['0_functional_localiser']['acc'],
            'Path': subjData['0_functional_localiser']['path'],
            'State': subjData['0_functional_localiser']['state']
        })
        export_csv = FL.to_csv (r'output/' + str(subjectID) + '_fl.csv', index = None, header=True)
