// function loadParameters() --> creates 'parameters' variable containing aspects of interest
/*jshint esversion: 6 */

window.loadParameters = function(sesstype,structype) {

  var parameters = {
    "exp_variables": {
      "debug_mode": false,
      "meg_mode": sesstype=="meg",
      "structure_type": structype
    },
    "key_responses": {
      "left_right_buttons": [["3", "1"],[["leftarrow","ArrowLeft"],["rightarrow","ArrowRight"]]], // MEG (button-box), not MEG (keyboard)
      "up_down_buttons": [["4", "2"],[["uparrow","ArrowUp"],["downarrow","ArrowDown"]]], // MEG (button-box), not MEG (keyboard)
    },
    "trial_numbers": {
      "nTrlLocaliser": 120, // no. of trials in each functional localiser block
      "nBlockLocaliser": 4, // no. of functional localiser blocks
      "nRepsT1": [2,1], // how many times each state is repeated in exploration learning trials [non-debug, debug]
      "nRepsT2": [2,1], // how many times each state is repeated in value learning trials [non-debug, debug]
      "memory_test_reps": [3,1], // how many times each sequence is probed in exploration memory test [non-debug, debug]
      "value_test_reps": [2,1] // how many times each sequence is probed in value memory test [non-debug, debug]
    },
    "timing": {
      "localiser_img": 1, // duration of localiser image (in seconds)
      "localiser_isi": [0.5, 1.5], // range of ISI for functional localiser (in seconds)
      "exploration_isi": [0], // in seconds (same for tutorial 1 value learning)
      "exploration_img_duration": [3,4], // time image is on screen,
      "exploration_response_type": "time-out", // "time-out" = ends at trial_duration, "time-or-response" = ends at trial_duration unless pressed earlier, "response" = ends only when button is pressed (same for tutorial 1 value learning)
      "value_learning_response_type": "time-out",
      "choice_dur": 30, // duration of planning period (in seconds)
      "resp_window": 1, // IF participant's have to wait until end of planning window to respond, how long is the response window? (in seconds)
      "feedback_window": 1.2, // how long the "too slow!" message is shown for (in seconds)
      "test_animation_dur": [2, 2.5]
    },
    "thresholds": { // score cut-offs to be able to continue with experiment
      "learning_threshold": 0.8, // minimum score in exploration learning phase
      "max_tries": 3, // maximum number of times they can do the exploration learning phase
      "practice_threshold": 0.8 // minimum percentage of NON-MISSED trials in practice to continue to main experiment
    },
    "style": {
      "overall_background_colour": "radial-gradient(ellipse at bottom, #1B2735 0%, #090A0F 100%)"
    },
    "values": {
      "safe_val": [1], // array of possible values for the safe state (i.e. the supply room)
      "miss_penalty": 1 // point(s) lost if fail to make a response at choice phase (in practice & main test)
    },
    "state_images":{
      "all_states": ['baby', 'backpack', 'bicycle', 'bowtie', 'car', 'cat', 'cupcake', 'hourglass', 'house', 'lamp', 'toothbrush', 'zebra'],
      "file_ext": ".png"
    }
  };

  // change key responses depending on whether this is MEG (button box) or not (keyboard)
  parameters.exp_variables.questionnaires = parameters.exp_variables.meg_mode ? false : true;

  parameters.key_responses.left_right_buttons = parameters.exp_variables.meg_mode ? parameters.key_responses.left_right_buttons[0] : parameters.key_responses.left_right_buttons[1];
  parameters.key_responses.up_down_buttons = parameters.exp_variables.meg_mode ? parameters.key_responses.up_down_buttons[0] : parameters.key_responses.up_down_buttons[1];

  // change trial numbers depending on whether debug mode is selected or not
  parameters.trial_numbers.nRepsT1 = parameters.exp_variables.debug_mode ? parameters.trial_numbers.nRepsT1[1] : parameters.trial_numbers.nRepsT1[0];
  parameters.trial_numbers.nRepsT2 = parameters.exp_variables.debug_mode ? parameters.trial_numbers.nRepsT2[1] : parameters.trial_numbers.nRepsT2[0];
  parameters.trial_numbers.memory_test_reps = parameters.exp_variables.debug_mode ? parameters.trial_numbers.memory_test_reps[1] : parameters.trial_numbers.memory_test_reps[0];
  parameters.trial_numbers.value_test_reps = parameters.exp_variables.debug_mode ? parameters.trial_numbers.value_test_reps[1] : parameters.trial_numbers.value_test_reps[0];

  return parameters;

};
