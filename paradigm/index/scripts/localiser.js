/**
 * based on jspsych-html-button-response
 *
 **/
/*jshint esversion: 6 */

jsPsych.plugins.localiser = (function() {

  var plugin = {};

  plugin.info = {
    name: 'localiser',
    description: '',
    parameters: {
      stimuli: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Stimuli',
        default: "0",
        array: true,
        description: 'Image file names for the stimuli'
      },
      stimuli_fname: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Stimuli filename',
        default: "0",
        array: true,
        description: 'Image file names for the stimuli (with full filepath)'
      },
      stim_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Stimulus number',
        default: null,
        description: 'Stimulus number (0 or 1)'
      },
      path: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Path number',
        default: null,
        description: 'Path number (0 or 1)'
      },
      trial_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Trial number',
        default: 0,
        description: 'Trial number in localiser block.'
      },
      block_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Block number',
        default: 0,
        description: 'Block number in localiser.'
      },
      img_duration: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Image duration',
        array: true,
        default: null,
        description: 'Duration of image presentation (uniform distribution)'
      },
      isi: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Inter-stimulus interval',
        default: 0,
        description: 'Time between stimuli, where screen is blank.'
      },
      trial_end_type: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Trial end type',
        default: "response",
        description: '"time-out" = ends at trial_duration, "time-or-response" = ends at trial_duration unless pressed earlier, "response" = ends only when button is pressed, "min-time" = ends at img_duration, or longer if a response is not made before img_duration'
      },
      words: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Words',
        array: true,
        default: null,
        description: 'The two word options to describe the image'
      },
      choices: {
        type: jsPsych.plugins.parameterType.KEYCODE,
        array: true,
        pretty_name: 'Choices',
        default: jsPsych.ALL_KEYS,
        description: 'The keys the subject is allowed to press to respond to the stimulus.'
      }
    }
  };

  plugin.trial = function(display_element, trial) {

    document.getElementById('jspsych-content').classList = 'jspsych-content my-container';
    document.getElementsByClassName("stars")[0].style.opacity = "0";
    document.getElementsByClassName("twinkling")[0].style.opacity = "0";

    console.log("LOCALISER: BLOCK " + (trial.block_num+1) + ", TRIAL " + (trial.trial_num+1));

    const parameters = loadParameters();
    trial.choices = Array.isArray(trial.choices) ? trial.choices.flat(Infinity) : trial.choices;

    var this_isi = trial.isi.length > 1 ? jsPsych.randomization.shuffle(window.linspace(trial.isi[0], trial.isi[1], 0.1))[0] : trial.isi;
    var this_img_dur = trial.img_duration > 1 ? jsPsych.randomization.shuffle(window.linspace(trial.img_duration[0], trial.img_duration[1], 0.1))[0] : trial.img_duration;

    // img information
    var this_state = trial.stim_num;
    var this_path = trial.path;
    var start_time = null;

    var triggers = [];

    var response = {
      rt: null,
      trial: trial.trial_num,
      block: trial.block_num,
      img: trial.stimuli,
      state: this_state,
      path: this_path,
      isi: this_isi,
      img_duration: this_img_dur,
      words: trial.words,
      accuracy: null
    };

    // display image
    display_element.innerHTML = '<img src="' + trial.stimuli_fname + '" class="localiser-img" alt="' + trial.stimuli + '"><div id="photodiode" style="animation-duration: 0.5s;"></div>';
    triggers.push(["image_onset",performance.now()]);

    // switch to image
    jsPsych.pluginAPI.setTimeout(function() {
      show_cue();
    }, this_img_dur * 1000);

    var html = '';
    function show_cue() {

      html = '<div style="width:80vw; display:flex; flex-direction:row; justify-content:space-evenly;"><p class="localiser-cue">' + trial.words[0] + '</p><p class="localiser-cue">' + trial.words[1] + '</p></div>';
      if (parameters.exp_variables.meg_mode) {
        html += '<div id="photodiode" style="animation-duration: 0.1s;"></div>';
      }

      display_element.innerHTML = html;
      triggers.push(["cue_onset",performance.now()]);
      start_time = performance.now();

      var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
        callback_function: after_response,
        valid_responses: trial.choices,
        rt_method: 'performance',
        persist: false,
        allow_held_key: false
      });

    }

    function after_response(info) {

      jsPsych.pluginAPI.cancelAllKeyboardResponses();
      response.rt = info.rt;
      response.button = info.key;

      // display feedback
      var correct_key_name = trial.words[0] == trial.stimuli ? trial.choices[0] : trial.choices[1];
      var this_key = jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(info.key);
      var answer_correct = correct_key_name == this_key ? true : false;

      response.acc = answer_correct;

      var this_delay = 0;
      if (trial.trial_end_type == "min-time") {
        var current_offset = performance.now() - start_time;
        this_delay = (trial.img_duration * 1000) - current_offset;
        this_delay = this_delay < 0 ? 0 : this_delay;
      }

      display_element.innerHTML = answer_correct ? '<p class="localiser-cue" style="color:rgb(6, 252, 46)">+</p>' : '<p class="localiser-cue" style="color:rgb(255, 0, 0)">+</p><div id="photodiode"  style="animation-duration: 0.3s;"></div>';
      triggers.push(["feedback_onset",performance.now()]);

      jsPsych.pluginAPI.setTimeout(function() {
        end_trial();
      }, this_isi * 1000);

    }

    function end_trial() {

      response.triggers = triggers.flat(Infinity);

      // kill any remaining setTimeout handlers
      jsPsych.pluginAPI.clearAllTimeouts();
      // move on to the next trial
      jsPsych.finishTrial(response);

    }

  };

  return plugin;
})();
