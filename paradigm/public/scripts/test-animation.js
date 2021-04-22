/**
 * based on jspsych-html-button-response
 *
 **/
/*jshint esversion: 6 */

jsPsych.plugins["test-animation"] = (function() {

  var plugin = {};

  plugin.info = {
    name: 'test-animation',
    description: '',
    parameters: {
      stim_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Stimulus number',
        default: 0,
        description: 'Stimulus number (0 or 1)'
      },
      trial_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Trial number',
        default: 0,
        description: 'Trial number in experiment.'
      },
      stimuli_info: {
        pretty_name: 'Stimuli information',
        default: null,
        description: 'Object with image file names for each state in each path.'
      },
      is_practice: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Is practice',
        default: 0, // can be "1","2", or "0"
        description: 'Whether this is the practice (i.e. no "lights off") or not.'
      },
      this_val: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'This stimulus value',
        default: 0,
        description: 'How many points are gained or lost when visiting this state.'
      },
      negs: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Negator states',
        default: null,
        array: true,
        description: 'Which states in each path are negator states.'
      },
      trial_duration: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Trial duration',
        array: true,
        default: null,
        description: 'Time-out duration of the trial (in seconds).'
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
        description: '"time-out" = ends at trial_duration, "time-or-response" = ends at trial_duration unless pressed earlier, "response" = ends only when button is pressed'
      },
      block_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Block number',
        default: null,
        description: 'Block number in experiment.'
      },
      show_score: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Show score',
        default: false,
        description: 'Whether to show block score or not.'
      },
      last_trial: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: "Last trial",
        default: false,
        description: "Whether this is the last trial of the block or not."
      }
    }
  };

  plugin.trial = function(display_element, trial) {

    var door_names = ["DOOR 1", "DOOR 2", "SUPPLY ROOM"];

    const parameters = loadParameters();

    if (parameters.exp_variables.meg_mode) {
      document.getElementsByClassName("stars")[0].style.opacity = "0";
      document.getElementsByClassName("twinkling")[0].style.opacity = "0";
    } else {
      document.getElementsByClassName("stars")[0].style.opacity = "1";
      document.getElementsByClassName("twinkling")[0].style.opacity = "1";
    }

    var safe_val = jsPsych.randomization.shuffle(parameters.values.safe_val)[0];

    var save_response = true;

    trial.isi = trial.isi.length > 1 ? jsPsych.randomization.shuffle(window.linspace(trial.isi[0], trial.isi[1], 0.1))[0] : trial.isi;
    if (trial.stim_num > 0) {
      if (jsPsych.data.get().filter({
          trial_type: 'test-animation'
        }).select('path').values.slice(-1)[0] == 2) { // if this is the supply room
        trial.isi = jsPsych.data.get().select('isi').values.slice(-1)[0];
        save_response = false;
      }
    }

    var test_trial = trial.is_practice != 1 ? true : false;
    var outcome_val = null;
    var skip_animation = false;
    var negate_state = null;
    var cumsum_array = [];
    var previous_free = true;
    var change_score = true;
    var show_images = true;
    var plus_sign = '';
    var font_end = '';

    var animation_delay = [0.5, 1];

    // figure out which sequence this is
    trial.block_num = jsPsych.data.get().select('block').values.slice(-1)[0];
    var this_seq_name = "";
    var this_seq_num = null;
    var data = jsPsych.data.get().filter({
      trial: trial.trial_num,
      practice: trial.is_practice,
      block: trial.block_num,
      trial_type: 'test-choice'
    });

    if (data.select('rt').values.slice(-1)[0] != null) { // if they didn't MISS the choice phase

      trial.choice = data.select('door_name_transition').values.slice(-1)[0] == null ? data.select('door_name').values.slice(-1)[0] : data.select('door_name_transition').values.slice(-1)[0];
      trial.choice = door_names.indexOf(trial.choice);
      previous_free = data.select('choice_type').values.slice(-1)[0] == "free" ? true : false;

      if (previous_free) { // parameters.exp_variables.meg_mode == true
        show_images = false;
      }

      if (trial.choice == 2) {
        if (trial.stim_num == 0) {
          outcome_val = safe_val; // for displaying the individual state values
        } else {
          skip_animation = true;
        }
      } else if (trial.choice < 2) {

        trial.stimuli = trial.stimuli_info.seq_array[trial.choice];
        outcome_val = trial.this_val[trial.choice][trial.stim_num]; // for displaying the individual state values

        cumsum_array = trial.this_val[trial.choice];
        for (let st = 0; st < cumsum_array.length; st++){
          if (st > 0){
            cumsum_array[st] = cumsum_array[st-1] + cumsum_array[st];
          }
          if (trial.negs != null){
            if (trial.negs[trial.choice] == st && cumsum_array[st] % 2 != 0){
              cumsum_array[st] = cumsum_array[st] * [-1];
            }
          }
        }

        if (trial.is_practice != 1) {
          negate_state = trial.stim_num == trial.negs[trial.choice] && cumsum_array[trial.stim_num] % 2 != 0 ? true : false;
          if (negate_state) {
            trial.trial_duration[1] += 0.5;
          }
        }
      }


      if (trial.stim_num > 2) {
        if (trial.choice == 2) {
          outcome_val = safe_val;
        } else {
          outcome_val = cumsum_array[cumsum_array.length - 1];
        }
      }

    } else {
      skip_animation = true;
      save_response = false;
    }

    var trial_img = '';
    if (!skip_animation) {
      if (trial.choice == 2) {
        trial_img = "SUPPLY ROOM";
      } else if (trial.stim_num < 3) {
        trial_img = trial.stimuli[trial.stim_num];
      }
    }

    // store response
    var response = {
      trial: trial.trial_num,
      img: trial_img,
      outcome: outcome_val,
      practice: trial.is_practice,
      block: trial.block_num,
      isi: trial.isi,
      state: trial.stim_num,
      path: trial.choice,
      value: outcome_val
    };

    if (trial.stim_num == 3 || !save_response) {
      response = {
        practice: trial.is_practice
      };
    }

    response.last_trial = trial.last_trial;

    var triggers = [];

    if (skip_animation) { // if they chose the supply ROOM, then only run this ONCE - OR - if they missed the response, don't run
      try {
        document.getElementById("choice-container").style.animation = "fadeOut 0.3s ease-out 0s forwards";
      } catch(e) {

      }
      if (trial.stim_num == 3) {
        jsPsych.pluginAPI.setTimeout(function() {
          end_trial();
        }, (trial.isi * 1000)*1.5);
      } else {
        end_trial();
      }
    } else {

      // create prompt text
      var html = "";
      var button_names = ["Next Room"];
      var button_html = "";
      var plural = '';
      if (outcome_val == 0 || outcome_val > 1 || outcome_val < -1) {
        plural = 's';
      }
      if (trial.choice < 2 && trial.stim_num < 3) { // any presentation except the last (i.e. 'You return to control room...')

        var prompt = '';
        if (trial.is_practice == 1 || trial.is_practice == 2) {
          prompt = 'You enter the room.';
          if (trial.stim_num == 1) {
            prompt = 'You enter the next room.';
          } else if (trial.stim_num == 2) {
            prompt = 'You enter the last room.';
          }
        }

        var state_img = '';
        if (show_images) {
          state_img = '<div id="state-img-container"><img id="state-img" src="' + trial.stimuli_info.seq_array_filenames[trial.choice][trial.stim_num] + '" width="auto" height="100%" alt="' + trial.stimuli[trial.stim_num] + '" style="animation: fadeIn 0.3s ease-in forwards;"></div>';
        } else {
          // state_img = '<p class="neg-cue-animation" style="animation: fadeIn 0.3s ease-in forwards;">' + trial.stimuli[trial.stim_num] + '</p>';
          state_img = '<p class="neg-cue-animation" style="animation: fadeIn 0.3s ease-in forwards;">?</p>';
        }

        var container_fade = trial.stim_num == 0 | trial.stim_num == 3 ? 'animation: fadeIn 0.3s ease-in forwards;' : '';
        html = '<div class="choice-background" style="grid-template-rows: 4; ' + container_fade + '; text-align:center;">' +
          '<h2 style="text-align:center;">' + door_names[trial.choice] + '</h2 style="text-align:center;">';
        if (trial.is_practice == 1 || trial.is_practice == 2) {
          html += '<p id="top-text" style="animation: fadeIn 0.3s ease-in forwards; text-align:center;">' + prompt + '</p>';
        }
        html += state_img;

        // create value text
        var outcome_text1 = '';
        var outcome_text2 = '';
        var plus_sign = '';
        var outcome_color = '';
        if (!negate_state) {
          if (outcome_val > 0){
            plus_sign = '+';
            outcome_color = 'rgb(14,204,33)';
          } else if (outcome_val == 0){
            outcome_color = 'rgb(240,240,240)';
          } else if (outcome_val < 0){
            outcome_color = 'rgb(255,20,20)';
          }
          outcome_text1 = '<p id="bottom-text" style="opacity: 0; animation: fadeIn 0.3s ease-in ' + animation_delay[0] + 's forwards; color:' + outcome_color + '; text-align:center;"><strong>' + plus_sign + outcome_val + ' point' + plural + '.</strong></p>';
        } else {
          var tmp_color = 'rgb(255,255,255)';
          if (outcome_val > 1) {
            tmp_color = 'rgb(14,204,33)';
          }
          if (outcome_val > 0){
            plus_sign = '+';
            outcome_color = 'rgb(14,204,33)';
          } else if (outcome_val == 0){
            outcome_color = 'rgb(240,240,240)';
          } else if (outcome_val < 0){
            outcome_color = 'rgb(255,20,20)';
          }

          trial.trial_duration[0] = trial.trial_duration[0] + 2;

          // trial.trial_duration += 2; // add time the trial duration
          outcome_text1 = '<p id="bottom-text" style="opacity: 0; animation: fadeIn 0.3s ease-in ' + animation_delay[0] + 's forwards; color:' + outcome_color + '; text-align:center;"><strong>' + plus_sign + outcome_val + ' point' + plural + '.</strong></p>';
          outcome_text2 = '<p id="outcome-text" style="opacity: 0; animation: fadeIn 0.3s ease-in ' + (animation_delay[0]+1) + 's forwards; color:rgb(255, 216, 0); text-align:center;"><strong>HAZARD ROOM</strong>: Total sum is an <strong>odd number</strong>. <i><strong>Reversing sum...</strong></i></p>';
          trial.trial_duration[1] = trial.trial_duration[1] + 1;
        }
        html += outcome_text1 + outcome_text2 + '</div>';
      } else if (trial.choice == 2) {
        var plus_sign = outcome_val > 0 ? '<font color="#38ff46">+' : '';
        var font_end = outcome_val > 0 ? '</font>' : '';
        if (trial.is_practice != 0) {
          html = '<div class="choice-background"><h2 style="text-align:center;">SUPPLY ROOM</h2>' +
            '<p id="top-text" style="animation: fadeIn 0.3s ease-in 0s forwards; text-align:center;">You go to the supply room.</p><p id="outcome-text" style="opacity: 0; animation: fadeIn 0.3s ease-in ' + animation_delay[1] + 's forwards; text-align:center;"><strong>' + plus_sign + outcome_val + ' point' + plural + '.' + font_end + '</strong></p>' +
            '</div>';
          trial.trial_duration[1] = trial.trial_duration[1] + 6;
        } else {
          if (!previous_free && data.select('door_name').values.slice(-1)[0] != "SUPPLY ROOM") {
            html = '<div class="choice-background" style="border: 8px solid rgb(154, 72, 241);"><h2 style="text-align:center; color:rgb(154, 72, 241);">Forced Choice: SUPPLY ROOM</h2>';
            trial.trial_duration[1] = trial.trial_duration[1] + 1;
          } else {
            html = '<div class="choice-background"><h2 style="text-align:center;">SUPPLY ROOM</h2 style="text-align:center;">';
          }
          html += '<p id="top-text" style="animation: fadeIn 0.3s ease-in 0s forwards; text-align:center;"></p><p id="outcome-text" style="opacity: 0; animation: fadeIn 0.3s ease-in forwards; text-align:center;"><strong>You collected ' + plus_sign + outcome_val + ' point' + plural + '.' + font_end + '</strong></p>' +
            '</div>';
        }
        button_names = ["Return"];
      } else {
        plus_sign = outcome_val > 0 ? '<font color="#38ff46">+' : '';
        font_end = outcome_val > 0 ? '</font>' : '';
        html = '<div class="choice-background">' +
          '<p id="top-text" style="animation: fadeIn 0.3s ease-in 0s forwards; text-align:center;">You return to the control room with: <strong>' + plus_sign + outcome_val + ' point' + plural + '.' + font_end + '</strong></p>' +
          '</div>';
        change_score = false;
        if (trial.is_practice == 2){
          trial.trial_duration[0] = trial.trial_duration[0] + 1;
        }
      }

      if (trial.trial_end_type == "time-or-response" || trial.trial_end_type == "response") {

        button_html = trial.stim_num == 0 ? '<div id="buttons-container" style="animation: fadeIn 0.3s ease-in 0s forwards;">' : '<div id="buttons-container">';

        for (let i = 0; i < button_names.length; i++) {
          if (i == 0) {
            button_html += '<div>'; // enclose the 'LEFT' and 'RIGHT' buttons together
          }
          button_html += '<button id="btn-continue" class="instructions-btn" style="font-family:\'Open Sans Condensed\'; text-transform:uppercase; margin: 5% 2% 0 2%; font-size:1.5em;" disabled="disabled">' + button_names[i] + '</button>';
          if (i == 1) {
            button_html += '</div>';
          }
        }
        button_html += '</div>';
      }

      // check cumulative sum
      var cumulative_acc = 0;
      if (trial.stim_num > 0) {
        cumulative_acc = cumsum_array[trial.stim_num - 1];
      }
      var score_colour = cumulative_acc > 0 ? 'rgb(34, 255, 63)' : cumulative_acc < 0 ? 'rgb(240, 0, 0)' : '';

      var little_points = '';
      var little_sign = '';
      var little_val = '';
      if (change_score){
        little_points = '';
        little_sign = '';
        little_val = outcome_val;
        if (outcome_val != null) {
          little_points = '<p style="color:';
          if (little_val > 0) {
            little_points += 'rgb(34,255,63)';
            little_sign = '+';
          } else if (little_val < 0) {
            little_points += 'rgb(240,0,0)';
          } else if (little_val == 0) {
            little_points += 'rgb(0, 255, 194)';
          }
          little_points += '; animation: pointFlash 1s ease-out ' + animation_delay[1] + 's forwards; opacity:0;"><big>' + little_sign + little_val + '</big></p>';
        }
      }

      var score_html = '';
      score_html = trial.stim_num == 0 ? '<div id="carry-score" style="animation: fadeIn 0.3s ease-in forwards;">' : '<div id="carry-score">';
      var score_text = 'Total sum';
      if (trial.stim_num == 3){
        score_text = 'Outcome';
      }
      score_html += '<p>' + score_text + ': </p><p id="big-points" style="color:' + score_colour + ';"><strong><big> ' + cumulative_acc + '</big></strong></p>' + little_points + '</div>';

      var total_score = '';
      if (trial.show_score) {
        // keep total oxygen score (unchanged) in the top right corner
        var total_acc = "0";
        var acc_array = jsPsych.data.get().filter({
          trial_type: 'test-animation',
          practice: trial.is_practice
        }).select('outcome').values;
        if (trial.is_practice == 0) {
          acc_array = jsPsych.data.get().filter({
            block: trial.block_num
          }).select('outcome').values;
        }
        if (acc_array.length > 0) {
          total_acc = acc_array.filter(x => x != null).reduce((a, b) => a + b, 0);
        }
        score_colour = total_acc > 0 ? 'rgb(34, 255, 63)' : total_acc < 0 ? 'rgb(240, 0, 0)' : '';
        total_score += '<div id="oxygen-score"><p>Points collected = </p><p style="color:' + score_colour + ';"><strong><big> ' + total_acc + '</big></strong></p></div>';
      }

      // display stimulus
      var dhtml = '';
      if (button_html == ''){
        dhtml = '<div id="choice-container">' + html + button_html + score_html + '</div></div>' + total_score;
      } else {
        dhtml = '<div id="choice-container">' + html + score_html + button_html + '</div></div>' + total_score;
      }

      if (parameters.exp_variables.meg_mode) {
        dhtml += '<div id="photodiode" style="animation-duration: 0.5s;"></div>';
      }

      if (trial.is_practice == 0) {
        display_element.innerHTML = dhtml;
      } else {
        display_element.innerHTML = dhtml;
      }
      triggers.push(["image_onset", performance.now()]);

      // set timer to change big points when little points flash up after delay
      jsPsych.pluginAPI.setTimeout(function() {
        var this_sum = parseInt(cumulative_acc) + little_val;
        document.getElementById('big-points').innerHTML = '<strong><big> ' + this_sum + '</big></strong>';
        if (this_sum > 0) {
          document.getElementById('big-points').style.color = 'rgb(34, 255, 63)';
        } else if (this_sum < 0) {
          document.getElementById('big-points').style.color = 'rgb(240, 0, 0)';
        } else {
          document.getElementById('big-points').style.color = 'rgb(0, 255, 194)';
        }
      }, animation_delay[1] * 1000);

      if (negate_state){
        jsPsych.pluginAPI.setTimeout(function() {
          document.getElementById('big-points').innerHTML = '<strong><big>' + cumsum_array[trial.stim_num] + '</big></strong>';
          if (cumsum_array[trial.stim_num] > 0) {
            document.getElementById('big-points').style.color = 'rgb(34, 255, 63)';
          } else if (cumsum_array[trial.stim_num] < 0) {
            document.getElementById('big-points').style.color = 'rgb(240, 0, 0)';
          } else {
            document.getElementById('big-points').style.color = 'rgb(255, 255, 255)';
          }
        }, (animation_delay[1]+1.5) * 1000);
      }

      // add event listeners to button (after delay)
      if ((trial.trial_end_type == "time-or-response" || trial.trial_end_type == "response") || (trial.is_practice != 0 && trial.choice == 2)) {

        jsPsych.pluginAPI.setTimeout(function() {
          // document.getElementById('btn-continue').removeAttribute("disabled");
          // document.getElementById('btn-continue').addEventListener('click', btnListener);
          if (trial.choices != null) {
            var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
              callback_function: after_response,
              valid_responses: trial.choices,
              rt_method: 'performance',
              persist: false,
              allow_held_key: false
            });
          }
        }, animation_delay[1] * 1000);
      }
      if (trial.is_practice != 0 && trial.choice == 2){
        trial.trial_duration[1] = 10;
      }

      jsPsych.pluginAPI.setTimeout(function() {
        if (trial.stim_num == 3){
          document.getElementById('choice-container').style.animation = 'fadeOut 0.2s ease-out forwards';
        } else {
          try {
            document.getElementById('state-img').style.animation = 'fadeOut 0.2s ease-out forwards';
            document.getElementById('bottom-text').style.animation = 'fadeOut 0.2s ease-out forwards';
            document.getElementById('outcome-text').style.animation = 'fadeOut 0.2s ease-out forwards';
          } catch(e){}
        }
      }, trial.trial_duration[0] * 1000);

      if (trial.trial_end_type == "time-or-response" || trial.trial_end_type == "time-out") {
        var endtime = null;
        if (negate_state){
          endtime = (trial.trial_duration[1] + trial.isi) * 1000;
        } else {
          endtime = (trial.trial_duration[0] + trial.isi) * 1000;
        }
        jsPsych.pluginAPI.setTimeout(function() {
          end_trial();
        }, endtime);
      }
    }

    function btnListener(e) {
      after_response();
    }

    // function to handle responses by the subject
    function after_response() {

      // after a valid response, the stimulus will have the CSS class 'responded'
      // which can be used to provide visual feedback that a response was recorded
      document.getElementById('btn-continue').className += ' responded';

      // disable button after response
      document.getElementById('btn-continue').setAttribute('disabled', 'disabled');

      // hide image & text
      var fade_out = "fadeOut 0.3s ease-out 0s forwards";
      document.getElementById('choice-container').style.animation = fade_out;

      jsPsych.pluginAPI.setTimeout(function() {
        end_trial(); // waits for fade out from previous function to end
      }, trial.isi * 1000);
    }

    // function to end trial when it is time
    function end_trial() {

      // remove button listeners
      if (trial.trial_end_type == "time-or-response" || trial.trial_end_type == "response") {
        document.getElementById('btn-continue').removeEventListener('click', btnListener);
        jsPsych.pluginAPI.cancelAllKeyboardResponses();
      }

      if (trial.stim_num < 3){
        response.outcome = trial.choice == 2 ? 1 : cumsum_array[trial.stim_num];
      } else {
        response.outcome = trial.choice == 2 ? 1 : cumsum_array.pop();
      }

      if (parameters.exp_variables.meg_mode) {
        response.triggers = triggers.flat(Infinity);
      }

      // kill any remaining setTimeout handlers
      jsPsych.pluginAPI.clearAllTimeouts();
      // move on to the next trial
      jsPsych.finishTrial(response);

    }

  };

  return plugin;
})();
