/**
 * based on jspsych-html-button-response
 *
 **/
/*jshint esversion: 6 */

jsPsych.plugins["test-choice"] = (function() {

  var plugin = {};

  plugin.info = {
    name: 'test-choice',
    description: '',
    parameters: {
      choice: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Choice',
        default: "0",
        description: 'Forced choice (0 or 1) or free choice.'
      },
      imgs_or_words: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Images or words',
        default: "images",
        description: 'Whether images or words are shown for the negators.'
      },
      transition: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Transition',
        default: null,
        description: 'Success or Failure (based on first choice)'
      },
      open_prob: {
        type: jsPsych.plugins.parameterType.INT,
        array: true,
        pretty_name: 'Opening probabilities',
        default: null,
        description: 'The probabilities of each door opening.'
      },
      expected_value: {
        type: jsPsych.plugins.parameterType.INT,
        array: true,
        pretty_name: 'Expected value',
        default: null,
        description: 'The final value of each path (negated)'
      },
      is_practice: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Is practice',
        default: 0, // can be "1","2", or "0"
        description: 'Whether this is the practice (i.e. no "lights off") or not.'
      },
      trial_duration: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Trial duration',
        default: null,
        description: 'How long to show the trial.'
      },
      trial_negs: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Trial negators',
        array: true,
        default: null,
        description: 'Array (length = 2) of which state is the negator state for each path.'
      },
      end_pause: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'End pause',
        default: 0,
        description: 'How long to freeze frame after response.'
      },
      response_ends_trial: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Response ends trial',
        default: true,
        description: 'If true, then trial will end when user responds.'
      },
      trial_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Trial number',
        default: 0,
        description: 'Trial number in experiment.'
      },
      block_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Block number',
        default: null,
        description: 'Block number in experiment.'
      },
      feedback_duration: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Feedback duration',
        default: 0,
        description: 'How long to show the "too slow" screen for (in seconds).'
      },
      response_window: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Response window',
        default: 0,
        description: 'How long to wait for response after planning time (in seconds).'
      },
      miss_penalty: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Miss penalty',
        default: 1,
        description: 'How much oxygen points go down by if they are too slow to respond.'
      },
      control_alert: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Control alert',
        default: 0,
        description: 'Whether there has been a change to PROBABILITY or NEGATORS since the last trial.'
      },
      show_score: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Show score',
        default: false,
        description: 'Whether to show block score or not.'
      },
      choices: {
        type: jsPsych.plugins.parameterType.KEYCODE,
        array: true,
        pretty_name: 'Choices',
        default: null,
        description: 'The keys the subject is allowed to press to respond to the stimulus.'
      },
      stimuli_info: {
        pretty_name: 'Stimuli information',
        default: null,
        description: 'Object with image file names for each state in each path.'
      },
      meg_mode: {
        type: jsPsych.plugins.parameterType.BOOLEAN,
        pretty_name: 'MEG mode',
        array: true,
        default: false,
        description: 'Whether this is the MEG block or the behavioural block.'
      }
    }
  };

  plugin.trial = function(display_element, trial) {

    const tmp_parameters = trial.meg_mode ? loadParameters("meg","") : loadParameters("behav","");
    var door_names = ['DOOR 1', 'DOOR 2', 'SUPPLY ROOM'];
    var test_trial = trial.is_practice == 0 || trial.is_practice == 2 ? true : false;

    if (trial.is_practice == 0) {
      console.log("TEST");
    } else {
      console.log("PRACTICE");
    }
    console.log("Block ", trial.block_num);
    console.log("Trial ", trial.trial_num);
    console.log("EV: ", trial.expected_value);

    if (trial.background=="black" | trial.meg_mode==true) {
      document.getElementsByClassName("stars")[0].style.opacity = "0";
      document.getElementsByClassName("twinkling")[0].style.opacity = "0";
    } else {
      document.getElementsByClassName("stars")[0].style.opacity = "1";
      document.getElementsByClassName("twinkling")[0].style.opacity = "1";
    }

    // determine transition success
    if (trial.choice == "free") {
      trial.transition = jsPsych.randomization.sampleWithReplacement([0, 1], 1, trial.open_prob)[0];
    }
    console.log("Transition:", trial.transition);

    if (trial.is_practice != 1) {

      // get V from experiment history
      var pastV = [
        [0, 0, 0],
        [0, 0, 0]
      ];
      try {
        var pastStim = jsPsych.data.get().filter({
          trial_type: 'test-animation'
        }).select('img').values.length;
        for (let i = 0; i < pastStim; i++){
          let p = jsPsych.data.get().filter({trial_type: 'test-animation'}).select('path').values[i];
          let s = jsPsych.data.get().filter({trial_type: 'test-animation'}).select('state').values[i];
          let newV = jsPsych.data.get().filter({trial_type: 'test-animation'}).select('value').values[i];
          pastV[p][s] = pastV[p][s] + newV;
        }
      } catch (e) {
        // no previous trials
        // currentV = pastV;
      }

      // determine optimal action
      trial.best_choice = trial.expected_value > 1 ? 0 : trial.expected_value < 1 ? 1 : "either";

    }

    var val_try_count = null;
    if (trial.is_practice == 1) {
      val_try_count = jsPsych.data.get().filter({
        trial_type: "test-choice"
      }).select('block').values.slice(-1)[0]; // get last attempt_num
      if (trial.trial_num == 0) {
        val_try_count = val_try_count == undefined ? 0 : val_try_count + 1;
      }
      trial.block_num = val_try_count;
    }

    document.getElementById('jspsych-content').classList = 'jspsych-content my-container';

    // store response
    var response = {
      trial: trial.trial_num,
      door_name: null,
      practice: trial.is_practice,
      rt: null,
      button: null,
      block: trial.block_num,
      choice_type: trial.choice,
      best_choice: trial.best_choice,
      door_prob: trial.open_prob,
      negators: trial.trial_negs,
      path_outcomes: trial.expected_value,
      acc: null
    };

    var triggers = [];

    if (test_trial) {
      response.door_name_transition = null;
      response.outcome = null;
    }

    // make timer, if applicable
    var timer_html = '';
    if (trial.is_practice == 0) {
      timer_html = test_trial ? "<div class='timer-container'>" +
        "<div class='timer-spinner' style='animation-duration: " + trial.trial_duration + "s'></div>" +
        "<div class='timer-filler' style='animation-duration: " + trial.trial_duration + "s'></div>" +
        "<div class='timer-mask' style='animation-duration: " + trial.trial_duration + "s'></div>" +
        "<div class='timer-ring'></div>" +
        "<div id='timer-text'><p id='timer-text-p'>Go!</p></div>" +
        "</div>" : '0';
    }

    // create prompt text
    var prompt = '';
    if (trial.is_practice == 1 && val_try_count > 0 && trial.trial_num == 0) {
      prompt = 'Your memory hasn\'t quite recovered yet, so you try again to memorise how many points you gain/lose in each room.<br><br>';
    }
    if (trial.choice == "free") {
      prompt = '';
    } else {
      if (trial.is_practice == 1) {
        prompt += trial.choice < 2 ? 'You decide to go through ' + door_names[trial.choice] : 'You decide to go to the SUPPLY ROOM';
      }
    }

    var html = '<div class="choice-background">';

    if (trial.is_practice != 1) {

      // var neg_alert = trial.control_alert == 2 || trial.control_alert == 3 ? 'animation: alertFlash 1s linear forwards;' : '';
      // var prob_alert = trial.control_alert == 1 || trial.control_alert == 3 ? 'animation: alertFlash 1s linear forwards;' : '';
      var neg_alert = '';
      var prob_alert = '';

      var neg_html = '';
      var neg_flex = '';
      var neg_order = [0, 1]; //jsPsych.randomization.shuffle([[0,1],[1,0]])[0]; // randomise whether the door 1 or door 2 negs are shown left/right
      if ((trial.imgs_or_words == "images" & trial.choice != "free")) {
        neg_html = '<img src="' + trial.stimuli_info.seq_array_filenames[neg_order[0]][trial.trial_negs[neg_order[0]]] + '" width="45%" style="padding:1%;">' +
          '<img src="' + trial.stimuli_info.seq_array_filenames[neg_order[1]][trial.trial_negs[neg_order[1]]] + '" width="45%" style="padding:1%;">';
        neg_flex = 'row';
      } else {
        neg_html = '<p class="neg-cue">' + trial.stimuli_info.seq_array[0][trial.trial_negs[0]] + '</p><p class="neg-cue">' + trial.stimuli_info.seq_array[1][trial.trial_negs[1]] + '</p>';
        neg_flex = 'column';
      }

      html += '<div id="warning-screen" style="display:grid; grid-template-rows: auto auto; grid-template-columns: auto auto;">' +
        '<h3 style="text-align:center;">Warning indicators:</h3>' +
        '<h3 style="text-align:center;">Door probability:</h3>' +
        '<div style="' + neg_alert + ' display:flex; justify-content:center; flex-direction:' + neg_flex + ';">' + neg_html + '</div>' +
        '<div style="height:45%; width:80%; display:flex; align-items:center; flex-direction:column; margin:auto; justify-content:space-evenly; ' + prob_alert + '">' +
        '<div id="door-1-prob-text" style="margin: 7% 0 7% 0;"><p style="text-align:center; margin:10px;"><big>Door 1: ' + Math.round(trial.open_prob[0] * 100) + '%</big></p>' +
        '<div id="door-1-prob-bar">' +
        '<div style="width:100%; height:8px; border:3px solid rgb(255, 230, 0); padding:0; margin:0;">' +
        '<div style="width:' + (trial.open_prob[0] * 100) + '%; height:8px; background:rgb(255, 230, 0); padding:0; margin:0;">' + '</div>' +
        '</div>' +
        '</div>' +
        '</div>' +
        '<div id="door-2-prob-text" style="margin: 7% 0 7% 0;"><p style="text-align:center; margin:10px;"><big>Door 2: ' + Math.round(trial.open_prob[1] * 100) + '%</big></p>' +
        '<div id="door-2-prob-bar">' +
        '<div style="width:100%; height:8px; border:3px solid rgb(255, 230, 0); padding:0; margin:0;">' +
        '<div style="width:' + (trial.open_prob[1] * 100) + '%; height:8px; background:rgb(255, 230, 0); padding:0; margin:0;">' + '</div>' +
        '</div>' +
        '</div>' +
        '</div>' +
        '</div>' +
        '</div>';
    }

    if (trial.is_practice == 2) {
      if (trial.choice == "free") {
        html += '<p id="prompt" style="text-align:center;">' + prompt + '</p><p style="text-align:center;">(Before you make your choice, think about whether you are more likely to get more points by going through the Airlock or going to the Supply Room)</p></div>';
      } else {
        html += '<p id="prompt" style="text-align:center;">' + prompt + '</p><p style="text-align:center;">(Practice calculating how many points you would get from Door 1 vs. Door 2)</p></div>';
      }
    } else {
      html += '<p id="prompt" style="text-align:center;">' + prompt + '</p></div>';
    }

    var button_html = '<div id="buttons-container">';
    var button_label = ''; // to append <, v, >
    for (let i = 0; i < 3; i++) {
      var allow_button = '';
      if (test_trial || (trial.is_practice == 1 && i != trial.choice)) {
        allow_button = "disabled='disabled'";
      }
      var door_prob = '';
      if (trial.is_practice == 1) {
        button_label = i == 0 ? '&lt; ' + door_names[0] : i == 1 ? '^ ' + door_names[1] + ' ^' : door_names[2] + ' &gt;';
        button_html += '<button id="door-' + i + '" class="instructions-btn" ' + allow_button + ' data-choice=' + i + '>' + button_label + door_prob + '</button>';
      } else if (i < 2) {
        button_label = i == 0 ? '&lt; AIRLOCK' : 'SUPPLY ROOM &gt;';
        if (i == 1 && (trial.is_practice == 2 || trial.is_practice == 0) && trial.choice == 0) {
          allow_button = "disabled='disabled'";
        }
        button_html += '<button id="door-' + i + '" class="instructions-btn" ' + allow_button + ' data-choice=' + i + '>' + button_label + '</button>';
        if (i == 0) {
          button_html += timer_html;
        }
      }
    }
    button_html += '</div>';

    var total_score = '';
    // check cumulative sum
    var cumulative_acc = 0;
    var acc_array = [];
    var data;
    if (trial.trial_num > 0) {

      if (trial.is_practice != 0) {
        data = jsPsych.data.get().filter({
          practice: trial.is_practice,
          choice_type: "free"
        });
      } else if (trial.is_practice == 0) {
        data = jsPsych.data.get().filter({
          block: trial.block_num,
          choice_type: "free",
          practice: trial.is_practice
        });
      }

      // correct negative values (i.e. make them twice as large, so that it's a loss, not zero)
      for (let i = 0; i < data.filter({
          trial_type: 'test-choice'
        }).select('trial').values.length; i++) {
        let trl = data.filter({
          trial_type: 'test-choice',
          block: trial.block_num
        }).select('trial').values[i];
        let trl_outcome = data.filter({
          trial_type: 'test-choice'
        }).select('door_name').values[i];
        let tmp = JSON.parse(jsPsych.data.get().filter({
          trial_type: 'test-animation',
          trial: trl,
          block: trial.block_num
        }).json());
        if (tmp.length == 0) { // missed response
          acc_array.push(-tmp_parameters.values.miss_penalty);
        } else {
          acc_array.push(tmp[tmp.length - 1].outcome);
        }
      }

      if (acc_array.length > 0) {
        cumulative_acc = acc_array.filter(x => x != null).reduce((a, b) => a + b, 0);
      }
    }

    if (trial.show_score && trial.choice == "free") {
      var score_colour = cumulative_acc > 0 ? 'rgb(34, 255, 63)' : cumulative_acc < 0 ? 'rgb(255, 10, 10)' : '';

      var little_points = '<p id="little-points" style="animation: fadeOut 0.3s ease-out ' + (trial.trial_duration + trial.response_window + 0.7) + 's forwards;">&nbsp</p>';

      total_score = '<div id="oxygen-score"><p>Points collected = </p><p id="big-points" style="color:' + score_colour + ';"><strong><big> ' + cumulative_acc + '</big></strong></p>' + little_points + '</div>';
    }

    // display stimulus
    var dhtml = '<div id="choice-container" style="animation: fadeIn 0.3s ease-in 0s forwards">' + html + button_html + '</div>' + total_score;
    if (trial.meg_mode) {
      dhtml += '<div id="photodiode"  style="animation-duration: 0.1s;"></div>';
    }
    display_element.innerHTML = dhtml;
    triggers.push(["trial_onset", performance.now()]);

    // start time
    var start_time = performance.now();

    // add event listeners to button
    var key_choices = [];
    if (trial.is_practice != 1 && trial.choice != "free") {
      key_choices = trial.choices[trial.choice];
    } else if (trial.is_practice != 1 && trial.choice == "free") {
      key_choices = trial.choices.flat();
    } else if (trial.is_practice == 1) {
      key_choices = trial.choices.slice(0);
    }
    if (trial.choices != null) {
      var response_delay = trial.is_practice ? 0 : 5000;
      jsPsych.pluginAPI.setTimeout(function() {
        var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
          callback_function: key_response,
          valid_responses: key_choices,
          rt_method: 'performance',
          persist: false,
          allow_held_key: false
        });}, response_delay);
    }

    function btnListener(e) {

      var choice = e.currentTarget.getAttribute('data-choice'); // don't use dataset for jsdom compatibility
      response.button = choice;

      // log response
      response.rt = performance.now() - start_time;
      response.button = choice;
      response.door_name = door_names[choice];
      if (trial.is_practice != 1) {
        response.door_name = choice == 0 ? 'Airlock' : door_names[2];
        if (trial.choice == "free") {
          if (choice == 1) {
            response.door_name_transition = response.door_name;
          } else {
            response.door_name_transition = door_names[trial.transition];
          }
        } else {
          if (trial.choice == 1) {
            response.door_name_transition = door_names[2];
          } else {
            response.door_name_transition = door_names[trial.transition];
          }
        }
      }

      // move on to next function
      if (trial.is_practice == 1 || (trial.is_practice != 1 && trial.choice == "free" && response.door_name == door_names[2]) || (trial.is_practice != 1 && trial.choice != "free" && response.door_name_transition == door_names[2])) {
        after_response();
      } else {
        after_response_transition();
      }
    }

    // set up button listeners
    if (trial.is_practice != 1) {
      jsPsych.pluginAPI.setTimeout(function() {
        for (let i = 0; i < 2; i++) {
          document.getElementById('door-' + i).addEventListener('click', btnListener);
          if ((i == trial.choice) || (trial.choice == 'free')) {
            document.getElementById('door-' + i).disabled = false;
          }
        }
      }, 5000);
    } else {
      for (let i = 0; i < 3; i++) {
        if (trial.choice == "free" || trial.choice == i) {
          document.getElementById('door-' + i).addEventListener('click', btnListener);
          document.getElementById('door-' + i).disabled = false;
        }
      }
    }


    // make timer disappear and 'GO!' text appear
    if (trial.is_practice == 0) {

      // 'GO!' timer
      jsPsych.pluginAPI.setTimeout(function() {
        try {
          document.getElementById('timer-text').style.opacity = 1;
        } catch (err) {}
      }, trial.trial_duration * 1000);

      // 'TOO SLOW!' timer (first)
      jsPsych.pluginAPI.setTimeout(function() {
        if (response.rt == null) {
          document.getElementById('timer-text-p').innerHTML = 'Too slow! <font color="red">-' + Math.abs(trial.miss_penalty) + '</font>';
          document.getElementsByClassName('timer-ring')[0].style.border = 'red solid 6px';
          response.outcome = trial.miss_penalty;
          if (trial.show_score && trial.choice == "free") {
            document.getElementById('little-points').innerHTML = '<font color="red"><big>' + trial.miss_penalty + '</big></font>';
            document.getElementById('big-points').innerHTML = '<strong><big> ' + (parseInt(cumulative_acc) + trial.miss_penalty) + '</big></strong>';
            if ((parseInt(cumulative_acc) + trial.miss_penalty) > 0) {
              document.getElementById('big-points').style.color = 'rgb(34, 255, 63)';
            } else if ((parseInt(cumulative_acc) + trial.miss_penalty) < 0) {
              document.getElementById('big-points').style.color = 'rgb(240, 0, 0)';
            } else {
              document.getElementById('big-points').style.color = 'rgb(0, 255, 194)';
            }
          }
          for (let i = 0; i < 2; i++) {
            document.getElementById('door-' + i).removeEventListener('click', btnListener);
            document.getElementById('door-' + i).setAttribute('disabled', 'disabled');
          }
        }
      }, trial.trial_duration * 1000 + trial.response_window * 1000);

      jsPsych.pluginAPI.setTimeout(function() {
        if (response.rt == null) {
          end_trial();
        }
      }, trial.trial_duration * 1000 + trial.response_window * 1000 + trial.feedback_duration * 1000);

    }

    function key_response(choice) {

      // remove button listeners
      for (let i = 0; i < 2; i++) {
        document.getElementById('door-' + i).removeEventListener('click', btnListener);
      }

      jsPsych.pluginAPI.cancelAllKeyboardResponses();
      var button = jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(choice.key);
      if (trial.is_practice == 1) {
        if (trial.choices.includes(button)) {
          button = trial.choice;
        }
      } else {
        for (let i = 0; i < trial.choices.length; i++) {
          if (trial.choices[i].includes(button)) {
            button = i;
          }
        }
      }
      response.button = button;

      // log response
      response.rt = choice.rt;
      if (trial.is_practice != 1) {
        response.door_name = button == 0 ? 'Airlock' : door_names[2];
        if (trial.choice == "free") {
          if (button == 1) {
            response.door_name_transition = response.door_name;
          } else {
            response.door_name_transition = door_names[trial.transition];
          }
        } else {
          if (trial.choice == 1) {
            response.door_name_transition = door_names[2];
          } else {
            response.door_name_transition = door_names[trial.transition];
          }
        }
      } else {
        response.door_name = door_names[trial.choice];
      }

      // move on to next function
      if (trial.is_practice == 1 || (trial.is_practice != 1 && response.door_name == door_names[2])) {
        after_response();
      } else {
        after_response_transition();
      }

    }

    // function to handle responses by the subject
    function after_response() {

      // remove button listeners
      for (let i = 0; i < 2; i++) {
        document.getElementById('door-' + i).removeEventListener('click', btnListener);
      }

      document.getElementById("choice-container").style.animation = "fadeOut 0.3s ease-out 0s forwards";
      jsPsych.pluginAPI.setTimeout(function() {
        end_trial();
      }, 300);
    }

    // function to handle responses if there is a probabilistic transition
    function after_response_transition() {

      // remove button listeners
      for (let i = 0; i < 2; i++) {
        document.getElementById('door-' + i).removeEventListener('click', btnListener);
      }

      document.getElementById("choice-container").style.animation = "fadeOut 0.3s ease-out 0s forwards";
      jsPsych.pluginAPI.setTimeout(function() {
        reveal_transition();
      }, 300);
    }

    function reveal_transition() {
      var reveal_html = '';
      if (trial.choice != "free") {
        reveal_html += '<div class="choice-background" style="animation: fadeInOut 2.5s linear forwards; border: 8px solid rgb(154,72,214);"><h3 style="text-align:center; color:rgb(154, 72, 241);">FORCED CHOICE</h3>';
      } else {
        reveal_html += '<div class="choice-background" style="animation: fadeInOut 2.5s linear forwards;">';
      }
      reveal_html += '<p style="text-align:center;">You go through the airlock.</p>' +
        '<p style="text-align:center;"><big><strong>' + door_names[trial.transition] + '</strong> is open.</big></p>' +
        '</div>' + total_score;
      reveal_html += trial.meg_mode ? '<div id="photodiode"  style="animation-duration: 0.3s;"></div>' : '';
      display_element.innerHTML = reveal_html;
      triggers.push(["doorReveal_onset", performance.now()]);
      jsPsych.pluginAPI.setTimeout(function() {
        end_trial();
      }, 2500);
    }

    // function to end trial when it is time
    function end_trial() {

      if (trial.meg_mode) {
        response.triggers = triggers.flat(Infinity);
      }

      response.acc = trial.best_choice == response.button;

      // kill any remaining setTimeout handlers
      jsPsych.pluginAPI.clearAllTimeouts();
      // move on to the next trial
      jsPsych.finishTrial(response);

    }

  };

  return plugin;
})();
