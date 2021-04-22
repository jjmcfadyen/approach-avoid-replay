/**
 * based on jspsych-html-button-response
 *
 **/
/*jshint esversion: 6 */

jsPsych.plugins["value-test"] = (function() {

  var plugin = {};

  plugin.info = {
    name: 'value-test',
    description: '',
    parameters: {
      trial_duration: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Trial duration',
        default: null,
        description: 'How long to show the trial, in the event of no response.'
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
      last_test: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Last test',
        default: false,
        description: 'Last memory test trial for this rep'
      },
      this_val: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'This stimulus value',
        default: 0,
        description: 'How many points are gained or lost when visiting this state.'
      },
      this_seq: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'This sequence number',
        default: 0,
        description: 'Which sequence (0 or 1) from which to draw a probe state from.'
      },
      stimuli_info: {
        pretty_name: 'Stimuli information',
        default: null,
        description: 'Object with image file names for each state in each path.'
      },
      choices: {
        type: jsPsych.plugins.parameterType.KEYCODE,
        array: true,
        pretty_name: 'Choices',
        default: null,
        description: 'The keys the subject is allowed to press to respond to the stimulus.'
      }
    }
  };

  plugin.trial = function(display_element, trial) {

    console.log("-------------------------------");
    console.log("Tutorial 1 TEST trial " + trial.trial_num + 1);
    console.log("--- Probe door = " + trial.test_seq);

    document.getElementById('jspsych-content').classList = 'jspsych-content my-container exploration-test-main';
    document.getElementsByClassName("stars")[0].style.opacity = "1";
    document.getElementsByClassName("twinkling")[0].style.opacity = "1";

    const parameters = loadParameters();

    var door_names = ['DOOR 1', 'DOOR 2', 'SUPPLY ROOM'];

    // empty global variables
    var start_time_first = null;
    var start_time_second = null;
    var value_options = null;

    function generate_values(correct_val){
      var value_options = null;
      var generate_values = true;
      while (generate_values) {
        value_options = jsPsych.randomization.shuffle([correct_val, Math.round(Math.random() * 10)-5, Math.round(Math.random() * 10)-5]);
        // check for duplicates
        var duplicates = 0;
        for (let i = 0; i < value_options.length; i++){
          duplicates = value_options.filter(x => x == value_options[i]).length > 1 ? duplicates + 1 : duplicates;
        }
        generate_values = duplicates == 0 ? false : true;
      }
      return value_options;
    }

    var val_try_count = jsPsych.data.get().filter({
      trial_type: 'value-test'
    }).select('attempt_num').values.slice(-1)[0]; // get last attempt_num
    if (trial.trial_num == 0) {
      val_try_count = val_try_count == undefined ? 0 : val_try_count + 1;
    }

    var response = {
      trial: trial.trial_num,
      door_probe: door_names[trial.this_seq],
      attempt_num: val_try_count,
      first_img: null,
      first_seq: trial.this_seq,
      first_state: null,
      first_rt: null,
      first_resp: null,
      first_acc: null,
      first_btn: null,
      first_value_probes: null,
      second_seq: trial.this_seq,
      second_img: null,
      second_rt: null,
      second_state: null,
      second_resp: null,
      second_acc: null,
      second_btn: null,
      second_value_probes: null,
      rep_acc: null
    };

    var cumulative_vals = null;
    var top_html = '';
    var html = '';
    var correct_val = null;

    var cumulative_acc = '--';
    var acc_array = [jsPsych.data.get().filter({
      attempt_num: val_try_count,
      trial_type: 'value-test'
    }).select('first_acc').values, jsPsych.data.get().filter({
      attempt_num: val_try_count,
      trial_type: 'value-test'
    }).select('second_acc').values].flat(Infinity);
    if (acc_array.length > 0) {
      cumulative_acc = (acc_array.map(x => x == true).filter(x => x).length / acc_array.length * 100).toFixed(0) + '%';
    }

    var text_header = '';
    if (parameters.exp_variables.meg_mode == false){
      text_header = '<h1>Memory test</h1>';
    }

    // create HTML
    function question_html(question_number) { // question number is 1 or 2

      var probe_question = '';
      var this_state = '';
      if (question_number == 1) {
        this_state = jsPsych.randomization.shuffle([0, 1, 2])[0]; // which state to ask about
        response.first_state = this_state;
        probe_question = 'How many oxygen points did you gain/lose in this room?';
      } else {
        this_state = jsPsych.randomization.shuffle([0, 1, 2].filter(x => x != response.first_state))[0]; // don't ask about state already probed
        response.second_state = this_state;
        probe_question = 'After leaving this room, how many oxygen points would you have?<br>(i.e. what would the CUMULATIVE SUM be?)';
      }
      let this_img = trial.stimuli_info.seq_array[trial.this_seq][this_state]; // name of that state
      let this_img_filename = trial.stimuli_info.seq_array_filenames[trial.this_seq][this_state];
      let probe_image = '<img src="' + this_img_filename + '" alt="' + this_img + '" height="100%" width="auto">';
      if (question_number == 1){
        response.first_img = this_img;
      } else {
        response.second_img = this_img;
      }

      // work out correct answer
      if (question_number == 1) {
        correct_val = trial.this_val[trial.this_seq][this_state];
        response.first_value_probes = generate_values(correct_val);
        value_options = response.first_value_probes;
      } else {
        cumulative_vals = trial.this_val[trial.this_seq].reduce(function(r, a) {
          if (r.length > 0)
            a += r[r.length - 1];
          r.push(a);
          return r;
        }, []);
        correct_val = cumulative_vals[this_state];
        response.second_value_probes = generate_values(correct_val);
        value_options = response.second_value_probes;
      }

      var this_animation = 'animation: fadeIn 0.3s ease-in 0s forwards;';

      top_html = question_number == 1 ? '<div class="instructions-background" ' + this_animation + '><div id="value-content">' : '<div class="instructions-background"><div id="value-content" style="' + this_animation + '">';
      html = text_header +
        '<p id="question-text">' + probe_question + '</p>' +
        '<div id="value-image-container" style="display:flex; justify-content:center; animation:fadeIn 0.3s ease-in forwards;">' + probe_image + '</div><div class="value-btn-container">';

      for (let i = 0; i < value_options.length; i++) {
        html = value_options[i] == correct_val ? html + '<button id="btn-' + i + '" class="value-btn" data-choice="' + i + '">' : html + '<button id="btn-' + i + '" class="value-btn" data-choice="' + i + '">';
        html += value_options[i] + '</button>';
      }
      html += '</div><div class="value-btn-container"><p>&lt;</p><p>^</p><p>&gt;</p></div>';

      let bottom_html = '<p id="outcome-text" style="opacity:0; font-weight:600; text-align:center;">&nbsp</p></div></div>';

      // check whether to display memory performance or not
      if (question_number > 1) {
        acc_array.push(response.first_acc);
        cumulative_acc = (acc_array.map(x => x == true).filter(x => x).length / acc_array.length * 100).toFixed(0) + '%';
      }
      bottom_html += '<div id="memory-score"><p>Memory performance = ' + cumulative_acc + '</p></div>';

      // display on screen
      display_element.innerHTML = top_html + html + bottom_html;

      if (question_number == 1) {
        start_time_first = performance.now();
      } else {
        start_time_second = performance.now();
      }

      if (question_number == 1){
        // add event listeners to button
        if (trial.choices != null) {
          var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
            callback_function: after_response,
            valid_responses: trial.choices.flat(Infinity),
            rt_method: 'performance',
            persist: false,
            allow_held_key: false
          });
        }

        for (let i = 0; i < value_options.length; i++) {
          document.getElementById('btn-' + i).addEventListener('click', firstBtnListener);
        }
      } else {
        // add event listeners to button
        if (trial.choices != null) {
          var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
            callback_function: second_question,
            valid_responses: trial.choices.flat(Infinity),
            rt_method: 'performance',
            persist: false,
            allow_held_key: false
          });
        }

        for (let i = 0; i < value_options.length; i++) {
          document.getElementById('btn-' + i).addEventListener('click', secondBtnListener);
        }
      }

    }

    function firstBtnListener(e) {
      var choice = e.currentTarget.getAttribute('data-choice');
      response.first_rt = performance.now() - start_time_first;
      after_response(choice);
    }

    function secondBtnListener(e) {
      var choice = e.currentTarget.getAttribute('data-choice');
      response.second_rt = performance.now() - start_time_first;
      second_question(choice);
    }

    question_html(1);

    // collect answer from first 'submit' button
    function after_response(info) {

      jsPsych.pluginAPI.cancelAllKeyboardResponses();

      // measure response time
      response.first_rt = response.first_rt == null ? info.rt : response.first_rt;

      // save response
      if (info.length == 1) {
        response.first_resp = response.first_value_probes[info];
        response.first_btn = info;
      } else {
        var btn_code = jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(info.key);
        response.first_btn = trial.choices[0].includes(btn_code) ? "0" : trial.choices[1].includes(btn_code) ? "1" : "2";
        response.first_resp = response.first_value_probes[parseInt(response.first_btn)];
      }
      document.getElementById("btn-" + response.first_btn).style.background = 'rgb(0, 203, 178)';

      // determine accuracy
      response.first_acc = response.first_resp == correct_val;

      // show feedback
      var element = document.getElementById("outcome-text");
      if (response.first_acc) {
        element.style.color = 'rgb(73, 235, 59)';
        element.innerHTML = 'Correct';
      } else {
        element.style.color = 'rgb(255, 60, 60)';
        element.innerHTML = correct_val == 1 ? 'Incorrect. There was ' + correct_val : 'Incorrect. There were ' + correct_val;
      }
      element.style.animation = 'fadeIn 0.3s ease-in forwards';

      // fade out box
      document.getElementById("value-content").style.animation = "fadeOut 0.3s ease-out 1.5s forwards";

      jsPsych.pluginAPI.setTimeout(function() {
        question_html();
      }, trial.trial_duration * 1000);

    }


    function second_question(info) {

      jsPsych.pluginAPI.cancelAllKeyboardResponses();

      // measure response time
      response.second_rt = response.second_rt == null ? info.rt : response.second_rt;

      // save response
      if (info.length == 1) {
        response.second_resp = response.second_value_probes[info];
        response.second_btn = info;
      } else {
        var btn_code = jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(info.key);
        response.second_btn = trial.choices[0].includes(btn_code) ? "0" : trial.choices[1].includes(btn_code) ? "1" : "2";
        response.second_resp = response.second_value_probes[parseInt(response.second_btn)];
      }
      document.getElementById("btn-" + response.second_btn).style.background = 'rgb(0, 203, 178)';

      // determine accuracy
      response.second_acc = response.second_resp == correct_val;

      // show feedback
      var element = document.getElementById("outcome-text");
      if (response.second_acc) {
        element.style.color = 'rgb(73, 235, 59)';
        element.innerHTML = 'Correct';
        element.style.animation = 'fadeIn 0.3s ease-in forwards';
        document.getElementById("value-content").style.animation = "fadeOut 0.3s ease-out 1.5s forwards";
        jsPsych.pluginAPI.setTimeout(function() {
          endTrial();
        }, trial.trial_duration * 1000);
      } else {
        second_feedback();
      }

    }

    function second_feedback() {

      var probe_images = '';
      var these_state_values = '';
      var these_cumulative_values = '';
      for (let i = 0; i < 3; i++) {
        probe_images += '<img src="' + trial.stimuli_info.seq_array_filenames[trial.this_seq][i] + '" alt="' + trial.stimuli_info.seq_array[trial.this_seq][i] + '">';
        if (i == 0) {
          these_state_values += '<div style="display:flex; flex-wrap:nowrap; flex-direction:row; justify-content:space-between; width:33%;"><p style="text-align:left;">Gained/lost: </p><p style="width:100%; text-align:center;">' + trial.this_val[trial.this_seq][i] + '</p></div>';
          these_cumulative_values += '<div style="display:flex; flex-wrap:nowrap; flex-direction:row; justify-content:space-between; width:33%;"><p style="text-align:left;">Total sum: </p><p style="width:100%; text-align:center;"><strong>' + cumulative_vals[i] + '</strong></p></div>';
        } else {
          these_state_values += '<p style="width:33%; text-align:center;">' + trial.this_val[trial.this_seq][i] + '</p>';
          these_cumulative_values += '<p style="width:33%; text-align:center;"><strong>' + cumulative_vals[i - 1] + ' + ' + trial.this_val[trial.this_seq][i] + ' = ' + cumulative_vals[i] + '</strong></p>';
        }
      }

      acc_array.push(response.second_acc);
      cumulative_acc = (acc_array.map(x => x == true).filter(x => x).length / acc_array.length * 100).toFixed(0) + '%';

      var this_img_filename = trial.stimuli_info.seq_array[0].indexOf(response.second_img) != -1 ? trial.stimuli_info.seq_array_filenames[0][trial.stimuli_info.seq_array[0].indexOf(response.second_img)] : trial.stimuli_info.seq_array_filenames[1][trial.stimuli_info.seq_array[1].indexOf(response.second_img)];

      var small_buttons1 = '';
      var small_buttons2 = '';
      if (parameters.exp_variables.meg_mode == true){
        small_buttons1 = 'font-size:1.5rem;';
        small_buttons2 = 'style="font-size:1.5rem;"';
      }

      var btn_html = '<div class="value-btn-container">';
      for (let i = 0; i < response.second_value_probes.length; i++) {
        if (i == response.second_btn){
          btn_html += '<button id="btn-' + i + '" class="value-btn" data-choice="' + i + '" disabled style="background:rgb(0, 203, 178);' + small_buttons1 + '">' + response.second_value_probes[i] + '</button>';
        } else {
          btn_html += '<button id="btn-' + i + '" class="value-btn" data-choice="' + i + '" disabled ' + small_buttons2 + '>' + response.second_value_probes[i] + '</button>';
        }
      }
      btn_html += '</div>';

      var next_button = 'Next';
      if (parameters.exp_variables.meg_mode){
        next_button = '<span style="text-transform:lowercase;"><small> v </small></span> Next <span style="text-transform:lowercase;"><small> v </small></span>';
      }

      var small_image = '';
      if (parameters.exp_variables.meg_mode == true){
        small_image = 'style="height:150px; min-height:150px;"';
      }

      display_element.innerHTML = '<div class="instructions-background"><div id="value-content">' +
        text_header + '<p>After leaving this room, how many oxygen points would you have?<br>(i.e. what would the CUMULATIVE sum be?)</p>' +
        '<div id="value-image-container"' + small_image + '><img src="' + this_img_filename + '" alt="' + response.second_img + '" height="100%" width="auto"></div>' +
        btn_html +
        '<p id="outcome-text" style="opacity:0; animation: fadeIn 0.3s ease-in forwards; font-weight:600; color:rgb(255, 60, 60); text-align:center;">Incorrect. You would have ' + correct_val + '</p>' +
        '<div id="value-image-container2" style="animation: fadeIn 0.3s ease-in forwards;">' + probe_images + '</div>' +
        '<div id="value-container" style="display:flex; flex-direction:row; flex-wrap:nowrap; justify-content:space-around; opacity:0; animation: fadeIn 0.3s ease-in forwards;">' + these_state_values + '</div>' +
        '<div id="cumulative-value-container" style="display:flex; flex-direction:row; flex-wrap:nowrap; justify-content:space-around; opacity:0; animation: fadeIn 0.3s ease-in 1.5s forwards;">' + these_cumulative_values + '</div>' +
        '</div></div>' +
        '<button class="instructions-btn" id="next-button">' + next_button + '</button>' +
        '<div id="memory-score"><p>Memory performance = ' + cumulative_acc + '</p></div>';

        if (trial.choices != null) {
          if (parameters.exp_variables.meg_mode){
            var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
              callback_function: endTrial,
              valid_responses: parameters.key_responses.up_down_buttons[1],
              rt_method: 'performance',
              persist: false,
              allow_held_key: false
            });
          } else {
            var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
              callback_function: endTrial,
              valid_responses: trial.choices[2].flat(Infinity),
              rt_method: 'performance',
              persist: false,
              allow_held_key: false
            });
          }
        }

      document.getElementById('next-button').addEventListener('click', function(e) {
        endTrial()
      });

    }

    function endTrial() {
      jsPsych.pluginAPI.cancelAllKeyboardResponses();
      if (trial.last_test) {
        var rep_acc = [jsPsych.data.get().filter({
          attempt_num: val_try_count, trial_type: "value-test"
        }).select('first_acc').values, jsPsych.data.get().filter({
          attempt_num: val_try_count, trial_type: "value-test"
        }).select('second_acc').values].flat(Infinity);
        if (acc_array.length > 0) {
          rep_acc = acc_array.map(x => x == true).filter(x => x).length / acc_array.length;
        }
        response.rep_acc = rep_acc;
      }
      jsPsych.pluginAPI.clearAllTimeouts(); // kill any remaining setTimeout handlers
      jsPsych.finishTrial(response); // move on to the next trial
    }


  };

  return plugin;
})();
