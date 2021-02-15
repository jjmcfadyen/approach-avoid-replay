/**
 * based on jspsych-html-button-response
 *
 **/
/*jshint esversion: 6 */

jsPsych.plugins["exploration-test"] = (function() {

  var plugin = {};

  plugin.info = {
    name: 'exploration-test',
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
        default: undefined,
        description: 'Trial number in experiment.'
      },
      test_seq: {
        type: jsPsych.plugins.parameterType.STRING,
        array: true,
        pretty_name: 'Test sequence',
        description: 'Which sequence is being tested'
      },
      last_test: {
        pretty_name: 'Last test',
        default: false,
        description: 'Last memory test trial for this rep'
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

    const parameters = loadParameters();
    document.getElementById('jspsych-content').classList = 'jspsych-content my-container exploration-test-main';
    if (parameters.exp_variables.meg_mode) {
      document.getElementsByClassName("stars")[0].style.opacity = "0";
      document.getElementsByClassName("twinkling")[0].style.opacity = "0";
    } else {
      document.getElementsByClassName("stars")[0].style.opacity = "1";
      document.getElementsByClassName("twinkling")[0].style.opacity = "1";
    }

    var door_names = ['DOOR 1', 'DOOR 2'];

    // empty global variables
    var start_time_first = null;
    var start_time_second = null;

    var try_count = jsPsych.data.get().filter({
      trial_type: 'exploration-animation'
    }).select('attempt_num').values.slice(-1)[0]; // get last attempt_num
    if (try_count == undefined) {
      try_count = jsPsych.data.get().filter({
        trial_type: 'exploration-test'
      }).select('attempt_num').values.slice(-1)[0];
      if (try_count == undefined && trial.trial_num == 0) {
        try_count = 0;
      } else if (try_count == undefined && trial.trial_num != 0) {
        throw ("Something went wrong with establishing the no. of attempts made with exploration (see exploration-test.js)");
      }
    }

    var response = {
      trial: trial.trial_num,
      door_probe: door_names[trial.test_seq[0]],
      attempt_num: try_count,
      first_prompts: [],
      first_rt: null,
      first_state: null,
      first_btn: null,
      first_acc: null,
      second_prompts: [],
      second_rt: null,
      second_state: null,
      second_btn: null,
      second_acc: null,
      rep_acc: null
    };

    // create HTML
    var prompt = '<span style="font-family:\'Open Sans Condensed\'; font-size: 1.5rem;">' + door_names[trial.test_seq[0]] + '</span>';

    // choose an example image from each path
    var probe_images = '';
    var probe_idx = jsPsych.randomization.shuffle([0, 1]);
    var correct_btn1 = probe_idx.indexOf(trial.test_seq[0]);

    var blank_counter = -1;
    var which_state = '';
    var first_state_names = [];
    var state_name = '';
    for (let i = 0; i < 2; i++) {

      which_state = trial.stimuli_info.seq_array_filenames[probe_idx[i]][trial.test_seq[1]];
      state_name = trial.stimuli_info.seq_array[probe_idx[i]][trial.test_seq[1]];
      first_state_names.push(state_name);

      let is_correct = i == correct_btn1;

      blank_counter = blank_counter + 2;

      probe_images += '<button class="state-btn" id="btn-blank-' + (blank_counter - 1) + '">' +
        '<img src="' + which_state + '" alt="blank" style="opacity:0;">' +
        '</button>' +
        '<button class="state-btn" id="btn1-' + i + '" state-name="' + state_name + '" acc="' + is_correct + '">' +
        '<img src="' + which_state + '" alt=""></button>' +
        '<button class="state-btn" id="btn-blank-' + blank_counter + '">' +
        '<img src="' + which_state + '" alt="blank" style="opacity:0;">' +
        '</button>';
      }


    // var blank_counter = -1;
    // for (let i = 0; i < trial.stimuli_info.seq_array.length; i++) {
    //   let which_state = jsPsych.randomization.shuffle([0, 1, 2])[0];
    //   let state_name = trial.stimuli_info.seq_array[probe_idx[i]][which_state];
    //   which_state = trial.stimuli_info.seq_array_filenames[probe_idx[i]][which_state];
    //   let is_correct = probe_idx[i] == trial.test_seq[0] ? true : false;
    //   if (is_correct) {
    //     correct_btn1 = i;
    //   }
    //   blank_counter = blank_counter + 2;
    //   probe_images += '<button class="state-btn" id="btn-blank-' + (blank_counter - 1) + '">' +
    //     '<img src="' + which_state + '" alt="blank" style="opacity:0;">' +
    //     '</button>' +
    //     '<button class="state-btn" id="btn1-' + i + '" state-name="' + state_name + '" acc="' + is_correct + '">' +
    //     '<img src="' + which_state + '" alt=""></button>' +
    //     '<button class="state-btn" id="btn-blank-' + blank_counter + '">' +
    //     '<img src="' + which_state + '" alt="blank" style="opacity:0;">' +
    //     '</button>';
    //   first_state_names.push(state_name);
    // }

    var html = '<div id="exploration-test-container" class="instructions-background" style="animation: fadeIn 0.3s ease-in 0s forwards; overflow-y:hidden;">' +
      '<h1>Memory test</h1>' +
      '<p id="top-text-1">If you went through ' + prompt + ', which of the rooms below would it lead you to?</p>' +
      '<p id="top-text-2">&nbsp</p>' +
      probe_images +
      '<p id="bottom-text" style="text-align:center;">Click on the image that you think is correct<br>(or press the corresponding UP or DOWN arrow key).</p>' +
      '</div>';

    // check whether to display memory performance or not
    var cumulative_acc = "--";
    var acc_array = [jsPsych.data.get().filter({
      attempt_num: try_count,
      trial_type: 'exploration-test'
    }).select('first_acc').values, jsPsych.data.get().filter({
      attempt_num: try_count,
      trial_type: 'exploration-test'
    }).select('second_acc').values].flat(Infinity);
    if (acc_array.length > 0) {
      cumulative_acc = (acc_array.map(x => x == "true").filter(x => x).length / acc_array.length * 100).toFixed(0) + '%';
    }
    html += '<div id="memory-score"><p>Memory performance = ' + cumulative_acc + '</p></div>';

    // display on screen
    display_element.innerHTML = html;

    start_time_first = performance.now();

    // log info
    for (let i = 0; i < probe_idx.length; i++) {
      let element = document.getElementById('btn1-' + i);
      let this_img = element.getAttribute('state-name');
      response.first_prompts.push('stim_' + trial.stimuli_info.seq_array[probe_idx[i]].indexOf(this_img) + '_' + this_img);
    }

    var correct_state1 = "";
    for (let i = 0; i < probe_idx.length; i++) {
      if (document.getElementById("btn1-" + i).getAttribute("acc") == "true") {
        correct_state1 = document.getElementById("btn1-" + i).getAttribute("state-name");
      }
    }

    // add event listeners to image buttons
    function btnListener1(e) {
      response.first_rt = performance.now() - start_time_first;
      var choice = e.currentTarget.getAttribute('id');
      choice = choice.substr(-1);
      response.first_acc = e.currentTarget.getAttribute('acc');
      response.first_btn = e.currentTarget.getAttribute('id');
      response.first_state = e.currentTarget.getAttribute('state-name');
      after_first_response(choice);
    }

    for (let i = 0; i < probe_idx.length; i++) {
      document.getElementById('btn1-' + i).addEventListener('click', btnListener1);
    }

    // add event listeners to button
    if (trial.choices != null) {
      var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
        callback_function: after_first_response,
        valid_responses: trial.choices[0].flat(Infinity),
        rt_method: 'performance',
        persist: false,
        allow_held_key: false
      });
    }

    function btnListener2(e) {
      response.second_rt = performance.now() - start_time_second;
      var choice = e.currentTarget.getAttribute('id');
      choice = choice.substr(-1);
      response.second_acc = e.currentTarget.getAttribute('acc');
      response.second_btn = e.currentTarget.getAttribute('id');
      response.second_state = e.currentTarget.getAttribute('state-name');
      after_second_response(choice);
    }

    function after_first_response(choice) {

      if (choice.length != 1) { // keyboard response
        jsPsych.pluginAPI.cancelAllKeyboardResponses();
        response.first_rt = choice.rt;
        response.first_btn = jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(choice.key);
        response.first_btn = trial.choices[0][0].includes(response.first_btn) ? "btn1-" + 0 : "btn1-" + 1;
        response.first_state = document.getElementById(response.first_btn).getAttribute("state-name");
        response.first_acc = document.getElementById(response.first_btn).getAttribute("acc");
      }

      for (let i = 0; i < probe_idx.length; i++) {
        var element = document.getElementById('btn1-' + i);
        element.setAttribute('disabled', 'disabled');
        if (i != parseInt(response.first_btn.substr(-1))) {
          element.style.animation = 'fadeOut 0.3s ease-out 0s forwards';
        } else {
          let imgStr = element.innerHTML.slice(0,-1);
          element.innerHTML = imgStr + ' style="border: 3px solid rgb(2, 240, 176);">';
        }
      }
      document.getElementById('bottom-text').style.animation = 'fadeOut 0.3s ease-out 0s forwards';

      jsPsych.pluginAPI.setTimeout(function() {
        after_response_feedback(); // waits for fade out from previous function to end
      }, 300);

    }

    function after_response_feedback() {

      var timer = 0;

      // fade in feedback
      var got_it_right = response.second_btn == null ? response.first_acc == "true" : response.second_acc == "true";

      document.getElementById("bottom-text").innerHTML = got_it_right ? '<span style="color:rgb(0,255,0);"><strong>Correct</strong></span>' : '<span style="color:rgb(255,0,0);"><strong>Incorrect</strong></span>';
      document.getElementById("bottom-text").style.animation = "fadeInOut 1.5s linear 0s forwards";
      timer += 1500;

      // show correct state if they got it wrong
      if (!got_it_right) {
        jsPsych.pluginAPI.setTimeout(function() {
          show_correct(); // waits for fade out from previous function to end
        }, timer);
        timer += 2000;
      }

      // move on to next question
      if (response.second_btn == null) {
        jsPsych.pluginAPI.setTimeout(function() {
          second_question(); // waits for fade out from previous function to end
        }, timer);
      } else {
        jsPsych.pluginAPI.setTimeout(function() {
          end_trial_fade(); // waits for fade out from previous function to end
        }, timer);

        timer += 300;
        jsPsych.pluginAPI.setTimeout(function() {
          end_trial(); // waits for fade out from previous function to end
        }, timer);
      }

    }

    function end_trial_fade() {
      document.getElementById("exploration-test-container").style.animation = "fadeOut 0.3s ease-out 0s forwards";
    }

    function show_correct() {

      // grey out answer
      var answer_element = response.second_btn == null ? document.getElementById(response.first_btn) : document.getElementById(response.second_btn);
      answer_element.style.animation = "greyOutImage 0.3s ease-out 0s forwards";
      var new_answer_element = answer_element.cloneNode(true);
      answer_element.parentNode.replaceChild(new_answer_element, answer_element);

      // replace bottom text
      var text_element = document.getElementById("bottom-text");
      text_element.innerHTML = "This was the correct answer.";
      var new_text_element = text_element.cloneNode(true);
      text_element.parentNode.replaceChild(new_text_element, text_element);

      // show correct answer
      for (let i = 0; i < 3; i++) {
        var element = response.second_btn == null ? document.getElementById('btn1-' + i) : document.getElementById('btn2-' + i);
        if (element != null) {
          if (element.getAttribute('acc') == "true") {
            let imgStr = element.innerHTML.slice(0,-1);
            if (response.second_btn == null) {
              element.innerHTML = imgStr + ' style="border: 3px solid rgb(2, 240, 176);">';
            }
            element.style.animation = 'fadeIn 0.3s linear 0s forwards';
          }
        }
      }

    }

    function second_question() {

      // update score
      acc_array.push(response.first_acc);
      if (acc_array.length > 0) {
        cumulative_acc = (acc_array.map(x => x == "true").filter(x => x).length / acc_array.length * 100).toFixed(0) + '%';
      }
      document.getElementById("memory-score").innerHTML = '<p>Memory performance = ' + cumulative_acc + '</p>';

      // if previous answer was incorrect, fade it out
      if (response.first_acc == "false") {
        document.getElementById(response.first_btn).style.animation = "fadeOutImage 0.3s ease-out 0s forwards";
      }

      // grey out first question
      document.getElementById("top-text-1").style.color = "hsla(0, 0%, 70%,.25)";

      // decide whether to ask about the door BEFORE or AFTER
      var adjacent_state = "";
      var adjacent_type = "";
      if (correct_state1 == trial.stimuli_info.seq_array[0][0] || correct_state1 == trial.stimuli_info.seq_array[1][0]) {
        adjacent_state = 1; // if first state, has to probe the next state
        adjacent_type = "NEXT";
      } else if (correct_state1 == trial.stimuli_info.seq_array[0][2] || correct_state1 == trial.stimuli_info.seq_array[1][2]) {
        adjacent_state = 1; // if last state, have to probe the previous state
        adjacent_type = "BEFORE";
      } else {
        adjacent_state = jsPsych.randomization.shuffle([0, 2])[0];
        adjacent_type = adjacent_state == 0 ? "BEFORE" : "AFTER";
      }

      var correct_state2 = trial.stimuli_info.seq_array_filenames[trial.test_seq[0]][adjacent_state];

      // fade in new question
      document.getElementById("top-text-2").style.animation = 'fadeIn 0.3s ease-in 0s forwards';
      document.getElementById("top-text-2").innerHTML = "Which of the following rooms came " + adjacent_type + "?";

      // fade in bottom text
      document.getElementById("bottom-text").innerHTML = "Click on the image that you think is corect<br>(or press the corresponding LEFT, UP/DOWN, or RIGHT arrow key).";
      document.getElementById("bottom-text").style.animation = "fadeIn 0.3s ease-in 0s forwards";

      // generate prompt images (2 stimuli from same path + 1 random other stimulus)
      var prompt_stimuli = trial.stimuli_info.seq_array[trial.test_seq[0]];
      prompt_stimuli = prompt_stimuli.filter(x => x != correct_state1);
      var other_seqs = [0, 1].filter(x => x != trial.test_seq[0]);
      other_seqs = trial.stimuli_info.seq_array[other_seqs[0]];
      for (let i = 0; i < first_state_names.length; i++) {
        other_seqs = other_seqs.filter(x => x != first_state_names[i]);
      }
      prompt_stimuli.push(jsPsych.randomization.shuffle(other_seqs)[0]);

      /////// make sure the OTHER sequence isn't one that was on the screen before


      var prompt_stimuli_fname = Array(prompt_stimuli.length);
      for (let i = 0; i < prompt_stimuli.length; i++) {
        prompt_stimuli_fname[i] = trial.stimuli_info.seq_array_filenames.flat(1).filter(x => x.includes(prompt_stimuli[i]));
      }
      prompt_stimuli_fname = jsPsych.randomization.shuffle(prompt_stimuli_fname);

      // decide on which row to present the prompt stimuli
      var used_row = "";
      for (let i = 0; i < probe_idx.length; i++) {
        if (document.getElementById('btn1-' + i).getAttribute('state-name') == correct_state1) {
          used_row = i;
          break;
        }
      }

      var prompt_row = used_row == 0 ? 1 : used_row == 1 ? 0 : 1;
      var blank_buttons = prompt_row == 0 ? [0, 1] : [2, 3];
      var replace_this_button = prompt_row == 0 ? "btn1-0" : "btn1-1";
      for (let i = 0; i < prompt_stimuli.length; i++) {

        var replace_element = i == 0 ? document.getElementById("btn-blank-" + blank_buttons[0]) : i == 1 ? document.getElementById(replace_this_button) : document.getElementById("btn-blank-" + blank_buttons[1]);

        var new_prompt_element = document.createElement('button');
        new_prompt_element.id = "btn2-" + i;
        new_prompt_element.innerHTML = '<img src="' + prompt_stimuli_fname[i] + '" alt="" width="100%">';
        replace_element.parentNode.replaceChild(new_prompt_element, replace_element);


        let is_correct = prompt_stimuli_fname[i] == correct_state2 ? true : false;
        document.getElementById("btn2-" + i).classList = "state-btn";
        document.getElementById("btn2-" + i).setAttribute("state-name", prompt_stimuli[i]);
        document.getElementById("btn2-" + i).setAttribute("acc", is_correct);
        document.getElementById("btn2-" + i).style.animation = "fadeIn 0.3s ease-in 0s forwards";

      }

      // log info
      for (let i = 0; i < prompt_stimuli.length; i++) {
        let element = document.getElementById('btn2-' + i);
        var this_img = element.getAttribute('state-name');
        var new_probe_idx = "";
        for (let j = 0; j < trial.stimuli_info.seq_array.length; j++) {
          if (trial.stimuli_info.seq_array[j].includes(this_img)) {
            new_probe_idx = j;
            break;
          }
        }
        response.second_prompts.push('stim_' + trial.stimuli_info.seq_array[new_probe_idx].indexOf(this_img) + '_' + this_img);
      }

      start_time_second = performance.now();

      for (let i = 0; i < prompt_stimuli.length; i++) {
        document.getElementById('btn2-' + i).addEventListener('click', btnListener2);
      }

      // add event listeners to button
      if (trial.choices != null) {
        var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
          callback_function: after_second_response,
          valid_responses: trial.choices[1].flat(Infinity),
          rt_method: 'performance',
          persist: false,
          allow_held_key: false
        });
      }

    }

    function after_second_response(second_choice) {

      if (second_choice.length != 1) { // keyboard response
        jsPsych.pluginAPI.cancelAllKeyboardResponses();
        response.second_rt = second_choice.rt;
        response.second_btn = jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(second_choice.key);
        response.second_btn = trial.choices[1][0].includes(response.second_btn) ? "btn2-" + 0 : trial.choices[1][1].flat(Infinity).includes(response.second_btn) ? "btn2-" + 1 : trial.choices[1][2].includes(response.second_btn) ? "btn2-" + 2 : undefined;
        response.second_state = document.getElementById(response.second_btn).getAttribute("state-name");
        response.second_acc = document.getElementById(response.second_btn).getAttribute("acc");
      }

      // remove all event listeners & disable buttons
      for (let i = 0; i < trial.stimuli_info.seq_array[trial.test_seq[0]].length; i++) {
        document.getElementById('btn2-' + i).removeEventListener('click', btnListener2);
        document.getElementById('btn2-' + i).setAttribute('disabled', 'disabled');
      }

      // fade out other buttons
      for (let i = 0; i < trial.stimuli_info.seq_array[trial.test_seq[0]].length; i++) {
        if (i != parseInt(response.second_btn.substr(-1))) {
          document.getElementById('btn2-' + i).style.animation = 'fadeOut 0.3s ease-out 0s forwards';
        }
      }
      document.getElementById('bottom-text').style.animation = 'fadeOut 0.3s ease-out 0s forwards';

      // update score
      acc_array.push(response.second_acc);
      if (acc_array.length > 0) {
        cumulative_acc = (acc_array.map(x => x == "true").filter(x => x).length / acc_array.length * 100).toFixed(0) + '%';
      }
      document.getElementById("memory-score").innerHTML = '<p>Memory performance = ' + cumulative_acc + '</p>';

      jsPsych.pluginAPI.setTimeout(function() {
        after_response_feedback(); // waits for fade out from previous function to end
      }, 300);

    }

    // function to end trial when it is time
    function end_trial(answer_correct) {

      if (trial.last_test) {
        var rep_acc = [jsPsych.data.get().filter({
          rep: trial.rep_num,
          trial_type: "exploration-test"
        }).select('first_acc').values, jsPsych.data.get().filter({
          rep: trial.rep_num,
          trial_type: "exploration-test"
        }).select('second_acc').values].flat(Infinity);
        if (acc_array.length > 0) {
          rep_acc = acc_array.map(x => x == "true").filter(x => x).length / acc_array.length;
        }
        response.rep_acc = rep_acc;
      }

      // kill any remaining setTimeout handlers
      jsPsych.pluginAPI.clearAllTimeouts();
      jsPsych.pluginAPI.cancelAllKeyboardResponses();

      // move on to the next trial
      jsPsych.finishTrial(response);

    }

    // end trial if time limit is set
    if (trial.trial_duration !== null) {
      jsPsych.pluginAPI.setTimeout(function() {
        end_trial();
      }, trial.trial_duration);
    }

  };

  return plugin;
})();
