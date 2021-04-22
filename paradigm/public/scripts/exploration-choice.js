/**
 * based on jspsych-html-button-response
 *
 **/
 /*jshint esversion: 6 */

jsPsych.plugins["exploration-choice"] = (function() {

  var plugin = {};

  plugin.info = {
    name: 'exploration-choice',
    description: '',
    parameters: {
      choice: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Choice',
        default: "0",
        description: 'Forced choice (0 or 1)'
      },
      trial_duration: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Trial duration',
        default: null,
        description: 'How long to show the trial.'
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
    document.getElementById('jspsych-content').classList = 'jspsych-content my-container';
    if (parameters.exp_variables.meg_mode) {
      document.getElementsByClassName("stars")[0].style.opacity = "0";
      document.getElementsByClassName("twinkling")[0].style.opacity = "0";
    } else {
      document.getElementsByClassName("stars")[0].style.opacity = "1";
      document.getElementsByClassName("twinkling")[0].style.opacity = "1";
    }

    trial.choices = Array.isArray(trial.choices) ? trial.choices.flat(Infinity) : trial.choices;

    var door_names = ['DOOR 1','DOOR 2'];
    var button_prefix = ['',''];
    if (trial.choices != null){ // if a keyboard/button response is enabled (instead of just mouse clicks)
      button_prefix[0] = '\< ' + door_names[0];
      button_prefix[1] = door_names[1] + ' \>';
    }

    var try_count = jsPsych.data.get().select('attempt_num').values.slice(-1)[0]; // get last attempt_num
    if (trial.trial_num == 0){
      try_count = try_count == undefined ? 0 : try_count + 1;
    }

    // create prompt text
    var prompt = trial.trial_num == 0 ? 'You decide to go through ' : 'You arrive back in the <strong>CONTROL ROOM</strong> and go through the <strong>AIRLOCK</strong>.<br><br>This time, you decide to go through ';
    prompt += door_names[trial.choice] + '.';

    if (try_count > 0 && trial.trial_num == 0){
      prompt = 'Your memory hasn\'t quite recovered yet, so you try again to memorise the order of the rooms (represented by different images) behind ' + door_names[0] + ' and ' + door_names[1] + '.<br><br>' + prompt;
    }

    var html = '<div class="instructions-background"><p style="text-align:center;">' + prompt + '</p></div>';
    var button_html = '<div id="buttons-container">';
    for (let i = 0; i < door_names.length; i++){
      var allow_button = (i == trial.choice) ? '' : "disabled='disabled'";
      button_html += '<button id="door-' + i + '" class="instructions-btn" style="font-family:\'Open Sans Condensed\'; text-transform:uppercase; margin: 5% 2% 0 2%; font-size:1.5em;" ' + allow_button + ' data-choice=' + i + '>' + button_prefix[i] + '</button>';
    }
    button_html += '</div>';

    // display stimulus
    html = '<div id="choice-container" style="animation: fadeIn 0.3s ease-in 0s forwards">' + html + button_html + '</div>';

    if (parameters.exp_variables.meg_mode) {
      html += '<div id="photodiode"></div>';
    }
    display_element.innerHTML = html;

    // start time
    var start_time = performance.now();

    // add event listeners to button
    if (trial.choices != null){
      var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
        callback_function: after_response,
        valid_responses: trial.choices,
        rt_method: 'performance',
        persist: false,
        allow_held_key: false
      });
    }

    function btnListener(e){
      response.rt = performance.now() - start_time;
      var choice = e.currentTarget.getAttribute('data-choice'); // don't use dataset for jsdom compatibility
      response.button = choice;
      after_response(choice);
    }

    document.getElementById('door-' + trial.choice).addEventListener('click', btnListener);

    // store response
    var response = {
      trial: trial.trial_num,
      door_name: door_names[trial.choice],
      attempt_num: try_count,
      rt: null,
      button: null
    };

    // function to handle responses by the subject
    function after_response(choice) {

      var click_rt = performance.now() - start_time;

      // remove button listeners
      for (let i = 0; i < document.querySelectorAll('button').length; i++) {
        document.getElementById('door-' + i).removeEventListener('click', btnListener);
      }

      // after a valid response, the stimulus will have the CSS class 'responded'
      // which can be used to provide visual feedback that a response was recorded
      document.getElementById('door-' + trial.choice).className += ' responded';

      // disable button after response
      document.getElementById('door-' + trial.choice).setAttribute('disabled', 'disabled');

      if (choice.length != 1) { // keyboard response
        jsPsych.pluginAPI.cancelAllKeyboardResponses();
        response.rt = choice.rt;
        response.button = trial.choices.indexOf(jsPsych.pluginAPI.convertKeyCodeToKeyCharacter(choice.key));
      } else { // mouse click response
        response.rt = click_rt;
        response.button = choice;
      }

      end_trial();

    }

    // function to end trial when it is time
    function end_trial() {

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
