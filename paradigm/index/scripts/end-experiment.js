/**
 * simple-button-text
 * adapted from jsPsych's html-button-response
 *
 * documentation: docs.jspsych.org
 *
 **/
/*jshint esversion: 6 */

jsPsych.plugins["end-experiment"] = (function() {

  var plugin = {};

  plugin.info = {
    name: 'end-experiment',
    description: '',
    parameters: {
      background: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Background',
        default: 'stars',
        description: 'Whether the background is black or stars.'
      }
    }
  };

  plugin.trial = function(display_element, trial) {

    document.getElementById('jspsych-content').classList = 'jspsych-content my-container exploration-test-main';
    if (trial.background == 'stars') {
      document.getElementsByClassName("stars")[0].style.opacity = "1";
      document.getElementsByClassName("twinkling")[0].style.opacity = "1";
    } else {
      document.getElementsByClassName("stars")[0].style.opacity = "0";
      document.getElementsByClassName("twinkling")[0].style.opacity = "0";
    }

    // display form
    var html = '<div id="end-container">' +
      '<h1>Finished!</h1>' +
      '<p>Thank you for completing the experiment. Your data will be analysed in the next few days and your bonus payment will be finalised. Before you go, please submit some feedback.</p>' +
      '<form id="end-survey">' +
      '<p style="font-size: 1rem;">Did you use a particular strategy in this game?</p>' +
      '<textarea name="strategy"></textarea>' +
      '<p style="font-size: 1rem;">Do you have any feedback for us (e.g. errors, visualisation, difficulty, etc.)?</p>' +
      '<textarea name="feedback"></textarea></form>' +
      '<div style="display:flex; justify-content:center;"><button type="submit" id="submit-complete" style="font-size: 1.3rem;">Submit</button></div></div>';

    display_element.innerHTML = html;

    // start time
    var start_time = performance.now();

    // store response
    var response = {
      rt: null,
      button: null
    };

    document.getElementById('submit-complete').addEventListener('click', btnListener);

    function btnListener(e) {
      e.target.removeEventListener('click', btnListener);
      after_response();
    }

    // function to handle responses by the subject
    function after_response() {

      // measure rt
      var end_time = performance.now();
      var rt = end_time - start_time;
      response.feedback = document.getElementById("end-survey")[1].value;
      response.strategy = document.getElementById("end-survey")[0].value;
      response.rt = rt;

      end_trial();

      }

    // function to end trial when it is time
    function end_trial() {

      // kill any remaining setTimeout handlers
      jsPsych.pluginAPI.clearAllTimeouts();

      // gather the data to store for the trial

      // clear the display
      display_element.innerHTML = '';

      // move on to the next trial
      jsPsych.finishTrial(response);
    }

  };

  return plugin;
})();
