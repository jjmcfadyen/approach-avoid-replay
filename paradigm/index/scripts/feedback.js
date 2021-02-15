/**
 * feedback
 * adapted from jspsych-html-button-response
 *
 **/
/*jshint esversion: 6 */

jsPsych.plugins.feedback = (function() {

  var plugin = {};

  plugin.info = {
    name: 'feedback',
    description: '',
    parameters: {
      feedback_type: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Feedback type',
        default: 'exploration-test-success',
        description: 'What sort of feedback this screen is for (i.e. which part of the experiment) - "exploration-test-success" or "block" or "endpoints"'
      },
      block: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Block number',
        default: 0,
        description: 'Block number that just finished'
      },
      timer_duration: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: 'Timer duration',
        default: 1,
        description: 'Duration of timer for block rest, in minutes'
      },
    }
  };

  plugin.trial = function(display_element, trial) {

    // set up background
    document.getElementById('jspsych-content').classList = 'jspsych-content';
    document.getElementsByClassName("stars")[0].style.opacity = "1";
    document.getElementsByClassName("twinkling")[0].style.opacity = "1";

    var html = '';

    if (trial.feedback_type == "exploration-test-success") {

      // get accuracy from the 'exploration-test' phase
      var this_accuracy = jsPsych.data.get().filter({
        trial_type: 'exploration-test'
      }).select('overall_accuracy').values[0]; // in decimals (e.g. --> 0.6)

      // display html
      html += '<div class="instructions-background">Well done! You scored ' + (this_accuracy * 100) + '%</div>';
      html += '<button class="instructions-btn">Continue</button>';

      display_element.innerHTML = html;

      // add event listener to 'Continue' button
      display_element.querySelector('button').addEventListener('click', function(e) {
        end_trial();
      });

    } else if (trial.feedback_type == "block") {
      var timer_text = '';
      if (trial.timer_duration > 0 && trial.timer_duration < 1){
        timer_text = "<p>You can take a " + trial.timer_duration + "-minute break. After " + trial.timer_duration + " minute, the next block will automatically start. Or...</p>";
      } else if (trial.timer_duration > 1) {
        timer_text = "<p>You can take a " + trial.timer_duration + "-minute break. After " + trial.timer_duration + " minutes, the next block will automatically start. Or...</p>";
      }
      html += "<div class='instructions-background'><p>End of block " + (trial.block+1) + " of " + num_blocks + "</p>" +
        timer_text +
        "<p>Press the button below to continue to the next block.</p></div>" +
        "<button class='instructions-btn'>Start next block</button>" +
        "<div class='timer-container'>" +
        "<div class='timer-spinner' style='animation-duration: " + trial.timer_duration*60 + "s'></div>" +
        "<div class='timer-filler' style='animation-duration: " + trial.timer_duration*60 + "s'></div>" +
        "<div class='timer-mask' style='animation-duration: " + trial.timer_duration*60 + "s'></div>" +
        "<div class='timer-ring'></div>" +
        "</div>";

        display_element.innerHTML = html;

        // start time-out, if applicable
        if (trial.timer_duration != 0 && trial.timer_duration != null) {
          var timeout = jsPsych.pluginAPI.setTimeout(function() {
            end_trial(); // waits for fade out from previous function to end
          }, trial.timer_duration*1000*60);
        }

        // listen for clicks
        display_element.querySelector('button').addEventListener('click', function(e) {
          end_trial();
        });

    } else if (trial.feedback_type == "endpoints"){

      // Calculate overall points
      var overall_points = jsPsych.data.get().select("sum_points").values;
      var overall_sum = overall_points.reduce((a, b) => a + b, 0);

      // Write HTML
      html += '<div class="instructions-background">Well done!<br><br>Overall, you earnt <strong>' + overall_sum + ' points</strong> out of a possible ' + overall_points.length + '</div>';
      html += '<button class="instructions-btn">Finish</button>';

    }

    // function to end trial when it is time
    function end_trial() {

      // kill any remaining setTimeout handlers
      jsPsych.pluginAPI.clearAllTimeouts();

      // clear the display
      display_element.innerHTML = '';

      // move on to the next trial
      jsPsych.finishTrial();
    }

  };

  return plugin;
})();
