/**
 * based on jspsych-html-button-response
 *
 **/
/*jshint esversion: 6 */

jsPsych.plugins["exploration-animation"] = (function() {
  var plugin = {};

  plugin.info = {
    name: "exploration-animation",
    description: "",
    parameters: {
      stimuli: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: "Stimuli",
        default: null,
        array: true,
        description: "Image label [0] and full filename [1]."
      },
      stim_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: "Stimulus number",
        default: 0,
        description: "Stimulus number (0 or 1)"
      },
      trial_num: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: "Trial number",
        default: 0,
        description: "Trial number in experiment."
      },
      trial_duration: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: "Trial duration",
        array: true,
        default: null,
        description: "Time-out duration of the trial (in seconds)."
      },
      stimuli_info: {
        pretty_name: "Stimuli information",
        default: null,
        description: "Object with image file names for each state in each path."
      },
      isi: {
        type: jsPsych.plugins.parameterType.Int,
        pretty_name: "Inter-stimulus interval",
        default: 0,
        description: "Time between stimuli, where screen is blank."
      },
      trial_end_type: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: "Trial end type",
        default: "response",
        description:
          '"time-out" = ends at trial_duration, "time-or-response" = ends at trial_duration unless pressed earlier, "response" = ends only when button is pressed'
      },
      choices: {
        type: jsPsych.plugins.parameterType.KEYCODE,
        array: true,
        pretty_name: "Choices",
        default: null,
        description:
          "The keys the subject is allowed to press to respond to the stimulus."
      },
      background: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Background',
        default: 'stars',
        description: 'Whether the background is black or stars.'
      }
    }
  };

  plugin.trial = function(display_element, trial) {
    document.getElementById("jspsych-content").classList =
      "jspsych-content my-container";
    document.getElementsByClassName("stars")[0].style.opacity = "1";
    document.getElementsByClassName("twinkling")[0].style.opacity = "1";

    // set up background
    document.getElementById('jspsych-content').classList = 'jspsych-content';

    if (trial.background == 'stars') {
      document.getElementsByClassName("stars")[0].style.opacity = "1";
      document.getElementsByClassName("twinkling")[0].style.opacity = "1";
    } else {
      document.getElementsByClassName("stars")[0].style.opacity = "0";
      document.getElementsByClassName("twinkling")[0].style.opacity = "0";
    }

    trial.choices = Array.isArray(trial.choices)
      ? trial.choices.flat(Infinity)
      : trial.choices;
    trial.isi =
      trial.isi.length > 1
        ? jsPsych.randomization.shuffle(
            window.linspace(trial.isi[0], trial.isi[1], 0.1)
          )[0]
        : trial.isi;

    if (trial.stim_num == 3) {
      trial.trial_end_type = "response";
    }

    // create prompt text
    var html = "";
    var button_names = ["Next Room >"];
    var room_ids = ["first", "next", "final"];
    var button_html = "";
    if (trial.stim_num < 3) {
      var prompt =
        "You enter the " +
        room_ids[trial.stim_num] +
        " room. You associate it with this image:";
      var state_img =
        '<div id="state-img-container"><img id="state-img" src="' +
        trial.stimuli[1] +
        '" alt="' +
        trial.stimuli[0] +
        '" style="animation: fadeIn 0.3s ease-in 0s forwards;"></div>';

      var container_fade =
        (trial.stim_num == 0) | (trial.stim_num == 3)
          ? "animation: fadeIn 0.3s ease-in 0s forwards;"
          : "";
      html =
        '<div class="instructions-background" style="grid-template-rows: 4; ' +
        container_fade +
        '">' +
        "<h1>Door " +
        (trial.seq_num + 1) +
        "</h1>" +
        '<p id="top-text" style="text-align:center; animation: fadeIn 0.3s ease-in 0s forwards;">' +
        prompt +
        "</p>" +
        state_img;
      html += "    </div>";
    } else {
      html =
        '<div class="instructions-background">' +
        '<p id="top-text" style="text-align:center; animation: fadeIn 0.3s ease-in 0s forwards;">That was the end of the path.<br><br>You go down a corridor that leads directly back to the <strong>CONTROL ROOM</strong></p>' +
        "</div>";
      button_names = ["Return >"];
    }

    if (
      trial.trial_end_type == "time-or-response" ||
      trial.trial_end_type == "response"
    ) {
      button_html =
        trial.stim_num == 0
          ? '<div id="buttons-container" style="animation: fadeIn 0.3s ease-in 0s forwards;">'
          : '<div id="buttons-container">';

      for (let i = 0; i < button_names.length; i++) {
        button_html +=
          '<button id="btn-continue" class="instructions-btn" style="font-family:\'Open Sans Condensed\'; text-transform:uppercase; margin: 5% 2% 0 2%; font-size:1.5em;">' +
          button_names[i] +
          "</button>";
      }
      button_html += "</div>";
    }

    // display stimulus
    html = '<div id="choice-container">' + html + button_html + "</div>";

    if (
      trial.trial_end_type == "time-out" ||
      trial.trial_end_type == "time-or-response"
    ) {
      jsPsych.pluginAPI.setTimeout(function() {
        end_trial(); // waits for fade out from previous function to end
      }, trial.trial_duration[1] * 1000);
    }

    display_element.innerHTML = html;
    // start time
    var start_time = performance.now();

    // add event listeners to button
    if (
      trial.trial_end_type == "time-or-response" ||
      trial.trial_end_type == "response"
    ) {
      if (trial.choices != null) {
        var keyboardListener = jsPsych.pluginAPI.getKeyboardResponse({
          callback_function: after_response,
          valid_responses: trial.choices,
          rt_method: "performance",
          persist: false,
          allow_held_key: false
        });
      }
      document
        .getElementById("btn-continue")
        .addEventListener("click", btnListener);
    }

    function btnListener(e) {
      response.rt = performance.now() - start_time;
      after_response();
    }

    var try_count = jsPsych.data
      .get()
      .select("attempt_num")
      .values.slice(-1)[0]; // get last attempt_num

    // store response
    var response = {
      trial: trial.trial_num,
      attempt_num: try_count,
      rt: null,
      img: trial.stimuli[0],
      isi: trial.isi,
      state: trial.stimuli_info.seq_array[trial.seq_num].indexOf(
        trial.stimuli[0]
      ),
      path: trial.seq_num
    };

    // function to handle responses by the subject
    function after_response(choice) {
      // clear all listeners
      jsPsych.pluginAPI.cancelAllKeyboardResponses();

      if (response.rt == null && choice.length != 1) {
        response.rt = choice.rt;
      }

      // remove all button listeners and disable buttons
      try {
        document
          .getElementById("btn-continue")
          .removeEventListener("click", btnListener);
        document
          .getElementById("btn-continue")
          .setAttribute("disabled", "disabled");
      } catch (err) {
        // console.log("catch!");
      }

      // hide image & text
      var fade_out = "fadeOut 0.3s ease-out 0s forwards";
      if (trial.stim_num >= 2) {
        // if it's the last stimulus, fade out the whole text window & button
        document.getElementById(
          "choice-container"
        ).style.animation = fade_out;
      } else {
        // otherwise, just fade out the text and the image
        document.getElementById("top-text").style.animation = fade_out;
        if (trial.stim_num < 3) {
          document.getElementById("state-img").style.animation = fade_out;
        }
      }

      jsPsych.pluginAPI.setTimeout(function() {
          end_trial(); // waits for fade out from previous function to end
        }, trial.isi * 1000);
    }

    // function to end trial when it is time
    function end_trial() {
      // kill any remaining setTimeout handlers
      jsPsych.pluginAPI.clearAllTimeouts();
      jsPsych.pluginAPI.cancelAllKeyboardResponses();

      // kill any key/click listeners
      if (
        trial.trial_end_type == "time-or-response" ||
        trial.trial_end_type == "response"
      ) {
        document
          .getElementById("btn-continue")
          .removeEventListener("click", btnListener);
        jsPsych.pluginAPI.cancelAllKeyboardResponses();
      }
      // move on to the next trial
      jsPsych.finishTrial(response);
    }
  };

  return plugin;
})();
