/*
 * edited version of jspsych-instructions.js
 */
 /*jshint esversion: 6 */

jsPsych.plugins["custom-instructions"] = (function() {

  var plugin = {};

  plugin.info = {
    name: 'custom-instructions',
    description: '',
    parameters: {
      pages: {
        type: jsPsych.plugins.parameterType.HTML_STRING,
        pretty_name: 'Pages',
        default: undefined,
        array: true,
        description: 'Each element of the array is the content for a single page.'
      },
      allow_meg: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Allow MEG continue',
        default: false,
        description: 'If in MEG mode, allow participant to move on from the instructions section'
      },
      key_forward: {
        type: jsPsych.plugins.parameterType.KEYCODE,
        pretty_name: 'Key forward',
        default: 'rightarrow',
        description: 'The key the subject can press in order to advance to the next page.'
      },
      key_backward: {
        type: jsPsych.plugins.parameterType.KEYCODE,
        pretty_name: 'Key backward',
        default: 'leftarrow',
        description: 'The key that the subject can press to return to the previous page.'
      },
      allow_backward: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Allow backward',
        default: true,
        description: 'If true, the subject can return to the previous page of the instructions.'
      },
      allow_keys: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Allow keys',
        default: true,
        description: 'If true, the subject can use keyboard keys to navigate the pages.'
      },
      button_label_previous: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Button label previous',
        default: 'Previous',
        description: 'The text that appears on the button to go backwards.'
      },
      button_label_next: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Button label next',
        default: 'Next',
        description: 'The text that appears on the button to go forwards.'
      },
      map: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Map',
        default: '',
        description: 'Visualisations of the map to put at the top-left of the screen.'
      },
      fade_from_black: {
        type: jsPsych.plugins.parameterType.BOOLEAN,
        pretty_name: 'Fade from black',
        default: false,
        description: 'Fade in the entire screen (including starry background) from black (useful at start of experiment, or after a section change).'
      },
      save_data: {
        type: jsPsych.plugins.parameterType.BOOLEAN,
        pretty_name: 'Save data',
        default: false,
        description: 'Save the reaction time data for each page or not'
      },
      background: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Background',
        default: 'stars',
        description: 'Whether the background is black or stars.'
      },
      stimuli: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Stimuli',
        array: true,
        default: '',
        description: 'Extra stimuli (array refers to page number) that is trial-dependent.'
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

    const tmp_parameters = loadParameters();
    if (trial.meg_mode){
      trial.key_forward = tmp_parameters.key_responses.left_right_buttons[1];
      trial.key_backward = tmp_parameters.key_responses.left_right_buttons[0];
    }

    // set up background
    document.getElementById('jspsych-content').classList = 'jspsych-content';

    if (trial.background == 'stars') {
      document.getElementsByClassName("stars")[0].style.opacity = "1";
      document.getElementsByClassName("twinkling")[0].style.opacity = "1";
    } else {
      document.getElementsByClassName("stars")[0].style.opacity = "0";
      document.getElementsByClassName("twinkling")[0].style.opacity = "0";
    }

    var current_page = 0;
    var view_history = [];
    var start_time = performance.now();
    var last_page_update_time = start_time;

    function btnListener(evt) {
      evt.target.removeEventListener('click', btnListener);
      if (this.id === "jspsych-instructions-back") {
        back();
      } else if (this.id === 'jspsych-instructions-next') {
        next();
      }
    }

    function show_current_page() {

      var vis_html = '';
      if (trial.map.length != 0) {
        if (trial.map[current_page] == ''){
          vis_html = '';
        } else {
          if (trial.map[current_page].length != 2){
            vis_html = '<div id="map"><img src="assets/img/map/map_' + trial.map[current_page] + '.svg"></div>';
          } else {
            vis_html = '<div id="map-blinking2"><img src="assets/img/map/map_' + trial.map[current_page][0] + '.svg"></div><div id="map-blinking1"><img src="assets/img/map/map_' + trial.map[current_page][1] + '.svg"></div>';
          }
        }
      }

      var html = '';
      if (!trial.stimuli[current_page]){
        html = '<div class="instructions-background">' + trial.pages[current_page] + '</div>';
      } else {
        html = '<div class="instructions-background">' + trial.stimuli[current_page] + trial.pages[current_page] + '</div>';
      }

      if (current_page == 0 & trial.fade_from_black) {
        html = '<div class="start-fade" style="animation-name: fullScreenFadeOut; animation-duration: 1s; background:black;"></div>' + html;
      }

      var pagenum_display = "";
      if (trial.show_page_number) {
        pagenum_display = "<span style='margin: 0 1em;' class='" +
          "jspsych-instructions-pagenum'>Page " + (current_page + 1) + "/" + trial.pages.length + "</span>";
      }

      var nav_html = "<div id='instructions-nav' class='jspsych-instructions-nav' style='padding: 10px 0px;'>";
      if (trial.allow_backward) {
        var allowed = (current_page > 0) ? '' : "disabled='disabled'";
        nav_html += "<button id='jspsych-instructions-back' class='instructions-btn' style='margin-right: 5px;' " + allowed + ">&lt; " + trial.button_label_previous + "</button>";
      }
      if (trial.pages.length > 1 && trial.show_page_number) {
        nav_html += pagenum_display;
      }
      if (trial.meg_mode && current_page == trial.pages.length-1 && !trial.allow_meg)
      {
        nav_html += "</div>";
        trial.key_forward = "enter";
        if (trial.pages.length > 1)
        {
          jsPsych.pluginAPI.cancelKeyboardResponse(keyboard_listener);
          keyboard_listener = jsPsych.pluginAPI.getKeyboardResponse({
            callback_function: after_response,
            valid_responses: [trial.key_forward, trial.key_backward],
            rt_method: 'performance',
            persist: false
          });
        }
      } else {
        nav_html += "<button id='jspsych-instructions-next' class='instructions-btn'" +
          "style='margin-left: 5px;'>" + trial.button_label_next +
          " &gt;</button></div>";
          if (trial.meg_mode && !trial.allow_meg){
            trial.key_forward = tmp_parameters.key_responses.left_right_buttons[1];
            trial.key_backward = tmp_parameters.key_responses.left_right_buttons[0];
          }
      }

      var cwidth = '';
      if (trial.meg_mode){
        cwidth = 'style="width: 75vw; max-height:90vh;"';
      }

      html += nav_html;
      display_element.innerHTML = '<div id="instructions-container"' + cwidth + '>' + html + '</div>' + vis_html;

      if (current_page != 0 && trial.allow_backward) {
        display_element.querySelector('#jspsych-instructions-back').addEventListener('click', btnListener);
      }

      if (!(trial.meg_mode && current_page == trial.pages.length-1 && !trial.allow_meg))
      {
        display_element.querySelector('#jspsych-instructions-next').addEventListener('click', btnListener);
      }

      if (current_page == 0 & trial.allow_backward) {
        display_element.querySelector('#jspsych-instructions-back').style.display = "none";
      } else if (current_page == trial.pages.length) {
        display_element.querySelector('#jspsych-instructions-next').innerText = "Begin";
      }

    }

    function next() {

      add_current_page_to_view_history();

      current_page++;

      // if done, finish up...
      if (current_page >= trial.pages.length) {

        endTrial();

      } else {
        show_current_page();
      }

    }

    function back() {

      add_current_page_to_view_history();

      current_page--;

      show_current_page();
    }

    function add_current_page_to_view_history() {

      var current_time = performance.now();

      var page_view_time = current_time - last_page_update_time;

      view_history.push({
        page_index: current_page,
        viewing_time: page_view_time
      });

      last_page_update_time = current_time;
    }

    function endTrial() {

      if (trial.allow_keys) {
        jsPsych.pluginAPI.cancelKeyboardResponse(keyboard_listener);
      }

      display_element.innerHTML += '';

      var trial_data = {
        "view_history": JSON.stringify(view_history),
        "rt": performance.now() - start_time
      };

      if (trial.save_data){
        jsPsych.finishTrial(trial_data);
      } else {
        jsPsych.finishTrial();
      }
    }

    var after_response = function(info) {

      // have to reinitialize this instead of letting it persist to prevent accidental skips of pages by holding down keys too long
      keyboard_listener = jsPsych.pluginAPI.getKeyboardResponse({
        callback_function: after_response,
        valid_responses: [trial.key_forward, trial.key_backward],
        rt_method: 'performance',
        persist: false,
        allow_held_key: false
      });
      // check if key is forwards or backwards and update page
      if (jsPsych.pluginAPI.compareKeys(info.key, trial.key_backward)) {
        if (current_page !== 0 && trial.allow_backward) {
          back();
        }
      }

      if (jsPsych.pluginAPI.compareKeys(info.key, trial.key_forward)) {
        next();
      }

    };

    show_current_page();

    if (trial.allow_keys) {
      var keyboard_listener = jsPsych.pluginAPI.getKeyboardResponse({
        callback_function: after_response,
        valid_responses: [trial.key_forward, trial.key_backward],
        rt_method: 'performance',
        persist: false
      });
    }
  };

  return plugin;
})();
