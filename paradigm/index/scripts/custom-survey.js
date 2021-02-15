/**
 * custom-survey
 *
 */
 /*jshint esversion: 6 */


jsPsych.plugins['custom-survey'] = (function() {

  var plugin = {};

  plugin.info = {
    name: 'custom-survey',
    description: '',
    parameters: {
      questions: {
        type: jsPsych.plugins.parameterType.COMPLEX,
        array: true,
        pretty_name: 'Questions',
        default: undefined,
        nested: {
          prompt: {
            type: jsPsych.plugins.parameterType.STRING,
            pretty_name: 'Prompt',
            default: undefined,
            description: 'Prompt for the subject to response'
          },
          placeholder: {
            type: jsPsych.plugins.parameterType.STRING,
            pretty_name: 'Value',
            default: "",
            description: 'Placeholder text in the textfield.'
          },
          prefilled: {
            type: jsPsych.plugins.parameterType.STRING,
            pretty_name: 'Prefilled',
            default: null,
            description: 'What to prefill the field with.'
          },
          rows: {
            type: jsPsych.plugins.parameterType.INT,
            pretty_name: 'Rows',
            default: 1,
            description: 'The number of rows for the response text box.'
          },
          columns: {
            type: jsPsych.plugins.parameterType.INT,
            pretty_name: 'Columns',
            default: 40,
            description: 'The number of columns for the response text box.'
          },
          required: {
            type: jsPsych.plugins.parameterType.BOOL,
            pretty_name: 'Required',
            default: false,
            description: 'Require a response'
          },
          name: {
            type: jsPsych.plugins.parameterType.STRING,
            pretty_name: 'Question Name',
            default: '',
            description: 'Controls the name of data values associated with this question'
          }
        }
      },
      preamble: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Preamble',
        default: null,
        description: 'HTML formatted string to display at the top of the page above all the questions.'
      },
      button_label: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Button label',
        default:  'Continue',
        description: 'The text that appears on the button to finish the trial.'
      },
      ing: {
        type: jsPsych.plugins.parameterType.BOOL,
        pretty_name: 'Reverse scoring',
        array: true,
        default: null,
        description: 'Which of the questions are reverse scored'
      },
      categories: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Categories',
        array: true,
        default: null,
        description: 'Which category each question belongs to'
      },
      scale: {
        type: jsPsych.plugins.parameterType.STRING,
        pretty_name: 'Scale',
        array: true,
        default: undefined,
        description: 'The words used to describe each end of the scale'
      },
      num_buttons: {
        type: jsPsych.plugins.parameterType.INT,
        pretty_name: 'Number of buttons',
        array: false,
        default: 3,
        description: 'How many radio buttons to display'
      }
    }
  };

  plugin.trial = function(display_element, trial) {

    // make background white
    document.getElementsByClassName("stars")[0].style.opacity = "1";
    document.getElementsByClassName("twinkling")[0].style.opacity = "1";
    document.getElementById('jspsych-content').classList = 'jspsych-content';

    trial.reverse_scoring = trial.reverse_scoring == undefined ? Array(trial.questions.length).fill(false) : trial.reverse_scoring;

    for (let i = 0; i < trial.questions.length; i++) {
      if (typeof trial.questions[i].rows == 'undefined') {
        trial.questions[i].rows = 1;
      }
    }
    for (let i = 0; i < trial.questions.length; i++) {
      if (typeof trial.questions[i].columns == 'undefined') {
        trial.questions[i].columns = 40;
      }
    }
    for (let i = 0; i < trial.questions.length; i++) {
      if (typeof trial.questions[i].value == 'undefined') {
        trial.questions[i].value = "";
      }
    }

    var html = '';
    // show preamble text
    if(trial.preamble !== null){
      html += '<div id="preamble"><p>'+trial.preamble+'</p></div>';
    }

    if (trial.scale.length != trial.num_buttons){
      // start form
      html += '<form id="main-form" style="background:hsl(0,0%,0%); display:grid; grid-template-columns:50% 50%; padding:10px;">';
      // headings
      html += '<p></p><div id="form-headings" style="display:flex; flex-direction:row; justify-content:space-between;"><p>' + trial.scale[0] + '</p><p>' + trial.scale[1] + '</p></div>';
    } else {
      var grid_text = '35% ';
      var scale_text = '<p></p>';
      for (let i = 0; i < trial.scale.length; i++){
        grid_text += (100-35)/trial.scale.length + '% ';
        scale_text += '<p>' + trial.scale[i] + '</p>';
      }
      // start form
      html += '<form id="main-form" style="background:hsl(0,0%,0%); display:grid; grid-template-columns:' + grid_text + '; padding:10px;">';
      // headings
      html += scale_text;
    }

    // generate question order
    var question_order = [];
    for(var i=0; i<trial.questions.length; i++){
      question_order.push(i);
    }
    if(trial.randomize_question_order){
      question_order = jsPsych.randomization.shuffle(question_order);
    }

    // add questions
    for (var i = 0; i < trial.questions.length; i++) {
      var bg_color = i % 2 ? '' : 'background:hsl(0,0%,15%);';
      var question = trial.questions[question_order[i]];
      var question_index = question_order[i];
      html += '<p id="q'+i+'-text" style="' + bg_color + ' margin:0; padding:10px;">' + question.prompt + '</p>';
      var autofocus = i == 0 ? "autofocus" : "";
      var req = question.required ? "required" : "";
      var value_text = '';
      if (question.prefilled != null){
        value_text = 'value="' + question.prefilled + '"';
      }

      var answer_html = '';
      for (let r = 0; r < trial.num_buttons; r++){
        if (trial.scale.length != trial.num_buttons){
          answer_html += '<label class="custom-label"';
        } else {
          answer_html += '<label class="custom-label" style="' + bg_color + '"';
        }
        answer_html += ' for="q'+i+'-r'+r+'" style="padding-right:5px;">' + (r+1) + '<input type="radio" id="q'+i+'-r'+r+'" value="' + (r+1) + '" name="q'+i+'"><span class="custom-radio"></span></label>';
      }
      if (trial.scale.length != trial.num_buttons){
        html += '<div style="display:flex; flex-direction:row; justify-content:space-around; align-items:center; ' + bg_color + '">' + answer_html + '</div>';
      } else {
        html += answer_html;
      }
    }

    // add submit button
    var submit_html = '<button id="submit-btn" class="instructions-btn" style="margin-top: 15px;">' + trial.button_label + '</button>';

    html += '</form>';
    display_element.innerHTML = '<div id="instructions-container" style="width:90vw;"><div class="instructions-background">' + html + '</div>' + submit_html + '</div>';

    // backup in case autofocus doesn't work
    document.getElementById('q0-r0').focus();

    document.getElementById('submit-btn').addEventListener('click', function(e) {

      e.preventDefault();

      // check for missed responses
      var response_checklist = [];
      for (let i = 0; i < trial.questions.length; i++){
        var these_responses = [];
        for (let r = 0; r < trial.num_buttons; r++){
          these_responses.push(document.getElementById("q"+i+"-r"+r).checked);
        }
        response_checklist.push(these_responses.filter(x => x).length);
        if (response_checklist.slice(-1) == 0){
          document.getElementById("q"+i+"-text").style.color = 'rgb(255,0,0)';
        } else {
          document.getElementById("q"+i+"-text").style.color = '';
        }
      }

      if (response_checklist.filter(x => x == 1).length == trial.questions.length){
        after_response();
      } else {
        alert("Please select an answer for every question (missed questions are highlighted in red).");
      }
    });

    function after_response() {

      // measure response time
      var response_time = performance.now() - startTime;

      // save answers
      var answer_list = [];
      var reverse_idx = [5,4,3,2,1];
      for (let i = 0; i < trial.questions.length; i++){
        for (let r = 0; r < trial.num_buttons; r++){
          if (document.getElementById("q"+i+"-r"+r).checked){
            let this_answer = parseInt(document.getElementById("q"+i+"-r"+r).value);
            if (trial.reverse_scoring[i]){
              answer_list.push(reverse_idx[this_answer-1]);
            } else {
              answer_list.push(this_answer);
            }
          }
        }
      }

      if (trial.categories != null){
        var unique_categories = trial.categories.filter(function(item, pos){
          return trial.categories.indexOf(item)== pos;
        });
        var category_scores = Array(unique_categories.length).fill(0);
        for (let i = 0; i < trial.questions.length; i++){
          category_scores[unique_categories.indexOf(trial.categories[i])] += answer_list[i];
        }
      }

      // calculate score
      var total_score = answer_list.reduce((a,b) => a + b,0);

      // save data
      var trialdata = {
        "rt": response_time,
        "answer_list": answer_list,
        "total_score": total_score
      };

      if (trial.categories != null){
        trialdata.categories = trial.categories;
        for (let i = 0; i < category_scores.length; i++){
          trialdata[unique_categories[i]] = category_scores[i];
        }
      }

      display_element.innerHTML = '';

      // next trial
      jsPsych.finishTrial(trialdata);

    }

    var startTime = performance.now();
  };

  return plugin;
})();
