/*jshint esversion: 6 */

window.loadTimeline = function(parameters, db, uid) {
  ////////////
  // SET UP //
  ////////////

  function getAllIndexes(arr, val) {
    var indexes = [],
      i = -1;
    while ((i = arr.indexOf(val, i + 1)) != -1) {
      indexes.push(i);
    }
    return indexes;
  }

  // Load experiment structure
  var exp_list;
  if (parameters.exp_variables.meg_mode == false & parameters.exp_variables.structure_type == "A") {
    exp_list = exp_list_behav_A;
  } else if (parameters.exp_variables.meg_mode == false & parameters.exp_variables.structure_type == "B") {
    exp_list = exp_list_behav_B;
  } else if (parameters.exp_variables.meg_mode == true & parameters.exp_variables.structure_type == "A") {
    exp_list = exp_list_MEG_A;
  } else if (parameters.exp_variables.meg_mode == true & parameters.exp_variables.structure_type == "B") {
    exp_list = exp_list_MEG_B;
  }

  var exp_structure = JSON.parse(exp_list)[0];
  for (let key in exp_structure) {
    exp_structure[key] = Object.values(exp_structure[key]);
  }

  // Set up state images
  var state_names = jsPsych.randomization.shuffle(
    parameters.state_images.all_states
  );

  var seq_array = [state_names.slice(0, 3), state_names.slice(3, 6)];
  var seq_array_filenames = [
    seq_array[0].map(
      x => "assets/img/states/" + x + parameters.state_images.file_ext
    ),
    seq_array[1].map(
      x => "assets/img/states/" + x + parameters.state_images.file_ext
    )
  ];

  var image_info = {
    state_names: seq_array.flat(1),
    seq_array: seq_array,
    seq_array_filenames: seq_array_filenames
  };

  var background = "stars"; // 'black' or 'stars'
  if (parameters.exp_variables.meg_mode) {
    background = "black"; // 'black' or 'stars'
  }

  ////////////////////
  // QUESTIONNAIRES //
  ////////////////////

  var questionnaire_timeline = [];
  if (parameters.exp_variables.questionnaires) {
    var q_instructions = {
      type: "custom-instructions",
      pages: [
        "<h1>Questionnaires</h1><p>Before we start the game, there are a series of short questionnaires for you to fill out.</p><p>Please answer as honestly and accurately as you can. All your answers will remain anonymous.</p>"
      ],
      allow_backward: false,
      button_label_next: "Begin",
      fade_from_black: true,
      background: background,
      meg_mode: parameters.exp_variables.meg_mode
    };
    questionnaire_timeline.push(q_instructions);

    // Intolerance of uncertainty
    var q_ius = {
      type: "custom-survey",
      preamble: "Please select the number that best corresponds to how much you agree with each.",
      scale: [
        "Not at all characteristic of me",
        "A little characteristic of me",
        "Somewhat characteristic of me",
        "Very characteristic of me",
        "Entirely characteristic of me"
      ],
      num_buttons: 5,
      questions: [{
          prompt: "Unforseen events upset me greatly.",
          required: true
        },
        {
          prompt: "It frustrates me not having all the information I need.",
          required: true
        },
        {
          prompt: "Uncertainty keeps me from living a full life.",
          required: true
        },
        {
          prompt: "One should always look ahead so as to avoid surprises.",
          required: true
        },
        {
          prompt: "A small unforseen event can spoil everything, even with the best of planning.",
          required: true
        },
        {
          prompt: "When it's time to act, uncertainty paralyses me.",
          required: true
        },
        {
          prompt: "When I am uncertain, I can't function very well",
          required: true
        },
        {
          prompt: "I always want to know what the future has in store for me.",
          required: true
        },
        {
          prompt: "I can't stand being taken by surprise.",
          required: true
        },
        {
          prompt: "The smallest doubt can stop me from acting.",
          required: true
        },
        {
          prompt: "I should be able to organise everything in advance.",
          required: true
        },
        {
          prompt: "I must get away from all uncertain situations.",
          required: true
        }
      ],
      on_finish: function() {
        let filtered_data = JSON.parse(
          jsPsych.data
          .get()
          .filter({
            trial_type: "custom-survey"
          })
          .json()
        ).slice(-1)[0];
        let trial_data = {
          trial_type: "questionnaire-ius",
          time_elapsed: filtered_data.time_elapsed,
          rt: filtered_data.rt,
          answer_list: filtered_data.answer_list,
          total_score: filtered_data.total_score
        };

        db.collection("iterations")
          .doc("pilot_v2.5")
          .collection("subjects")
          .doc(uid)
          .collection("questionnaires")
          .doc("intolerance_of_uncertainty")
          .set({
            trial_data
          });
      }
    };
    questionnaire_timeline.push(q_ius);

    // Worry questionnaire
    var q_worry = {
      type: "custom-survey",
      preamble: "Rate each of the following statements on a scale of 1 (“not at all typical of me”) to 5 (“very typical of me”). Please do not leave any items blank.",
      scale: ["Not at all typical of me", "Very typical of me"],
      num_buttons: 5,
      questions: [{
          prompt: "If I do not have enough time to do everything, I do not worry about it.",
          required: true
        },
        {
          prompt: "My worries overwhelm me.",
          required: true
        },
        {
          prompt: "I do not tend to worry about things.",
          required: true
        },
        {
          prompt: "Many situations make me worry.",
          required: true
        },
        {
          prompt: "I know I should not worry about things, but I just cannot help it.",
          required: true
        },
        {
          prompt: "When I am under pressure I worry a lot.",
          required: true
        },
        {
          prompt: "I am always worrying about something.",
          required: true
        },
        {
          prompt: "I find it easy to dismiss worrisome thoughts.",
          required: true
        },
        {
          prompt: "As soon as I finish one task, I start to worry about everything else I have to do.",
          required: true
        },
        {
          prompt: "I never worry about anything.",
          required: true
        },
        {
          prompt: "When there is nothing more I can do about a concern, I do not worry about it any more.",
          required: true
        },
        {
          prompt: "I have been a worrier all my life.",
          required: true
        },
        {
          prompt: "I notice that I have been worrying about things.",
          required: true
        },
        {
          prompt: "Once I start worrying, I cannot stop.",
          required: true
        },
        {
          prompt: "I worry all the time.",
          required: true
        },
        {
          prompt: "I worry about projects until they are all done.",
          required: true
        }
      ],
      button_label: "Next",
      reverse_scoring: [
        true,
        false,
        true,
        false,
        false,
        false,
        false,
        true,
        false,
        true,
        true,
        false,
        false,
        false,
        false,
        false
      ],
      on_finish: function() {
        let filtered_data = JSON.parse(
          jsPsych.data
          .get()
          .filter({
            trial_type: "custom-survey"
          })
          .json()
        ).slice(-1)[0];
        let trial_data = {
          trial_type: "questionnaire-worry",
          time_elapsed: filtered_data.time_elapsed,
          rt: filtered_data.rt,
          answer_list: filtered_data.answer_list,
          total_score: filtered_data.total_score
        };

        db.collection("iterations")
          .doc("pilot_v2.5")
          .collection("subjects")
          .doc(uid)
          .collection("questionnaires")
          .doc("worry")
          .set({
            trial_data
          });
      }
    };
    questionnaire_timeline.push(q_worry);

    // Domain-sepcific risk
    var q_risk = {
      type: "custom-survey",
      preamble: "For each of the following statements, please indicate the likelihood that you would engage in the described activity or behavior if you were to find yourself in that situation.",
      scale: [
        "Extremely Unlikely",
        "Moderately Unlikely",
        "Somewhat Unlikely",
        "Not Sure",
        "Somewhat Likely",
        "Moderately Likely",
        "Extremely Likely"
      ],
      num_buttons: 7,
      questions: [{
          prompt: "Admitting that your tastes are different from those of a friend.",
          required: true
        },
        {
          prompt: "Going camping in the wilderness.",
          required: true
        },
        {
          prompt: "Betting a day's income at the horse races.",
          required: true
        },
        {
          prompt: "Investing 10% of your annual income in a moderate growth mutual fund.",
          required: true
        },
        {
          prompt: "Drinking heavily at a social function.",
          required: true
        },
        {
          prompt: "Taking some questionable deductions on your income tax return.",
          required: true
        },
        {
          prompt: "Disagreeing with an authority figure on a major issue.",
          required: true
        },
        {
          prompt: "Betting a day's income at a high-stake poker game.",
          required: true
        },
        {
          prompt: "Having an affiar with a married man/woman.",
          required: true
        },
        {
          prompt: "Passing off somebody else's work as your own.",
          required: true
        },
        {
          prompt: "Going down a ski run that is beyond your ability.",
          required: true
        },
        {
          prompt: "Investing 5% of your annual income in a very speculative stock.",
          required: true
        },
        {
          prompt: "Going whitewater rafting at high water in the spring.",
          required: true
        },
        {
          prompt: "Betting a day's income on the outcome of a sporting event.",
          required: true
        },
        {
          prompt: "Engaging in unprotected sex.",
          required: true
        },
        {
          prompt: "Revealing a friend's secret to someone else.",
          required: true
        },
        {
          prompt: "Driving a car without wearing a seat belt.",
          required: true
        },
        {
          prompt: "Investing 10% of your annual income in a new business venture.",
          required: true
        },
        {
          prompt: "Taking a skydiving class.",
          required: true
        },
        {
          prompt: "Riding a motorcycle without a helmet.",
          required: true
        },
        {
          prompt: "Choosing a career that you truly enjoy over a more secure one.",
          required: true
        },
        {
          prompt: "Speaking your mind about an unpopular issue in a meeting at work.",
          required: true
        },
        {
          prompt: "Sunbathing without sunscreen.",
          required: true
        },
        {
          prompt: "Bungee jumping off a tall bridge.",
          required: true
        },
        {
          prompt: "Piloting a small plane.",
          required: true
        },
        {
          prompt: "Walking home alone at night in an unsafe area of town.",
          required: true
        },
        {
          prompt: "Moving to a city far away from your extended family.",
          required: true
        },
        {
          prompt: "Starting a new career in your mid-thirties.",
          required: true
        },
        {
          prompt: "Leaving your young children alone at home while running an errand.",
          required: true
        },
        {
          prompt: "Not returning a wallet you found that contains $200.",
          required: true
        }
      ],
      categories: [
        "social",
        "recreational",
        "financial",
        "financial",
        "healthsafety",
        "ethical",
        "social",
        "financial",
        "ethical",
        "ethical",
        "recreational",
        "financial",
        "recreational",
        "financial",
        "healthsafety",
        "ethical",
        "healthsafety",
        "financial",
        "recreational",
        "healthsafety",
        "social",
        "social",
        "healthsafety",
        "recreational",
        "recreational",
        "healthsafety",
        "social",
        "social",
        "ethical",
        "ethical"
      ],
      on_finish: function() {
        let filtered_data = JSON.parse(
          jsPsych.data
          .get()
          .filter({
            trial_type: "custom-survey"
          })
          .json()
        ).slice(-1)[0];
        let trial_data = {
          trial_type: "questionnaire-risk",
          time_elapsed: filtered_data.time_elapsed,
          rt: filtered_data.rt,
          answer_list: filtered_data.answer_list,
          total_score: filtered_data.total_score,
          categories: filtered_data.categories,
          social_score: filtered_data.social,
          recreational_score: filtered_data.recreational,
          financial_score: filtered_data.financial,
          healthsafety_score: filtered_data.healthsafety,
          ethical_score: filtered_data.ethical
        };

        db.collection("iterations")
          .doc("pilot_v2.5")
          .collection("subjects")
          .doc(uid)
          .collection("questionnaires")
          .doc("risk")
          .set({
            trial_data
          });
      }
    };
    questionnaire_timeline.push(q_risk);

    // questionnaire_timeline = jsPsych.randomization.shuffle(questionnaire_timeline); // randomly shuffle the order
  }

  //////////////////////////
  // FUNCTIONAL LOCALISER //
  //////////////////////////
  var functional_localiser = [];
  if (parameters.exp_variables.meg_mode) {
    // Parameters
    var word_descriptions = [
      state_names,
      "hat",
      "pen",
      "dog",
      "bed",
      "garden",
      "carrot",
      "phone",
      "bottle",
      "shoe",
      "sofa",
      "scissors"
    ].flat(Infinity);
    var stim_list = seq_array_filenames.flat(1);
    var fl_reps = parameters.trial_numbers.nTrlLocaliser / stim_list.length;

    // Instructions
    var fl_instructions = {
      type: "custom-instructions",
      pages: [
        "<h1>Image Viewing</h1><br><p style='text-align:center;'>You will be shown a series of different images.</p><p style='text-align:center;'>After each image is presented, you will be shown <strong>TWO WORDS</strong>.</p><p style='text-align:center;'>You must pick either the <strong>LEFT</strong> or the <strong>RIGHT</strong> word that best describes the image.</p><p style='text-align:center;'><big>Please wait for the experimenter to start the program.<br><br></big><strong>Remember to stay very still.</strong></p>"
      ],
      allow_backward: false,
      button_label_next: "Begin",
      meg_mode: true,
      fade_from_black: true,
      background: "black" // 'black' or 'stars'
    };
    functional_localiser.push(fl_instructions);

    // Trials
    for (
      var block = 0; block < parameters.trial_numbers.nBlockLocaliser; block++
    ) {
      // get stimulus order for this block
      var fl_stimuli = [];
      for (let i = 0; i < fl_reps; i++) {
        fl_stimuli.push(jsPsych.randomization.shuffle(stim_list));
      }
      fl_stimuli = fl_stimuli.flat(Infinity);

      // loop through trials for this block
      for (var trl = 0; trl < parameters.trial_numbers.nTrlLocaliser; trl++) {
        var this_stim_fname = fl_stimuli.flat(Infinity)[trl];
        var this_state = "";
        var this_path = "";
        var this_stim = "";
        for (let i = 0; i < 2; i++) {
          let idx = seq_array_filenames[i].indexOf(this_stim_fname);
          this_stim =
            (idx != -1) & (this_stim == "") ? seq_array[i][idx] : this_stim;
          this_state = (idx != -1) & (this_state == "") ? idx : this_state;
          this_path = (idx != -1) & (this_path == "") ? i : this_path;
        }
        var these_words = [
          this_stim,
          jsPsych.randomization.shuffle(word_descriptions.filter(x => x != this_stim))[0]
        ];
        these_words = jsPsych.randomization.shuffle(these_words);

        var fl_trial = {
          type: "localiser",
          stimuli: this_stim,
          stimuli_fname: this_stim_fname,
          path: this_path,
          stim_num: this_state,
          img_duration: parameters.timing.localiser_img,
          isi: parameters.timing.localiser_isi,
          trial_num: trl,
          block_num: block,
          trial_end_type: "min-time",
          choices: parameters.key_responses.left_right_buttons.flat(Infinity),
          words: these_words,
          on_finish: function() {
            let filtered_data = JSON.parse(
              jsPsych.data
              .get()
              .filter({
                trial_type: "localiser"
              })
              .json()
            ).slice(-1)[0];
            let trial_data = {
              trial_type: "localiser",
              trial: filtered_data.trial + 1,
              block: filtered_data.block + 1,
              time_elapsed: filtered_data.time_elapsed,
              img: filtered_data.img,
              state: filtered_data.state,
              path: filtered_data.path,
              words: filtered_data.words,
              rt: filtered_data.rt,
              isi: filtered_data.isi,
              acc: filtered_data.acc,
              triggers: filtered_data.triggers
            };

            db.collection("iterations")
              .doc("pilot_v2.5")
              .collection("subjects")
              .doc(uid)
              .collection("0_functional_localiser")
              .doc(
                "block" +
                padNumber(trial_data.block) +
                "_trial" +
                padNumber(trial_data.trial) +
                "_choice"
              )
              .set({
                trial_data
              });

            // save overall localiser behavioural accuracy to subject's main document
            if (
              trl == parameters.trial_numbers.nTrlLocaliser - 1 &&
              block == parameters.trial_numbers.nBlockLocaliser - 1
            ) {
              var localiser_acc = jsPsych.data
                .get()
                .filter({
                  trial_type: "localiser"
                })
                .select("acc")
                .mean();
              var mainDB = db
                .collection("iterations")
                .doc("pilot_v2.5")
                .collection("subjects")
                .doc(uid);

              var update_mainDB;
              update_mainDB = mainDB.set({
                functional_localiser: localiser_acc
              }, {
                merge: true
              });
            }
          }
        };
        functional_localiser.push(fl_trial);
      }

      // show end of block screen
      var fl_end_block = {};
      if (block != parameters.trial_numbers.nBlockLocaliser - 1) {
        fl_end_block = {
          type: "custom-instructions",
          pages: [
            "<p style='text-align:center;'>End of block " +
            (block + 1) +
            " (of " +
            parameters.trial_numbers.nBlockLocaliser +
            ")</p><p>Please wait for the experimenter.</p>"
          ],
          allow_backward: false,
          meg_mode: true,
          button_label_next: "Continue",
          fade_from_black: true,
          background: "black" // 'black' or 'stars'
        };
      } else {
        fl_end_block = {
          type: "custom-instructions",
          pages: [
            "<p style='text-align:center;'>End of final block</p><p>Please wait for the experimenter.</p>"
          ],
          allow_backward: false,
          meg_mode: true,
          button_label_next: "Continue",
          fade_from_black: true,
          background: "black" // 'black' or 'stars'
        };
      }
      functional_localiser.push(fl_end_block); // push trials for this block
    }
  }

  /////////////////////////////
  // TUTORIAL 1: EXPLORATION //
  /////////////////////////////

  // Parameters
  var exploration_timeline = [];
  var trl_counter = -1;
  var choice_list = [0, 1]; // path 1, path 2

  var instructions_pages = [];
  if (parameters.exp_variables.meg_mode) {
    instructions_pages = [
      "<h1>Exploration</h1><p>You will now be shown the different images that represent each room in the spaceship. Memorise the images and the order they are shown in for each door.</p><p style='text-align:center;'><strong>Ready?</strong></p><p style='text-align:center;'>Please wait for the experimenter to start the program.</p><p style='text-align:center;'><strong>Remember to stay very still during the scan.</strong></p>"
    ];
  } else {
    instructions_pages = [
      "<h1>Exploration</h1>" +
      "<p>Imagine that you are on a spaceship. You have just entered an <strong>AIRLOCK</strong>. In front of you are now <strong>2 DOORS</strong> that each lead to a passageway containing 3 rooms.</p><p>Your task is to memorise the order of rooms behind each door. Each room will be represented by a unique <strong>image</strong>.</p>",

      "<p style='text-align:center; color:#FFAE57'><strong>YOUR TASK:</strong> Memorise the order of images behind DOOR 1 and DOOR 2.</p><p style='text-align:center;'><strong>Ready?</strong></p><p style='text-align:center;'>Press 'Next' to begin</p><p style='text-align:center;'><small>(or press 'Previous' to look back at the instructions)</small></p><p style='text-align:center;'>(Note: You can either CLICK the buttons below using your mouse, or you can use the ARROW KEYS on your keyboard.)</p>"
    ];
  }

  // Instructions
  var instructions = {
    type: "custom-instructions",
    pages: instructions_pages,
    fade_from_black: true,
    background: background,
    meg_mode: parameters.exp_variables.meg_mode,
    map: []
  };

  // Loop through trials
  // (trial number = no. of paths * repeat count in parameters.trial_numbers.nRepsT1)
  for (let rep = 0; rep < parameters.trial_numbers.nRepsT1; rep++) {
    for (let trl = 0; trl < choice_list.length; trl++) {
      trl_counter++;
      var trl_stimuli = seq_array[choice_list[trl]];

      var exploration_choice = {
        type: "exploration-choice",
        choice: choice_list[trl],
        trial_num: trl_counter,
        choices: parameters.key_responses.left_right_buttons[choice_list[trl]], // either LEFT or RIGHT, depending on whether Door 1 or Door 2 is selected for this forced-choice trial
        background: background,
        on_finish: function() {
          let filtered_data = JSON.parse(
            jsPsych.data
            .get()
            .filter({
              trial_type: "exploration-choice"
            })
            .json()
          ).slice(-1)[0];
          let trial_data = {
            trial_type: "exploration-choice",
            trial: filtered_data.trial + 1,
            block: filtered_data.attempt_num + 1,
            time_elapsed: filtered_data.time_elapsed,
            rt: filtered_data.rt,
            choice: filtered_data.door_name,
            button: filtered_data.button
          };

          db.collection("iterations")
            .doc("pilot_v2.5")
            .collection("subjects")
            .doc(uid)
            .collection("1_exploration_data")
            .doc(
              "block" +
              padNumber(trial_data.block) +
              "_trial" +
              padNumber(trial_data.trial) +
              "_choice"
            )
            .set({
              trial_data
            });
        }
      };
      exploration_timeline.push(exploration_choice);

      for (let stim = 0; stim < 4; stim++) {
        var this_stim = stim < 3 ? trl_stimuli[stim] : null;
        var this_stim_filename =
          stim < 3 ? seq_array_filenames[choice_list[trl]][stim] : null;

        var exploration_animation = {
          type: "exploration-animation",
          stimuli: [this_stim, this_stim_filename],
          stim_num: stim,
          trial_num: trl_counter,
          stimuli_info: image_info,
          seq_num: choice_list[trl],
          isi: parameters.timing.exploration_isi,
          trial_duration: parameters.timing.exploration_img_duration,
          trial_end_type: parameters.timing.exploration_response_type,
          background: background,
          choices: parameters.key_responses.left_right_buttons[1], // only RIGHT button
          on_finish: function() {
            let filtered_data = JSON.parse(
              jsPsych.data
              .get()
              .filter({
                trial_type: "exploration-animation"
              })
              .json()
            ).slice(-1)[0];
            let trial_data = {
              trial_type: "exploration-animation",
              trial: filtered_data.trial + 1,
              block: filtered_data.attempt_num + 1,
              time_elapsed: filtered_data.time_elapsed,
              img: filtered_data.img,
              state: filtered_data.state,
              path: filtered_data.path,
              isi: filtered_data.isi
            };

            if (stim < 3) {
              db.collection("iterations")
                .doc("pilot_v2.5")
                .collection("subjects")
                .doc(uid)
                .collection("1_exploration_data")
                .doc(
                  "block" +
                  padNumber(trial_data.block) +
                  "_trial" +
                  padNumber(trial_data.trial) +
                  "_stim" +
                  stim
                )
                .set({
                  trial_data
                });
            }
          }
        };
        exploration_timeline.push(exploration_animation);
      }
    }
  }

  ////////////////////////////////////////
  // TEST 1: EXPLORATION MEMORY TEST ////
  ///////////////////////////////////////

  var memory_test_cue = {
    type: "custom-instructions",
    pages: [
      "<h1>Memory Test</h1><p>You decide to test your memory for the order of the rooms behind each door, using the images you associated with each room.</p>"
    ],
    allow_backward: false,
    button_label_next: "Begin",
    fade_from_black: true,
    background: background,
    meg_mode: parameters.exp_variables.meg_mode
  };
  exploration_timeline.push(memory_test_cue);

  // decide which path to probe (equal number of times) - first col = path, second col = state from that path to probe
  var test_seq = [];
  for (let i = 0; i < parameters.trial_numbers.memory_test_reps; i++) {
    test_seq.push([0, i]);
    test_seq.push([1, i]);
  }
  test_seq = jsPsych.randomization.shuffle(test_seq);

  for (let seq = 0; seq < test_seq.length; seq++) {
    var memory_test_seq = {
      type: "exploration-test",
      stimuli_info: image_info,
      test_seq: test_seq[seq],
      trial_num: seq,
      last_test: seq == test_seq.length - 1,
      background: background,
      choices: [
        parameters.key_responses.up_down_buttons,
        [
          parameters.key_responses.left_right_buttons[0],
          parameters.key_responses.up_down_buttons,
          parameters.key_responses.left_right_buttons[1]
        ]
      ],
      on_finish: function() {
        let filtered_data = JSON.parse(
          jsPsych.data
          .get()
          .filter({
            trial_type: "exploration-test"
          })
          .json()
        ).slice(-1)[0];
        let trial_data = {
          trial_type: "exploration-test",
          trial: filtered_data.trial + 1,
          block: filtered_data.attempt_num + 1,
          time_elapsed: filtered_data.time_elapsed,
          q1_probe_imgs: filtered_data.first_prompts.map(x =>
            x.split("stim_")[1].slice(2)
          ),
          q1_rt: filtered_data.first_rt,
          q1_choice: filtered_data.first_state,
          q1_button: filtered_data.first_btn,
          q1_acc: filtered_data.first_acc,
          q2_probe_imgs: filtered_data.second_prompts.map(x =>
            x.split("stim_")[1].slice(2)
          ),
          q2_rt: filtered_data.second_rt,
          q2_choice: filtered_data.second_state,
          q2_button: filtered_data.second_btn,
          q2_acc: filtered_data.second_acc
        };

        db.collection("iterations")
          .doc("pilot_v2.5")
          .collection("subjects")
          .doc(uid)
          .collection("2_exploration_test")
          .doc(
            "block" +
            padNumber(trial_data.block) +
            "_trial" +
            padNumber(trial_data.trial)
          )
          .set({
            trial_data
          });

        // save overall test accuracy to subject's main document
        if (trial_data.trial == parameters.trial_numbers.memory_test_reps * 2) {
          var exploration_test_acc = jsPsych.data
            .get()
            .filter({
              trial_type: "exploration-test",
              attempt_num: filtered_data.attempt_num
            })
            .select("rep_acc")
            .values.slice(-1)[0];
          var mainDB = db
            .collection("iterations")
            .doc("pilot_v2.5")
            .collection("subjects")
            .doc(uid);

          var update_mainDB;
          if (trial_data.block == 1) {
            update_mainDB = mainDB.set({
              exploration_test_1: exploration_test_acc
            }, {
              merge: true
            });
          } else if (trial_data.block == 2) {
            update_mainDB = mainDB.set({
              exploration_test_2: exploration_test_acc
            }, {
              merge: true
            });
          } else if (trial_data.block == 3) {
            update_mainDB = mainDB.set({
              exploration_test_3: exploration_test_acc
            }, {
              merge: true
            });
          }
        }
      }
    };
    exploration_timeline.push(memory_test_seq);
  }

  // /* end message & save data */
  // var end_message = {
  //   type: "end-experiment",
  //   background: "black",
  //   on_finish: function() {
  //     window.location.href =
  //       "https://app.prolific.co/submissions/complete?cc=RL53GSEY";
  //     jsPsych.endExperiment("Experiment complete");
  //   }
  // };
  var end_message = {
    type: "custom-instructions",
    background: background,
    pages: ['<h1>Finished!</h1><p>Thank you for completing the experiment! Please inform the experimenter.</p>'],
    allow_keys: false,
    meg_mode: parameters.exp_variables.meg_mode
  };

  /* check if they failed the learning phase, and repeat if the performance was below threshold AND the try count is still low enough */
  var repeat_exploration = {
    timeline: exploration_timeline,
    conditional_function: function() {
      let below_thresh =
        jsPsych.data
        .getLastTimelineData()
        .filter({
          trial_type: "exploration-test"
        })
        .select("rep_acc")
        .values.filter(x => x != null) <
        parameters.thresholds.learning_threshold;
      let too_many_tries =
        jsPsych.data
        .getLastTimelineData()
        .filter({
          trial_type: "exploration-test"
        })
        .select("attempt_num").values[0] >=
        parameters.thresholds.max_tries - 1;
      if (below_thresh && !too_many_tries) {
        return true; // EXECUTE timeline in this object
      } else {
        return false; // SKIP timeline in this object
      }
    }
  };

  /* check if they failed the learning phase, and abort if the performance was below threshold AND the try count was too high */
  var abort_exploration = {
    timeline: [end_message],
    conditional_function: function() {
      let below_thresh =
        jsPsych.data
        .getLastTimelineData()
        .filter({
          trial_type: "exploration-test"
        })
        .select("rep_acc")
        .values.filter(x => x != null) <
        parameters.thresholds.learning_threshold;
      let too_many_tries =
        jsPsych.data
        .getLastTimelineData()
        .filter({
          trial_type: "exploration-test"
        })
        .select("attempt_num").values[0] >=
        parameters.thresholds.max_tries - 1;
      if (below_thresh && too_many_tries) {
        return true; // EXECUTE timeline in this object
      } else {
        return false; // SKIP timeline in this object
      }
    }
  };

  ////////////////////////////////
  // TUTORIAL 2: VALUE LEARNING //
  ////////////////////////////////

  var value_learning_instructions_pages = [];
  value_learning_instructions_pages = [
    "<h1>Points</h1>" +
    "<p>Now that you have memorised the rooms behind Door 1 and Door 2 in the Airlock, you can now learn how many <strong>POINTS</strong> you can <span style='color:rgb(0,255,0);'><strong>GAIN</strong></span> or <span style='color:rgb(255,20,20);'><strong>LOSE</strong></span> in each room.</p><p>You will also be shown a third possibility called the <strong>'SUPPLY ROOM'</strong>. This is a room that you can always access, and you will always receive +1 point</p>",
    "<p style='text-align:center; color:#FFAE57'><strong>YOUR TASK:</strong> Memorise the number of points gained or lost in each room.</p><p style='text-align:center;'><strong>Ready?</strong></p><p style='text-align:center;'>Press 'Next' to begin</p><p style='text-align:center;'><small>(or press 'Previous' to look back at the instructions)</small></p>"
  ];

  // Instructions
  var value_learning_instructions = {
    type: "custom-instructions",
    pages: value_learning_instructions_pages,
    map: [],
    background: background,
    meg_mode: parameters.exp_variables.meg_mode
  };

  // Parameters
  var tutorial_structure = {};
  var idx = getAllIndexes(exp_structure.Practice, 1);
  for (let key in exp_structure) {
    tutorial_structure[key] = [];
    for (let i = 0; i < idx.length; i++) {
      tutorial_structure[key].push(exp_structure[key][idx[i]]);
    }
  }
  var state_vals = [
    tutorial_structure.V[0][0],
    tutorial_structure.V[0][1],
    parameters.values.safe_val
  ];

  var value_learning_list = [];
  for (let i = 0; i < parameters.trial_numbers.nRepsT1; i++) {
    // visit the left door, right door, or supply room, in a rotating order
    value_learning_list.push([0, 1, 2]);
  }
  value_learning_list = value_learning_list.flat(Infinity);

  var these_choices = [
    parameters.key_responses.left_right_buttons[0],
    parameters.key_responses.up_down_buttons[0],
    parameters.key_responses.left_right_buttons[1]
  ];

  var value_learning_timeline = [];

  // Trials
  for (let trl = 0; trl < value_learning_list.length; trl++) {
    var value_learning_choice = {
      type: "test-choice",
      choice: value_learning_list[trl],
      meg_mode: parameters.exp_variables.meg_mode,
      is_practice: 1,
      trial_num: trl,
      show_score: false, // whether or not to show oxygen score for block in top-right corner
      background: background,
      choices: these_choices[value_learning_list[trl]],
      on_finish: function() {
        let filtered_data = JSON.parse(
          jsPsych.data
          .get()
          .filter({
            trial_type: "test-choice"
          })
          .json()
        ).slice(-1)[0];
        let trial_data = {
          trial_type: "value-learning-choice",
          trial: filtered_data.trial + 1,
          block: filtered_data.block + 1,
          time_elapsed: filtered_data.time_elapsed,
          option_names: ["Door 1", "Door 2", "Supply Room"],
          choice_type: filtered_data.choice_type == "free" ?
            "free" : filtered_data.choice_type == 0 ?
            "forced_door1" : "forced_door2",
          rt: filtered_data.rt,
          choice: filtered_data.door_name,
          button: filtered_data.button
        };

        db.collection("iterations")
          .doc("pilot_v2.5")
          .collection("subjects")
          .doc(uid)
          .collection("3_value_learning")
          .doc(
            "block" +
            padNumber(trial_data.block) +
            "_trial" +
            padNumber(trial_data.trial) +
            "_choice"
          )
          .set({
            trial_data
          });
      }
    };
    value_learning_timeline.push(value_learning_choice);

    for (let stim = 0; stim < 3; stim++) {
      var value_learning_animation = {
        type: "test-animation",
        stimuli: seq_array[value_learning_list[trl]],
        stimuli_info: image_info,
        background: background,
        meg_mode: parameters.exp_variables.meg_mode,
        is_practice: 1,
        this_val: state_vals,
        stim_num: stim,
        show_score: false, // whether or not to show oxygen score for block in top-right corner
        trial_duration: parameters.timing.exploration_img_duration, // in seconds
        trial_num: trl,
        choices: parameters.key_responses.left_right_buttons[1],
        isi: parameters.timing.exploration_isi,
        trial_end_type: parameters.timing.value_learning_response_type, // "time-out" = ends at trial_duration, "time-or-response" = ends at trial_duration unless pressed earlier, "response" = ends only when button is pressed
        on_finish: function() {
          let filtered_data = JSON.parse(
            jsPsych.data
            .get()
            .filter({
              trial_type: "test-animation"
            })
            .json()
          ).slice(-1)[0];

          if (filtered_data.hasOwnProperty("trial")) {
            let trial_data = {
              trial_type: "value-learning-animation",
              trial: filtered_data.trial + 1,
              block: filtered_data.block + 1,
              time_elapsed: filtered_data.time_elapsed,
              img: filtered_data.img,
              value: filtered_data.outcome,
              state: filtered_data.state,
              path: filtered_data.path,
              isi: filtered_data.isi
            };

            db.collection("iterations")
              .doc("pilot_v2.5")
              .collection("subjects")
              .doc(uid)
              .collection("3_value_learning")
              .doc(
                "block" +
                padNumber(trial_data.block) +
                "_trial" +
                padNumber(trial_data.trial) +
                "_stim" +
                stim
              )
              .set({
                trial_data
              });
          }
        }
      };
      value_learning_timeline.push(value_learning_animation);
    }
  }

  ///////////////////////////////
  // TEST 2: VALUE MEMORY TEST //
  ///////////////////////////////

  var value_test_cue = {
    type: "custom-instructions",
    pages: [
      "<h1>Memory Test</h1><p>You decide to test your memory for how many oxygen points you could find in each room.</p>"
    ],
    allow_backward: false,
    button_label_next: "Begin",
    fade_from_black: true,
    background: background,
    meg_mode: parameters.exp_variables.meg_mode
  };
  value_learning_timeline.push(value_test_cue);

  var val_test_seq = [];
  for (let i = 0; i < parameters.trial_numbers.value_test_reps; i++) {
    val_test_seq.push([0, 1]);
  }
  val_test_seq = jsPsych.randomization.shuffle(val_test_seq.flat(Infinity));
  for (let seq = 0; seq < val_test_seq.length; seq++) {
    var value_test_seq = {
      type: "value-test",
      test_seq: val_test_seq[seq],
      last_test: seq == val_test_seq.length - 1,
      this_val: state_vals,
      this_seq: val_test_seq[seq],
      stimuli_info: image_info,
      meg_mode: parameters.exp_variables.meg_mode,
      choices: [
        parameters.key_responses.left_right_buttons[0],
        parameters.key_responses.up_down_buttons[0],
        parameters.key_responses.left_right_buttons[1]
      ],
      trial_duration: 2, // in seconds
      trial_num: seq,
      on_finish: function() {
        let filtered_data = JSON.parse(
          jsPsych.data
          .get()
          .filter({
            trial_type: "value-test"
          })
          .json()
        ).slice(-1)[0];
        let trial_data = {
          trial_type: "value-test",
          trial: filtered_data.trial + 1,
          block: filtered_data.attempt_num + 1,
          time_elapsed: filtered_data.time_elapsed,
          path: filtered_data.first_seq,
          q1_probe_img: filtered_data.first_img,
          q1_probe_state: filtered_data.first_state,
          q1_probe_values: filtered_data.first_value_probes,
          q1_rt: filtered_data.first_rt,
          q1_choice: filtered_data.first_resp,
          q1_button: filtered_data.first_btn,
          q1_acc: filtered_data.first_acc,
          q2_probe_img: filtered_data.second_img,
          q2_probe_state: filtered_data.second_state,
          q2_probe_values: filtered_data.second_value_probes,
          q2_rt: filtered_data.second_rt,
          q2_choice: filtered_data.second_resp,
          q2_button: filtered_data.second_btn,
          q2_acc: filtered_data.second_acc
        };

        db.collection("iterations")
          .doc("pilot_v2.5")
          .collection("subjects")
          .doc(uid)
          .collection("4_value_test")
          .doc(
            "block" +
            padNumber(trial_data.block) +
            "_trial" +
            padNumber(trial_data.trial)
          )
          .set({
            trial_data
          });

        // save overall test accuracy to subject's main document
        if (trial_data.trial == parameters.trial_numbers.value_test_reps * 2) {
          var value_test_acc = jsPsych.data
            .get()
            .filter({
              trial_type: "value-test",
              attempt_num: filtered_data.attempt_num
            })
            .select("rep_acc")
            .values.slice(-1)[0];
          var mainDB = db
            .collection("iterations")
            .doc("pilot_v2.5")
            .collection("subjects")
            .doc(uid);

          var update_mainDB;
          if (trial_data.block == 1) {
            update_mainDB = mainDB.set({
              value_test_1: value_test_acc
            }, {
              merge: true
            });
          } else if (trial_data.block == 2) {
            update_mainDB = mainDB.set({
              value_test_2: value_test_acc
            }, {
              merge: true
            });
          } else if (trial_data.block == 3) {
            update_mainDB = mainDB.set({
              value_test_3: value_test_acc
            }, {
              merge: true
            });
          }
        }
      }
    };
    value_learning_timeline.push(value_test_seq);
  }

  /* check if they failed the learning phase, and repeat if the performance was below threshold AND the try count is still low enough */
  var repeat_value_learning = {
    timeline: value_learning_timeline,
    conditional_function: function() {
      let below_thresh =
        jsPsych.data
        .getLastTimelineData()
        .filter({
          trial_type: "value-test"
        })
        .select("rep_acc")
        .values.filter(x => x != null) <
        parameters.thresholds.learning_threshold;
      let too_many_tries =
        jsPsych.data
        .getLastTimelineData()
        .filter({
          trial_type: "value-test"
        })
        .select("attempt_num").values[0] >=
        parameters.thresholds.max_tries - 1;
      if (below_thresh && !too_many_tries) {
        return true; // EXECUTE timeline in this object
      } else {
        return false; // SKIP timeline in this object
      }
    }
  };

  /* check if they failed the learning phase, and abort if the performance was below threshold AND the try count was too high */
  var abort_value_learning = {
    timeline: [end_message],
    conditional_function: function() {
      let below_thresh =
        jsPsych.data
        .getLastTimelineData()
        .select("rep_acc")
        .values.filter(x => x != null) <
        parameters.thresholds.learning_threshold;
      let too_many_tries =
        jsPsych.data.getLastTimelineData().select("attempt_num").values[0] >=
        parameters.thresholds.max_tries - 1;
      if (below_thresh && too_many_tries) {
        return true; // EXECUTE timeline in this object
      } else {
        return false; // SKIP timeline in this object
      }
    }
  };

  //////////////////////
  // NEGATOR LEARNING //
  //////////////////////

  var practice_struct = {};
  var idx = getAllIndexes(exp_structure.Practice, 1);
  for (let key in exp_structure) {
    practice_struct[key] = [];
    for (let i = 0; i < idx.length; i++) {
      practice_struct[key].push(exp_structure[key][idx[i]]);
    }
  }

  // check if odd number and apply odd rule
  function checkOdd(arr){
    let thisSum = arr.reduce((a,b) => a + b,0);
    let isOdd = Math.abs(thisSum) % 2;
    let newSum = isOdd ? thisSum*(-1) : thisSum;
    let thisObject = {
      "isOdd": isOdd,
      "oldSum": thisSum,
      "newSum": newSum
    };
    return thisObject;
  }

  // compute sums for example 1 (warning rooms: door 1 room 2, door 2 room 3)
  var sumColours = ['rgb(22,239,243)','rgb(92,156,255)','rgb(195,139,255)'];

  var example1 = new Array(2);
  var example1_arrays = [practice_struct.V[0][0].slice(0,2), practice_struct.V[0][1].slice(0,3)];
  var example1_objects = [checkOdd(example1_arrays[0]), checkOdd(example1_arrays[1])];
  for (let i = 0; i < 2; i++) {
    if (example1_objects[i].isOdd) {
      example1[i] = "<p style='margin:0'><strong><span style='color:" + sumColours[example1_arrays[i].length-1] + ";'>" + example1_objects[i].oldSum + "</span> is <span style='color:rgb(255,161,0);'>ODD</span></strong></p>" +
                    "<p style='margin:0'><strong>New Sum = <span style='color:" + sumColours[example1_arrays[i].length-1] + ";'>" + example1_objects[i].newSum + "</strong></span></p>";
    } else {
      example1[i] = "<p style='margin:0'><strong><span style='color:" + sumColours[example1_arrays[i].length-1] + ";'>" + example1_objects[i].oldSum + "</span> is EVEN</strong></p>" +
                    "<p style='margin:0'><strong>Sum remains unchanged: <span style='color:" + sumColours[example1_arrays[i].length-1] + ";'>" + example1_objects[i].newSum + "</strong></span></p>";
    }
  }

  // compute sums for example 1 (warning rooms: door 1 room 1, door 2 room 2)
  var example2 = new Array(2);
  var example2_arrays = [practice_struct.V[0][0].slice(0,1), practice_struct.V[0][1].slice(0,2)];
  var example2_objects = [checkOdd(example2_arrays[0]), checkOdd(example2_arrays[1])];
  for (let i = 0; i < 2; i++) {
    if (example2_objects[i].isOdd) {
      example2[i] = "<p style='margin:0'><strong><span style='color:" + sumColours[example2_arrays[i].length-1] + ";'>" + example2_objects[i].oldSum + "</span> is <span style='color:rgb(255,161,0);'>ODD</span></strong></p>" +
                    "<p style='margin:0'><strong>New Sum = <span style='color:" + sumColours[example2_arrays[i].length-1] + ";'>" + example2_objects[i].newSum + "</strong></span></p>";
    } else {
      example2[i] = "<p style='margin:0'><strong><span style='color:" + sumColours[example2_arrays[i].length-1] + ";'>" + example2_objects[i].oldSum + "</span> is EVEN</strong></p>" +
                    "<p style='margin:0'><strong>Sum remains unchanged: <span style='color:" + sumColours[example2_arrays[i].length-1] + ";'>" + example2_objects[i].newSum + "</strong></span></p>";
    }
  }

  var warning_screen = '<div id="warning-screen" style="display:grid; grid-template-rows: auto auto; grid-template-columns: 60% 40%;">' +
    '<h3 style="text-align:center;">Warning indicators:</h3>' +
    '<h3 style="text-align:center;">Door reliability:</h3>' +
    '<div style="display:flex; justify-content:space-around; align-items:center; flex-direction:row; width:100%;">' +
      '<img src="' + seq_array_filenames[0][practice_struct.N[0][0]-1] + '" width="45%" style="padding:1%;">' +
      '<img src="' + seq_array_filenames[1][practice_struct.N[1][0]-1] + '" width="45%" style="padding:1%;">' +
    '</div>' +
    '<div style="height:45%; width:80%; display:flex; align-items:center; flex-direction:column; margin:auto; justify-content:space-evenly;">' +
    '<div style="height:45%; width:80%; display:flex; align-items:centerHey; flex-direction:column; margin:auto; justify-content:space-evenly; "><div id="door-1-prob-text" style="margin: 7% 0 7% 0;"><p style="text-align:center; margin:10px;"><big>Door 1: 50%</big></p><div id="door-1-prob-bar"><div style="width:100%; height:8px; border:3px solid rgb(255, 230, 0); padding:0; margin:0;"><div style="width:50%; height:8px; background:rgb(255, 230, 0); padding:0; margin:0;"></div></div></div></div><div id="door-2-prob-text" style="margin: 7% 0 7% 0;"><p style="text-align:center; margin:10px;"><big>Door 2: 50%</big></p><div id="door-2-prob-bar"><div style="width:100%; height:8px; border:3px solid rgb(255, 230, 0); padding:0; margin:0;"><div style="width:50%; height:8px; background:rgb(255, 230, 0); padding:0; margin:0;"></div></div></div></div></div></div></div>';

  var negator_learning_instructions_pages = [];
  negator_learning_instructions_pages = [
    "<h1>Risky Decision Task</h1>" +
    "<p>You have now learnt the order of rooms (i.e., images) behind Door 1 and Door 2, as well as the points you can gain or lose in each room.</p>" +
    "<div style='display:flex; flex-direction:column; justify-content:center; align-items:center;'>" +
      "<h3>Door 1</h3>" +
      "<div style='display:flex; flex-direction:row'>" +
        "<div style='display:flex; flex-direction:column; align-items:center;'>" +
          "<p><strong>" + seq_array[0][0] + "</strong></p>" +
          "<img src='" + seq_array_filenames[0][0] + "' style='width:50%;'><p>" + practice_struct.V[0][0][0] + "</p></div>" +
        "<div style='display:flex; flex-direction:column; align-items:center;'>" +
          "<p><strong>" + seq_array[0][1] + "</strong></p>" +
          "<img src='" + seq_array_filenames[0][1] + "' style='width:50%;'><p>" + practice_struct.V[0][0][1] + "</p></div>" +
        "<div style='display:flex; flex-direction:column; align-items:center;'>" +
          "<p><strong>" + seq_array[0][2] + "</strong></p>" +
          "<img src='" + seq_array_filenames[0][2] + "' style='width:50%;'><p>" + practice_struct.V[0][0][2] + "</p></div>" +
      "</div>" +
      "<h3>Door 2</h3>" +
      "<div style='display:flex; flex-direction:row'>" +
        "<div style='display:flex; flex-direction:column; align-items:center;'>" +
          "<p><strong>" + seq_array[1][0] + "</strong></p>" +
          "<img src='" + seq_array_filenames[1][0] + "' style='width:50%;'><p>" + practice_struct.V[0][1][0] + "</p></div>" +
        "<div style='display:flex; flex-direction:column; align-items:center;'>" +
          "<p><strong>" + seq_array[1][1] + "</strong></p>" +
          "<img src='" + seq_array_filenames[1][1] + "' style='width:50%;'><p>" + practice_struct.V[0][1][1] + "</p></div>" +
        "<div style='display:flex; flex-direction:column; align-items:center;'>" +
          "<p><strong>" + seq_array[1][2] + "</strong></p>" +
          "<img src='" + seq_array_filenames[1][2] + "' style='width:50%;'><p>" + practice_struct.V[0][1][2] + "</p></div>" +
      "</div>" +
    "</div>",
    "<h1>Risky Decision Task</h1>" +
    "<p>For the remainder of the experiment, you will use the information you have memorised in a risky gambling-style task.</p><p>On each trial, you will be presented with a screen that looks like this:</p>" +
    warning_screen,

    '<h1>Risky Decision Task</h1><h2>Warning Indicators</h2>' +
    '<p>On each trial, there will be <strong>TWO</strong> rooms (one from each path) that will have a <strong>WARNING INDICATOR</strong>. This means that you must apply something called the <strong>ODD RULE</strong> to that room.</p>' +
    '<p style="color:rgb(255,161,0); font-size:1.5rem; text-align:center; line-height:2rem;"><strong>THE ODD RULE</strong>:<br><strong>If the <u>cumulative sum of points</u> <i>in that room</i> is an ODD number (e.g., 1, 3, 5, 7, 9, etc.), then you must multiply the <u>cumulative sum</u> <i>in that room</i> by -1.</strong></p><p>This essentially "flips" the sum of points at that particular stage of your journey along the path (i.e., a positive sum can become negative, or a negative sum can become positive). This sum is then taken to the next room. Hence, if a room with a warning indicator <i>does</i> have an odd sum, then this can significantly alter the course of the points you add up along a path.</p>',

    // EXAMPLE 1
    '<h1>Risky Decision Task</h1><h2>Warning Indicators</h2>' +
    '<p>For example, let`s say that the <strong>second</strong> room behind Door 1 and the <strong>third</strong> room behind Door 2 have warning indicators:</p>' +

    "<div style='display:flex; flex-direction:column; justify-content:center; align-items:center;'>" +
      "<h3>Door 1</h3>" +
      "<div style='display:flex; flex-direction:row'>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[0][0] + "' style='width:50%; border:8px solid transparent;'>" +
          "<p><strong>Value = " + practice_struct.V[0][0][0] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[0] + ";'>" + practice_struct.V[0][0][0] + "</span></strong></p>" +
        "</div>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[0][1] + "' style='width:50%; border:8px solid rgb(255, 161, 0);'>" +
          "<p><strong>Value = " + practice_struct.V[0][0][1] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[0] + ";'>" + practice_struct.V[0][0][0] + "</span>  +  " + practice_struct.V[0][0][1] + "  =  <span style='color:" + sumColours[1] + ";'>" + checkOdd(practice_struct.V[0][0].slice(0,2)).oldSum + "</span></strong></p>" +
          example1[0] +
        "</div>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[0][2] + "' style='width:50%; border:8px solid transparent;'>" +
          "<p><strong>Value = " + practice_struct.V[0][0][2] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[1] + ";'>" + checkOdd(practice_struct.V[0][0].slice(0,2)).newSum + "</span>  +  " + practice_struct.V[0][0][2] + "  =  <span style='color:" + sumColours[2] + ";'>" + [checkOdd(practice_struct.V[0][0].slice(0,2)).newSum,practice_struct.V[0][0][2]].reduce((a,b) => a+b,0)  + "</span></strong></p>" +
          "<p><strong><u>Door 1 Final Sum = " + [checkOdd(practice_struct.V[0][0].slice(0,2)).newSum,practice_struct.V[0][0][2]].reduce((a,b) => a+b,0) + "</u></strong></p>" +
        "</div>" +
      "</div>" +

      "<h3>Door 2</h3>" +
      "<div style='display:flex; flex-direction:row'>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[1][0] + "' style='width:50%; border:8px solid transparent;'>" +
          "<p><strong>Value = " + practice_struct.V[0][1][0] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[0] + ";'>" + practice_struct.V[0][1][0] + "</span></strong></p>" +
        "</div>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[1][1] + "' style='width:50%; border:8px solid transparent;'>" +
          "<p><strong>Value = " + practice_struct.V[0][1][1] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[0] + ";'>" + practice_struct.V[0][1][0] + "</span>  +  " + practice_struct.V[0][1][1] + "  =  <span style='color:" + sumColours[1] + ";'>" + checkOdd(practice_struct.V[0][1].slice(0,2)).oldSum + "</span></strong></p>" +
        "</div>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[1][2] + "' style='width:50%; border:8px solid rgb(255, 161, 0);'>" +
          "<p><strong>Value = " + practice_struct.V[0][1][2] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[1] + ";'>" + checkOdd(practice_struct.V[0][1].slice(0,2)).oldSum + "</span>  +  " + practice_struct.V[0][1][2] + "  =  <span style='color:" + sumColours[2] + ";'>" + [checkOdd(practice_struct.V[0][1].slice(0,2)).oldSum,practice_struct.V[0][1][2]].reduce((a,b) => a+b,0)  + "</span></strong></p>" +
          example1[1] +
          "<p><strong><u>Door 2 Final Sum = " + checkOdd(practice_struct.V[0][1].slice(0,3)).newSum + "</u></strong></p>" +
        "</div>" +
      "</div>" +
    "</div>",

    // EXAMPLE 2
    '<h1>Risky Decision Task</h1><h2>Warning Indicators</h2>' +
    '<p>As another example, let`s say that the <strong>first</strong> room behind Door 1 and the <strong>second</strong> room behind Door 2 have warning indicators:</p>' +

    "<div style='display:flex; flex-direction:column; justify-content:center; align-items:center;'>" +
      "<h3>Door 1</h3>" +
      "<div style='display:flex; flex-direction:row'>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[0][0] + "' style='width:50%; border:8px solid rgb(255, 161, 0);'>" +
          "<p><strong>Value = " + practice_struct.V[0][0][0] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[0] + ";'>" + practice_struct.V[0][0][0] + "</span></strong></p>" +
          example2[0] +
        "</div>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[0][1] + "' style='width:50%; border:8px solid transparent;'>" +
          "<p><strong>Value = " + practice_struct.V[0][0][1] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[0] + ";'>" + checkOdd(practice_struct.V[0][0].slice(0,1)).newSum + "</span>  +  " + practice_struct.V[0][0][1] + "  =  <span style='color:" + sumColours[1] + ";'>" + [checkOdd(practice_struct.V[0][0].slice(0,1)).newSum,practice_struct.V[0][0][1]].reduce((a,b)=>a+b,0) + "</span></strong></p>" +
        "</div>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[0][2] + "' style='width:50%; border:8px solid transparent;'>" +
          "<p><strong>Value = " + practice_struct.V[0][0][2] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[1] + ";'>" + [checkOdd(practice_struct.V[0][0].slice(0,1)).newSum,practice_struct.V[0][0][1]].reduce((a,b)=>a+b,0) + "</span>  +  " + practice_struct.V[0][0][2] + "  =  <span style='color:" + sumColours[2] + ";'>" + [checkOdd(practice_struct.V[0][0].slice(0,1)).newSum,practice_struct.V[0][0][1],practice_struct.V[0][0][2]].reduce((a,b)=>a+b,0)  + "</span></strong></p>" +
          "<p><strong><u>Door 1 Final Sum = " + [checkOdd(practice_struct.V[0][0].slice(0,1)).newSum,practice_struct.V[0][0][1],practice_struct.V[0][0][2]].reduce((a,b)=>a+b,0)  + "</u></strong></p>" +
        "</div>" +
      "</div>" +

      "<h3>Door 2</h3>" +
      "<div style='display:flex; flex-direction:row'>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[1][0] + "' style='width:50%; border:8px solid transparent;'>" +
          "<p><strong>Value = " + practice_struct.V[0][1][0] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[0] + ";'>" + practice_struct.V[0][1][0] + "</span></strong></p>" +
        "</div>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[1][1] + "' style='width:50%; border:8px solid rgb(255, 161, 0);'>" +
          "<p><strong>Value = " + practice_struct.V[0][1][1] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[0] + ";'>" + practice_struct.V[0][1][0] + "</span>  +  " + practice_struct.V[0][1][1] + "  =  <span style='color:" + sumColours[1] + ";'>" + checkOdd(practice_struct.V[0][1].slice(0,2)).oldSum + "</span></strong></p>" +
          example2[1] +
        "</div>" +
        "<div style='display:flex; flex-direction:column; align-items:center; justify-content:flex-start;'>" +
          "<img src='" + seq_array_filenames[1][2] + "' style='width:50%; border:8px solid transparent;'>" +
          "<p><strong>Value = " + practice_struct.V[0][1][2] + "</span></strong></p>" +
          "<p><strong>Sum  =  <span style='color:" + sumColours[1] + ";'>" + checkOdd(practice_struct.V[0][1].slice(0,2)).newSum + "</span>  +  " + practice_struct.V[0][1][2] + "  =  <span style='color:" + sumColours[2] + ";'>" + [checkOdd(practice_struct.V[0][1].slice(0,2)).newSum,practice_struct.V[0][1][2]].reduce((a,b) => a+b,0)  + "</span></strong></p>" +
          "<p><strong><u>Door 2 Final Sum = " + [checkOdd(practice_struct.V[0][1].slice(0,2)).newSum,practice_struct.V[0][1][2]].reduce((a,b) => a+b,0) + "</u></strong></p>" +
        "</div>" +
      "</div>" +
    "</div>",

    '<h1>Risky Decision Task</h1><h2>Warning Indicators</h2>' +
    '<p>As you can see, the <strong>final sum</strong> of points for Door 1 and Door 2 can significantly change depending on <i>which</i> rooms have a warning indicator.</p>' +
    '<p>NOTE: One room from each path will always have a warning indicator.</p>',

    '<h1>Risky Decision Task</h1><h2>Warning Indicators</h2>' +
    '<p>Your task is to use the <strong>warning indicators</strong> shown at the start of each trial to <strong>calculate the final sum of points</strong> for Door 1 and Door 2. You must do this <strong>in your head</strong> by remembering the image order and the value of each image.</p>' +
    warning_screen +
    '<p>You will then choose between TWO options:</p>' +
    '<ul>' +
      '<li><strong>1. AIRLOCK</strong>: If you choose to go to the Airlock, you will go to <i>either</i> Door 1 or Door 2, according to the PROBABILITIES shown on screen.</li>' +
      '<li><strong>2. SUPPLY ROOM</strong>: If you choose to go to the Supply Room, you will receive +1 point.' +
    '</ul>' +
    '<p>You will need to decide whether it is worth risking the Airlock, or if you would be better off going to the Supply Room.</p>',

    '<p style="text-align:center;"><strong>Ready for a practice?</strong></p><p>For the first 6 trials, we will only let you go through the <strong>Airlock</strong>. These are called <strong><span style="color:rgb(203,91,255);">"forced choice"</span></strong> trials. These happen so you get a chance to see how the journey plays out for Door 1 and Door 2. After the forced choice trials, the <strong>images will be hidden</strong>, so you will have to use your <strong>memory</strong></p><p style="text-align:center;">Press "Next" to begin</p><p style="text-align:center;"><small>(or press "Previous" to look back at the instructions)</small></p>'
  ].flat(1);

  var negator_learning_instructions = {
    type: "custom-instructions",
    pages: negator_learning_instructions_pages,
    meg_mode: parameters.exp_variables.meg_mode
  };

  // Trials
  var negator_learning_timeline = [];
  for (let trl = 0; trl < practice_struct.Trial.length; trl++) {

    console.log(practice_struct.Trial.length + ' practice trials');

    if (practice_struct.Forced[trl] == 0) {
      practice_struct.Forced[trl] = 'free';
    }

    // see if there has been a control room change
    let this_trl_prob = [practice_struct.P[trl], 1 - practice_struct.P[trl]];
    let this_trl_negs = [practice_struct.N[trl][0]-1,practice_struct.N[trl][1]-1];

    var control_change = 0;

    var choiceType;
    if (practice_struct.Forced[trl] != 'free')
    {
      choiceType = 0;
    } else {
      choiceType = 'free';
    }

    var practice_choice = {
      type: "test-choice",
      choice: choiceType,
      open_prob: this_trl_prob, // door opening probabilities
      trial_negs: this_trl_negs, // minus one is to get it to 0-indexing
      expected_value: practice_struct.EV[trl],
      transition: practice_struct.Forced[trl]-1,
      is_practice: 2,
      block_num: 0,
      imgs_or_words: "images",
      background: background,
      trial_duration: null, // in seconds - this is the PLANNING WINDOW
      response_window: parameters.timing.resp_window, // in seconds - this is the RESPONSE WINDOW
      feedback_duration: parameters.timing.feedback_window, // in seconds - this is the 'TOO SLOW' response
      miss_penalty: -parameters.values.miss_penalty, // penalty to points for missing the response window
      trial_num: trl,
      stimuli_info: image_info,
      show_score: false,
      control_alert: control_change, // whether there has been a change to PROBABILITY or NEGATORS since the last trial
      choices: parameters.key_responses.left_right_buttons,
      on_finish: function() {
        let filtered_data = JSON.parse(
          jsPsych.data
          .get()
          .filter({
            trial_type: "test-choice",
            practice: 2
          })
          .json()
        ).slice(-1)[0];

        let trial_data = {
          trial_type: "practice-choice",
          trial: filtered_data.trial + 1,
          block: filtered_data.block + 1,
          time_elapsed: filtered_data.time_elapsed,
          resp_name: filtered_data.door_name,
          resp_btn: filtered_data.button,
          rt: filtered_data.rt,
          choice_type: filtered_data.choice_type,
          transition: filtered_data.door_name_transition,
          probability: filtered_data.door_prob,
          negators: filtered_data.negators,
          EV: filtered_data.path_outcomes,
          acc: filtered_data.acc
        };

        if (trial_data.rt == null) {
          trial_data.value = -parameters.values.miss_penalty;
        }

        // if (parameters.exp_variables.meg_mode) {
        //   trial_data.triggers = filtered_data.triggers;
        // }

        // figure out if this was the best choice or not
        if (trial_data.rt == null) {
          trial_data.acc = false;
        } else {

          if (filtered_data.choice_type == 0 | filtered_data.choice_type == 1){
            trial_data.acc = '';
          }
        }

        db.collection("iterations")
          .doc("pilot_v2.5")
          .collection("subjects")
          .doc(uid)
          .collection("5_negator_learning")
          .doc(
            "block" +
            padNumber(trial_data.block) +
            "_trial" +
            padNumber(trial_data.trial) +
            "_choice"
          )
          .set({
            trial_data
          });
      }
    };
    negator_learning_timeline.push(practice_choice);

    for (let stim = 0; stim < 4; stim++) {
      var practice_animation = {
        type: "test-animation",
        choice: practice_struct.Forced[trl],
        meg_mode: parameters.exp_variables.meg_mode,
        stimuli_info: image_info,
        is_practice: 2,
        block_num: 0,
        stim_num: stim,
        background: background,
        trial_duration: parameters.timing.test_animation_dur, // in seconds
        trial_num: trl,
        show_score: false,
        this_val: practice_struct.V[trl],
        negs: this_trl_negs,
        isi: 0, //[0.5, 0.8],
        last_trial: trl == practice_struct.Block.length - 1 && stim == 3,
        trial_end_type: "time-out", // "time-out" = ends at trial_duration, "time-or-response" = ends at trial_duration unless pressed earlier, "response" = ends only when button is pressed
        on_finish: function() {
          let filtered_data = JSON.parse(
            jsPsych.data
            .get()
            .filter({
              trial_type: "test-animation",
              practice: 2
            })
            .json()
          ).slice(-1)[0];
          if (filtered_data != undefined) {
            // if missed the choice phase
            if (filtered_data.hasOwnProperty("trial")) {
              // if blank screens from supply room or there being no explosion
              let trial_data = {
                trial_type: "practice-animation",
                trial: filtered_data.trial + 1,
                block: filtered_data.block + 1,
                time_elapsed: filtered_data.time_elapsed,
                img: filtered_data.img,
                state: filtered_data.state,
                path: filtered_data.path,
                isi: filtered_data.isi,
                negvalue: filtered_data.outcome,
                value: filtered_data.value
              };

              if (parameters.exp_variables.meg_mode) {
                trial_data.triggers = filtered_data.triggers;
              }

              if (trial_data.state == 3) {
                trial_data.state = "explosion";
                delete trial_data.img;
              }

              db.collection("iterations")
                .doc("pilot_v2.5")
                .collection("subjects")
                .doc(uid)
                .collection("5_negator_learning")
                .doc(
                  "block" +
                  padNumber(trial_data.block) +
                  "_trial" +
                  padNumber(trial_data.trial) +
                  "_stim" +
                  stim
                )
                .set({
                  trial_data
                });
            }
          }

          if (
            jsPsych.data
            .get()
            .select("last_trial")
            .values.slice(-1)[0]
          ) {
            var best_answers = jsPsych.data
              .get()
              .filter({
                block: 0
              })
              .select("best_choice").values;
            var actual_answers = jsPsych.data
              .get()
              .filter({
                block: 0
              })
              .select("button")
              .values.map(x => parseInt(x));
            var block_acc = [];
            for (let i = 0; i < best_answers.length; i++) {
              block_acc.push(best_answers[i] == actual_answers[i]);
            }

            var mainDB = db
              .collection("iterations")
              .doc("pilot_v2.5")
              .collection("subjects")
              .doc(uid);
            var update_mainDB = mainDB.set({
              block0_performance: block_acc.filter(x => x == true).length / block_acc.length
            }, {
              merge: true
            });
          }
        }
      };
      negator_learning_timeline.push(practice_animation);
    }

    // detect end of forced practice trials
    var practice_idx = exp_structure.Practice.reduce((out, bool, index) => bool ? out.concat(index) : out, []);
    var end_forced = [];
    for (let i of practice_idx) {
      if (exp_structure.Forced[i] == false && exp_structure.Forced[i - 1] == true) {
        end_forced = i;
      }
    }

    if (trl == end_forced-1) {
      var main_instructions_checkpoint = {
        type: "custom-instructions",
        pages: ['<p style="text-align:center;">The supply room is now open, so we will now hide the images when you visit each room. Instead, a "?" will be displayed. Try to remember and imagine the images.</p>'],
        background: background,
        key_forward: parameters.key_responses.left_right_buttons[1],
        meg_mode: parameters.exp_variables.meg_mode
      };
      negator_learning_timeline.push(main_instructions_checkpoint);
    }
  }

  /* check if they failed the learning phase */
  var check_negator_learning = {
    timeline: end_message,
    conditional_function: function() {
      var data = jsPsych.data
        .get()
        .filter({
          practice: 2
        })
        .select("rt").values;
      if (
        data.filter(x => x == null).length / data.length >=
        parameters.thresholds.practice_threshold
      ) {
        return false;
      } else {
        return true;
      }
    }
  };

  ////////////////
  // MAIN TEST //
  ///////////////

  var bonus_money = '';
  if (parameters.exp_variables.meg_mode){
    bonus_money = '10';
  } else {
    bonus_money = '5';
  }

  var main_instructions_pages;
  if (parameters.exp_variables.meg_mode)
  {
    main_instructions_pages = [
    "<h1>Survival</h1>" +
    "<p>You are now ready to complete the experiment. Please wait for the experimenter. </p>"
  ];
} else {
  main_instructions_pages = [
  "<h1>Survival</h1>" +
  "<p>You are now ready to start the main experiment.</p>" +
  "<p><strong>Remember: you can earn a <span style='color:rgb(0,255,0);'>BONUS of up to £" + bonus_money + "</span> if you complete the experiment! The more oxygen points you earn each time, the higher your bonus will be.</strong><p>" +
  '<p>Please take note of some <strong><span style="color:rgb(255,10,10);">important changes</span></strong> for the main experiment:</p>' +
  "<ul style='line-height:2rem;'><li>You will only have <strong>" +
  parameters.timing.choice_dur +
  " SECONDS</strong> to PLAN which choice to make each time. A timer will be shown. If the timer runs out, you will LOSE 1 oxygen point.</li>" +
  "<li>The <strong>NUMBER OF POINTS</strong> you pick up in each room will <strong>CHANGE</strong> at the start of each block in the experiment. Keep track of this by paying attention to the value of each image when you go through the Airlock at the beginning of the block (this is what the <strong><span style='color:rgb(203,91,255);'>'forced choice'</span></strong> trials are for).</li>" +
  "</ul>",

  "<p style='text-align:center; color:#FFAE57'><strong>YOUR TASK:</strong> Calculate the total sum of points for Door 1 and Door 2, according to which rooms are HAZARDOUS. Then, depending on whether Door 1 or Door 2 is more likely to be open, decide if it's worth the risk to go to the airlock or to get 1 point from the supply room.</p><p style='text-align:center'><strong>Ready?</strong></p><p style='text-align:center'>Press 'Next' to begin</p><p style='text-align:center'><small>(or press 'Previous' to look back at the instructions)</small></p>"
];
}

  var main_instructions = {
    type: "custom-instructions",
    pages: main_instructions_pages,
    background: background,
    meg_mode: parameters.exp_variables.meg_mode
  };

  var nblocks = exp_structure.Block[exp_structure.Block.length - 1];

  var test_timeline = [];
  var ntrials = 0;
  for (let block = 0; block < nblocks; block++) {

    var test_structure = {};
    for (let key in exp_structure) {
      test_structure[key] = [];
      for (let i = 0; i < exp_structure[key].length; i++) {
        if (exp_structure.Block[i] == block+1 && exp_structure.Practice[i] == false) {
          test_structure[key].push(exp_structure[key][i]);
        }
      }
    }

      var ntrials = test_structure.Trial.length; // trials per block (using block 1 as example)
    console.log(ntrials + ' test trials');

    for (let trl = 0; trl < ntrials; trl++) {

      if (test_structure.Forced[trl] == 0) {
        test_structure.Forced[trl] = 'free';
      }

      // see if there has been a control room change
      let this_trl_prob = [test_structure.P[trl], 1 - test_structure.P[trl]];
      let this_trl_negs = [test_structure.N[trl][0]-1,test_structure.N[trl][1]-1];

      var control_change = 0;

      var choiceType;
      if (test_structure.Forced[trl] != 'free')
      {
        choiceType = 0;
      } else {
        choiceType = 'free';
      }

      var test_choice = {
        type: "test-choice",
        background: background,
        choice: choiceType,
        imgs_or_words: "words",
        meg_mode: parameters.exp_variables.meg_mode,
        expected_value: test_structure.EV[trl],
        transition: test_structure.Forced[trl]-1,
        show_score: false,
        open_prob: this_trl_prob, // door opening probabilities
        trial_negs: this_trl_negs,
        is_practice: 0,
        stimuli_info: image_info,
        trial_duration: parameters.timing.choice_dur, // in seconds - this is the PLANNING WINDOW
        response_window: parameters.timing.resp_window, // in seconds - this is the RESPONSE WINDOW
        feedback_duration: parameters.timing.feedback_window, // in seconds - this is the 'TOO SLOW' response
        miss_penalty: -parameters.values.miss_penalty, // penalty to points for missing the response window
        trial_num: trl,
        block_num: block,
        control_alert: control_change,
        choices: parameters.key_responses.left_right_buttons,
        on_finish: function() {
          let filtered_data = JSON.parse(
            jsPsych.data
            .get()
            .filter({
              trial_type: "test-choice",
              practice: 0
            })
            .json()
          ).slice(-1)[0];

          let trial_data = {
            trial_type: "test-choice",
            trial: filtered_data.trial + 1,
            block: filtered_data.block + 1,
            time_elapsed: filtered_data.time_elapsed,
            resp_name: filtered_data.door_name,
            resp_btn: filtered_data.button,
            rt: filtered_data.rt,
            choice_type: filtered_data.choice_type,
            transition: filtered_data.door_name_transition,
            probability: filtered_data.door_prob,
            negators: filtered_data.negators,
            EV: filtered_data.path_outcomes,
            acc: filtered_data.acc
          };

          if (trial_data.rt == null) {
            trial_data.value = -parameters.values.miss_penalty;
          }

          if (parameters.exp_variables.meg_mode) {
            trial_data.triggers = filtered_data.triggers;
          }

          // figure out if this was the best choice or not
          if (trial_data.rt == null) {
            trial_data.acc = false;
          } else {
            if (filtered_data.choice_type == 0 | filtered_data.choice_type == 1){
              trial_data.acc = '';
            }
          }

          db.collection("iterations")
            .doc("pilot_v2.5")
            .collection("subjects")
            .doc(uid)
            .collection("6_test")
            .doc(
              "block" +
              padNumber(trial_data.block) +
              "_trial" +
              padNumber(trial_data.trial) +
              "_choice"
            )
            .set({
              trial_data
            });
        }
      };
      test_timeline.push(test_choice);

      for (let stim = 0; stim < 4; stim++) {
        var test_animation = {
          type: "test-animation",
          background: background,
          choice: test_structure.Forced[trl],
          meg_mode: parameters.exp_variables.meg_mode,
          is_practice: 0,
          stimuli_info: image_info,
          stim_num: stim,
          trial_duration: parameters.timing.test_animation_dur, // in seconds
          trial_num: trl,
          this_val: test_structure.V[trl],
          negs: this_trl_negs,
          last_trial: trl == ntrials - 1 && stim == 3,
          isi: 0, //[0.5, 0.8],
          block_num: block,
          trial_end_type: "time-out", // "time-out" = ends at trial_duration, "time-or-response" = ends at trial_duration unless pressed earlier, "response" = ends only when button is pressed
          on_finish: function() {
            let filtered_data = JSON.parse(
              jsPsych.data
              .get()
              .filter({
                trial_type: "test-animation",
                practice: 0
              })
              .json()
            ).slice(-1)[0];

            if (filtered_data != undefined) {
              // if missed the choice phase
              if (filtered_data.hasOwnProperty("trial")) {
                // if blank screens from supply room or there being no explosion
                let trial_data = {
                  trial_type: "test-animation",
                  trial: filtered_data.trial + 1,
                  block: filtered_data.block + 1,
                  time_elapsed: filtered_data.time_elapsed,
                  img: filtered_data.img,
                  state: filtered_data.state,
                  path: filtered_data.path,
                  isi: filtered_data.isi,
                  negvalue: filtered_data.outcome,
                  value: filtered_data.value
                };

                if (parameters.exp_variables.meg_mode) {
                  trial_data.triggers = filtered_data.triggers;
                }

                db.collection("iterations")
                  .doc("pilot_v2.5")
                  .collection("subjects")
                  .doc(uid)
                  .collection("6_test")
                  .doc(
                    "block" +
                    padNumber(trial_data.block) +
                    "_trial" +
                    padNumber(trial_data.trial) +
                    "_stim" +
                    stim
                  )
                  .set({
                    trial_data
                  });
              }
            }

            if (
              jsPsych.data
              .get()
              .select("last_trial")
              .values.slice(-1)[0]
            ) {
              var this_block = jsPsych.data
                .get()
                .select("block")
                .values.slice(-1)[0];

              var best_answers = jsPsych.data
                .get()
                .filter({
                  block: this_block
                })
                .select("best_choice").values;
              var actual_answers = jsPsych.data
                .get()
                .filter({
                  block: this_block
                })
                .select("button")
                .values.map(x => parseInt(x));
              var block_acc = [];
              for (let i = 0; i < best_answers.length; i++) {
                block_acc.push(best_answers[i] == actual_answers[i]);
              }

              var obj = {};
              obj["block" + this_block + "_performance"] =
                block_acc.filter(x => x == true).length / block_acc.length;

              var mainDB = db
                .collection("iterations")
                .doc("pilot_v2.5")
                .collection("subjects")
                .doc(uid);
              var update_mainDB = mainDB.set(obj, {
                merge: true
              });
            }
          }
        };
        test_timeline.push(test_animation);
      }
    }

    var continueText = parameters.exp_variables.meg_mode ? 'The experimenter will be with you shortly.' : 'Press the button below when you\'re ready to continue.';
    if (block < nblocks) {
      var end_block = {
        type: "custom-instructions",
        meg_mode: parameters.exp_variables.meg_mode,
        on_start: function(trial) {
          var last_choice = jsPsych.data.get().filter({
            trial_type: 'test-choice'
          });
          var last_block = jsPsych.data.get().select('block').values.pop();
          var block_acc = jsPsych.data.get().filter({block: last_block, choice_type: "free"}).select('acc').mean();
          trial.stimuli = ['<h3 style="text-align:center;">End of block ' + (last_block+1) + '</h3>' +
          '<p style="text-align:center;">You made the best choice <strong>' + (Math.round(block_acc*100)) + '%</strong> of the time.</p>' +
          '<p style="text-align:center; color:rgb(0,255,0);"><strong><big>You earnt £' + (Math.ceil(block_acc*10)/10).toFixed(2) + '!</big></strong></p>'];

          var mainDB = db
            .collection("iterations")
            .doc("pilot_v2.5")
            .collection("subjects")
            .doc(uid);

          var bonuses = {};
          bonuses['block_' + (last_block+1) + '_bonus'] = (Math.ceil(block_acc*10)/10).toFixed(2);

          var update_mainDB;
          update_mainDB = mainDB.set({bonuses}, {merge: true});
        },
        pages: [
          '<p style="text-align:center;">Please take a short break.</p><p style="text-align:center;">' + continueText + '</p>'
        ],
        button_label_next: "Continue",
        background: background
      };
      test_timeline.push(end_block);
    }
  }

  ///////////////////////
  // FINALISE TIMELINE //
  ///////////////////////

  var timeline = [];
  if (parameters.exp_variables.meg_mode) {
    timeline = [
      functional_localiser,
      instructions,
      exploration_timeline,
      repeat_exploration,
      abort_exploration,
      value_learning_instructions,
      value_learning_timeline,
      repeat_value_learning,
      abort_value_learning,
      negator_learning_instructions,
      negator_learning_timeline,
      check_negator_learning,
      main_instructions,
      test_timeline,
      end_message
    ].flat(Infinity);
  } else {
    timeline = [
      questionnaire_timeline,
      instructions,
      exploration_timeline,
      repeat_exploration,
      abort_exploration,
      value_learning_instructions,
      value_learning_timeline,
      repeat_value_learning,
      abort_value_learning,
      negator_learning_instructions,
      negator_learning_timeline,
      check_negator_learning,
      main_instructions,
      test_timeline,
      end_message
    ].flat(Infinity);
  }

  return [timeline, seq_array];
};

// Auxillary Functions

// add zeros in front of number less than 10
window.padNumber = function(number) {
  var numberStr = "";
  if (number < 10) {
    numberStr = "0" + number;
  } else {
    numberStr = number.toString();
  }
  return numberStr;
};

// calculate timings
window.linspace = function(start, stop, increment) {
  var arr = [start];
  var this_dec =
    Math.floor(increment) === increment ?
    0 :
    increment.toString().split(".")[1].length || 0;

  var keep_going = true;
  while (keep_going) {
    var new_val = arr[arr.length - 1] + increment;
    new_val = parseFloat(new_val.toFixed(this_dec));
    if (new_val <= stop) {
      arr.push(new_val);
    } else {
      keep_going = false;
    }
  }
  return arr;
};
