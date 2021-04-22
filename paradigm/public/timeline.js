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
      "<h1>The Awakening</h1>" +
      "<p>You are a passenger aboard a spaceship making the long journey to another planet. You and the other passengers have all been cryogenically frozen and set to wake up on the day of landing.</p><p>You have woken up earlier than you were supposed to. The doors of your cryogenic chamber open and you crawl out. You look around. The other passengers are all still asleep.</p><p>The spaceship's control panel is in the middle of the room. You walk up to it and inspect the screen. It shows that the date is 5th November, 2173. You have woken up several days before your scheduled arrival date.</p><p>You decide to explore the rooms of the spaceship to see what you can find.</p>",

      "<h1>The Awakening</h1>" +
      "<p>In front of you is a large door. You walk up to it, haul it open, and walk through. The door slides shut behind you and locks.</p><p>It seems that you have entered into a small airlocked room with <strong>2 DOORS</strong> that lead to other parts of the ship.</p><p>You now have no choice but to explore the paths behind each door. You will memorise the order of the rooms along each path. Each room will be represented by an <strong>IMAGE</strong>.</p><p>Later, you will test your memory for the order of images along each path.</p>",

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
        end_pause: 500,
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
      //exploration_timeline.push(exploration_choice);

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
        //exploration_timeline.push(exploration_animation);
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
        var tmp_params = loadParameters();
        if (trial_data.trial == tmp_params.trial_numbers.memory_test_reps * 2) {
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
  if (parameters.exp_variables.meg_mode) {
    value_learning_instructions_pages = [
      "<h1>Points</h1>" +
      "<p style='text-align:center;'>You must <strong>MEMORISE THE NUMBER OF OXYGEN POINTS</strong> that you <span style='color:rgb(0,255,0);'><strong>GAIN</strong></span> or <span style='color:rgb(255,20,20);'><strong>LOSE</strong></span> in each room.</p><p>Keep an eye on the number of oxygen points you are <strong>CURRENTLY CARRYING</strong> as you progress through each room. This reflects the 'cumulative sum' of points, which you will need to be able to recall easily later on. This will be shown at the bottom of the screen as your journey plays out.</p>",
      "<p style='text-align:center; color:#FFAE57'><strong>YOUR TASK:</strong> Memorise the number of points gained or lost in each room.</p><p style='text-align:center;'><strong>Ready?</strong></p><p style='text-align:center;'>Press 'Next' to begin</p><p style='text-align:center;'><small>(or press 'Previous' to look back at the instructions)</small></p>"
    ];
  } else {
    value_learning_instructions_pages = [
      "<h1>Red Alert</h1>" +
      "<p>You are sitting in the <strong>CONTROL ROOM</strong> when suddenly the room flashes red and alarms start wailing.</p><p>You go to the control panel and read the message...</p>",
      "<div id='warning-screen'><h2>WARNING</h2><p>ERROR #0624: INSUFFICIENT OXYGEN SUPPLY</p><p>Oxygen levels are dangerously low. Fit spacesuit and oxygen canister IMMEDIATELY. Ensure all passengers are CRYOGENICALLY FROZEN.</p></div><p>On the control panel, the light next to <strong>SUPPLY ROOM</strong> is flashing. You remember that the supply room is located at the end of a long, winding corridor. You make your way there...</p>",
      "<p>You move as fast as you can down the long corridor towards the <strong>SUPPLY ROOM</strong>, gasping for breath as the oxygen levels plummet.</p><p>Finally, you reach the room. Inside, the space suit is hanging on a hook on the wall. You desperately pull it off the hook and climb into the suit, zipping it up and fitting the helmet over your head.</p><p>Around you there are stacks of crates, full of supplies and - oxygen canisters!</p><p>You grab the canister nearest to you and fit it into the slot on the space suit. It clicks into place and a little light next to it turns green.</p>",
      "<p>Warmth rushes over you as you breathe in the steady flow of oxygen. You check the indicator light on the oxygen pack.</p><p>It has already turned yellow.</p><p>You must get back to the control room as fast as you can, where the oxygen levels are more stable.</p>",
      "<p>You arrive back at the control room, just as that oxygen canister runs out and the light turns red.<p>You could go back to the <strong>SUPPLY ROOM</strong>... but the room is far away.</p><p>Instead, you decide to explore the paths behind the two doors in the airlocked room that you learnt earlier to collect more oxygen.</p>",
      "<p style='text-align:center;'>You must <strong>MEMORISE THE NUMBER OF OXYGEN POINTS</strong> that you <span style='color:rgb(0,255,0);'><strong>GAIN</strong></span> or <span style='color:rgb(255,20,20);'><strong>LOSE</strong></span> in each room.</p><p>Keep an eye on the number of oxygen points you are <strong>CURRENTLY CARRYING</strong> as you progress through each room. This reflects the 'cumulative sum' of points, which you will need to be able to recall easily later on. This will be shown at the bottom of the screen as your journey plays out.</p>",
      "<p style='text-align:center; color:#FFAE57'><strong>YOUR TASK:</strong> Memorise the number of points gained or lost in each room.</p><p style='text-align:center;'><strong>Ready?</strong></p><p style='text-align:center;'>Press 'Next' to begin</p><p style='text-align:center;'><small>(or press 'Previous' to look back at the instructions)</small></p>"
    ];
  }

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
      is_practice: 1,
      end_pause: 500,
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
        var tmp_params = loadParameters();
        if (trial_data.trial == tmp_params.trial_numbers.value_test_reps * 2) {
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

  var neg_filenames = [
    seq_array_filenames[0][practice_struct.N[0][0]-1],
    seq_array_filenames[1][practice_struct.N[0][1]-1]
  ];
  var neg_names = [
    seq_array[0][practice_struct.N[0][0]-1],
    seq_array[1][practice_struct.N[0][1]-1]
  ];

  var prev_prob = [0.5, 0.5];
  var prev_negs = [-1, -1];
  var control_change = null;

  var indicator_screen =
    "<p>As an example, you may see the following screen before making your choice. Use this to mentally plan each possible outcome and then make the best decision (i.e. the choice that is the most likely to get you the most points):</p>" +
    '<div id="warning-screen" style="display:grid; grid-template-rows: auto auto; grid-template-columns: 60% 40%;">' +
    '<h3 style="text-align:center;">Warning indicators:</h3>' +
    '<h3 style="text-align:center;">Door reliability:</h3>' +
    '<div style="display:flex; justify-content:space-around; align-items:center; flex-direction:row; width:100%;">' +
    //'<div><p class="neg-cue">' +
    //seq_array[0][practice_struct.N[0][0]-1] + '</p></div>' +
    //'<div><p class="neg-cue">' +
    //seq_array[1][practice_struct.N[1][0]-1] + '</p></div>' +
    //'</div>' +
    '<img src="' + seq_array_filenames[0][practice_struct.N[0][0]-1] + '" width="45%" style="padding:1%;">' +
    '<img src="' + seq_array_filenames[1][practice_struct.N[1][0]-1] + '" width="45%" style="padding:1%;"></div>' +
    '<div style="height:45%; width:80%; display:flex; align-items:center; flex-direction:column; margin:auto; justify-content:space-evenly;">' +
    '<div style="height:45%; width:80%; display:flex; align-items:centerHey; flex-direction:column; margin:auto; justify-content:space-evenly; "><div id="door-1-prob-text" style="margin: 7% 0 7% 0;"><p style="text-align:center; margin:10px;"><big>Door 1: 50%</big></p><div id="door-1-prob-bar"><div style="width:100%; height:8px; border:3px solid rgb(255, 230, 0); padding:0; margin:0;"><div style="width:50%; height:8px; background:rgb(255, 230, 0); padding:0; margin:0;"></div></div></div></div><div id="door-2-prob-text" style="margin: 7% 0 7% 0;"><p style="text-align:center; margin:10px;"><big>Door 2: 50%</big></p><div id="door-2-prob-bar"><div style="width:100%; height:8px; border:3px solid rgb(255, 230, 0); padding:0; margin:0;"><div style="width:50%; height:8px; background:rgb(255, 230, 0); padding:0; margin:0;"></div></div></div></div></div></div></div>' +
    '<p>In this example, the <strong>"' + seq_array[0][practice_struct.N[0][0]-1] + '"</strong> room and the <strong>"' + seq_array[1][practice_struct.N[1][0]-1] + '"</strong> room can potentially reverse your points (but only if the total number of points you have in the room is an ODD number).</p><p>There is also a 50-50 chance that Door 1 or Door 2 will be open (see next page for more details).</p>';

  function oddeven_text(path) {
    var these_vals = practice_struct.V[0][path].slice(0);
    var these_neg_vals = practice_struct.nV[0][path].slice(0);

    var brackets0 = these_neg_vals[0];
    var brackets1 = these_neg_vals[1];
    var brackets2 = these_neg_vals[2];
    if (practice_struct.N[0][path]-1 == 0) {
      brackets0 = these_vals[0] + '<span style="color:rgb(248,161,2);"> ... ' + these_neg_vals[0] + '</span>';
    } else if (practice_struct.N[0][path]-1 == 1) {
      brackets1 = (these_vals[0] + these_vals[1]) + '<span style="color:rgb(248,161,2);"> ... ' + these_neg_vals[1] + '</span>';
    } else if (practice_struct.N[0][path]-1 == 2) {
      brackets2 = (these_vals[0] + these_vals[1] + these_vals[2]) + '<span style="color:rgb(248,161,2);"> ... ' + these_neg_vals[2] + '</span>';
    }

    var images = "";
    var text =
      '<div id="value-container" style="display:flex; flex-direction:row; flex-wrap:nowrap; justify-content:space-around;"><div style="display:flex; flex-wrap:nowrap; flex-direction:row; justify-content:space-between; width:33%; margin:auto;"><p style="text-align:left;">Room points: </p><p style="width:100%; text-align:center;">' +
      these_vals[0] +
      '</p></div><p style="width:33%; text-align:center;">' +
      these_vals[1] +
      '</p><p style="width:33%; text-align:center;">' +
      these_vals[2] +
      "</p></div>" +
      '<div id="cumulative-value-container" style="display:flex; flex-direction:row; flex-wrap:nowrap; justify-content:space-around;"><div style="display:flex; flex-wrap:nowrap; flex-direction:row; justify-content:space-between; width:33%;"><p style="text-align:left;">Gained/Lost: </p><p style="width:100%; text-align:center;"><strong>' +
      brackets0 +
      '</strong></p></div><p style="width:33%; text-align:center;">' +
      these_neg_vals[0] +
      " + " +
      these_vals[1] +
      " = <strong>" +
      brackets1 +
      '</strong></p><p style="width:33%; text-align:center;">' +
      these_neg_vals[1] +
      " + " +
      these_vals[2] +
      " = <strong>" +
      brackets2 +
      "</strong></p></div>";

    var oddeven = "";
    var badpath = these_neg_vals[practice_struct.N[path]-1] % 2;
    for (let i = 0; i < 3; i++) {
      var odd_label = Math.abs(these_neg_vals[i]) % 2 ? "ODD" : "EVEN";

      if (i == practice_struct.N[0][path]-1) {
        images +=
          '<img src="' +
          seq_array_filenames[path][i] +
          '" alt="' +
          seq_array[path][i] +
          '" style="border:8px solid rgb(255, 161, 0);">';
        if (odd_label == "ODD") {
          if (i == 0) {
            oddeven +=
              '<p style="width:100%; text-align:center; color:rgb(255, 161, 0);"><strong>ODD</strong>: total sum of points will be reversed!</p></div>';
          } else {
            oddeven +=
              '<p style="width:33%; text-align:center; color:rgb(255, 161, 0);"><strong>ODD</strong>: total sum of points will be reversed!</p>';
          }
        } else {
          if (i == 0) {
            oddeven +=
              '<p style="width:100%; text-align:center; color:rgb(0, 255, 0);"><strong>' +
              odd_label +
              "</strong></p></div>";
          } else {
            oddeven +=
              '<p style="width:33%; text-align:center; color:rgb(0, 255, 0);"><strong>' +
              odd_label +
              "</strong></p>";
          }
        }
      } else {
        images +=
          '<img src="' +
          seq_array_filenames[path][i] +
          '" alt="' +
          seq_array[path][i] +
          '">';
        if (i == 0) {
          oddeven +=
            '<p style="width:100%; text-align:center;"><strong>' +
            odd_label +
            "</strong></p></div>";
        } else {
          oddeven +=
            '<p style="width:33%; text-align:center;"><strong>' +
            odd_label +
            "</strong></p>";
        }
      }
    }

    var outcome_colour =
      these_neg_vals[2] < 0 ? "rgb(255,0,0)" : "rgb(0,255,0)";

    text +=
      '<div id="result-container" style="display:flex; flex-direction:row; flex-wrap:nowrap; justify-content:space-around;"><div style="display:flex; flex-wrap:nowrap; flex-direction:row; justify-content:space-between; width:33%;"><p style="text-align:left;">Reversal: </p>' +
      oddeven +
      "</div>" +
      '<div id="outcome-container" style="display:flex; flex-direction:row; flex-wrap:nowrap; justify-content:space-around;"><div style="display:flex; flex-wrap:nowrap; flex-direction:row; justify-content:space-between; width:33%;"><p style="text-align:left;"><strong><big>Final outcome:</strong></big></p><p style="width:100%; text-align:center; color:rgb(255, 161, 0);">&nbsp</p></div><p style="width:33%; text-align:center;">&nbsp</p><p style="width:33%; text-align:center;"><span style="color:' +
      outcome_colour +
      ';"><strong><big>' +
      these_neg_vals[2] +
      "</big></strong></span></p></div>";

    return [images, text];
  }

  var imagestext1 = oddeven_text(0);
  var imagestext2 = oddeven_text(1);

  var negator_learning_instructions_pages = [];
  if (parameters.exp_variables.meg_mode) {
    negator_learning_instructions_pages = [
      "<h1>Main Task Practice</h1>" +
      "<p>Please wait for the experimenter to explain the next part of the task.</p>",

      indicator_screen,
      "<p>To plan for Door 1, you would need to calculate the following:</p>" +
      '<div id="value-image-container2" style="opacity:1;">' +
      imagestext1[0] +
      "</div>" +
      imagestext1[1],

      "<p>To plan for Door 2, you would need to calculate the following:</p>" +
      '<div id="value-image-container2" style="opacity:1;">' +
      imagestext2[0] +
      "</div>" +
      imagestext2[1],

      '<div><p style="align-text:center;"><strong>You can earn a <span style="color:rgb(0,255,0);">BONUS of up to £10</span> if you make the best possible choices.</strong></p><p>To do this, you will have to calculate the total points for Door 1 and Door 2, and then decide whether you should risk going through the airlock, or instead go to the supply room.<p>' +
      '<p style="text-align:center;">&nbsp</p><p style="text-align:center;">Please wait for the experimenter to begin the program.</p><p style="text-align:center;"><strong>Remember to stay very still during scanning.</strong></p></div>'
    ].flat(1);
  } else {
    negator_learning_instructions_pages = [
      "<h1>Survival</h1>" +
      "<p>Here is a recap of what you have learnt so far.</p>" +
      "<h2>Aim</h2><p>Your aim is to collect as many oxygen points as possible so that you can maximise your chances of survival when the ship finally reaches its destination.</p><p>You will receive <strong><span style='color:rgb(0,255,0);'>real bonus money</span></strong> based on how well you maximise your points!</p>" +
      "<h2>Choices</h2><p>You have learnt that there are two paths - Door 1 and Door 2 (which are accessed via the AIRLOCK) - and the SUPPLY ROOM. Each day, you can choose to go to either the:</p>" +
      '<ul style="text-align:left;">' +
      "<li><strong>Airlock</strong> (where you can access Door 1 or Door 2)</li>" +
      "<li><strong>Supply Room</strong> (where you will always return with 1 oxygen point)</li>" +
      "</ul>",

      "<h1>New Rules</h1><h2>Door Probabilities</h2>" +
      '<p>If you choose to go through the airlock...' +
      '<ul style="padding-left:25px;"><li>Only <strong>ONE</strong> of the doors will ever be open at a time.</li>' +
      "<li>You will not know WHICH door is open until you get to the airlock.</li>" +
      "<li>To help you decide, you will be told the <strong>PROBABILITY</strong> that Door 1 vs Door 2 will be open (e.g., 90% chance Door 1, 10% chance Door 2)</li>" +
      "<li>The door to the supply room will <i>always</i> be open, and you will <i>always</i> retrieve 1 oxygen point.</li></ul>" +
      "</ul></p>",

      "<h1>New Rules</h1><h2>Hazardous Rooms</h2>" +
      '<p>Each day, you will receive a warning that some of the rooms can <strong>REVERSE</strong> your points:</p>' +
      '<ul style="text-align:left;">' +
      '<li>There will always be TWO "hazardous" rooms; one behind Door 1 and one behind Door 2.</li>' +
      "<li>If you <strong>LEAVE</strong> one of these rooms while carrying an <strong>ODD</strong> number of oxygen points, the points will <strong>REVERSE</strong>. This means that if they were positive then they will become <strong><span style='color:rgb(255,0,0);'>NEGATIVE</span></strong>, and if they were negative then they will become <strong><span style='color:rgb(0,255,0);'>POSITIVE</span></strong>.</li>" +
      "</ul>" +
      '<p>Proceed to the next page to see a worked-through example. Please ask the experimenter for any help.</p>',


      indicator_screen,
      "<p>To plan for Door 1, you would need to calculate the following:</p>" +
      '<div id="value-image-container2" style="opacity:1;">' +
      imagestext1[0] +
      "</div>" +
      imagestext1[1],

      "<p>To plan for Door 2, you would need to calculate the following:</p>" +
      '<div id="value-image-container2" style="opacity:1;">' +
      imagestext2[0] +
      "</div>" +
      imagestext2[1],

      '<div><p><strong>You can earn a <span style="color:rgb(0,255,0);">BONUS of up to £5</span> if you make the best possible choices!</strong></p><p>To do this, you will have to calculate the total points for Door 1 and Door 2, and then decide whether you should risk going through the airlock, or instead go to the supply room.<p>' +
      '<p>Notes:</p>' +
      '<ul style="text-align:left;">' +
      '<li><strong>The supply room will be CLOSED for the first few trials</strong>, so you will have to go through the AIRLOCK. These are called <span style="color:rgb(154, 72, 241)"><strong>forced choice</strong></span> trials. After the supply room OPENS, you will be able to choose between the airlock and the supply room.</li>' +
      '<li>We want you to <i>mentally imagine</i> the pictures for each room as much as possible. To encourage this, we will only SHOW you the images on the first few trials (i.e., the <span style="color:rgb(154, 72, 241)"><strong>forced choice</strong></span> trials). The rest of the time, the images will be replaced by "?" symbols.</li></ul>' +
      '<p>Before we do the main experiment, we will let you do a practice so that you are as familiar as possible with the game.</p>',

      '<p style="text-align:center;"><strong>Ready for a practice?</strong></p><p style="text-align:center;">Press "Next" to begin</p><p style="text-align:center;"><small>(or press "Previous" to look back at the instructions)</small></p>'
    ].flat(1);
  }


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

    let prob_change = [
      this_trl_prob[0] != prev_prob[0],
      this_trl_prob[1] != prev_prob[1]
    ];
    let neg_change = [
      this_trl_negs[0] != prev_negs[0],
      this_trl_negs[1] != prev_negs[1]
    ];

    // if (prob_change.includes(true) && neg_change.includes(true)) {
    //   control_change = 3;
    // } else if (prob_change.includes(true)) {
    //   control_change = 1;
    // } else if (neg_change.includes(true)) {
    //   control_change = 2;
    // } else {
      control_change = 0;
    // }

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
      end_pause: 500,
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
        isi: [0.5, 0.8],
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

    prev_prob = this_trl_prob.slice(0);
    prev_negs = this_trl_negs.slice(0);
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
  "<li>The <strong>NUMBER OF OXYGEN POINTS</strong> you pick up in each room will gradually <strong>CHANGE</strong> across the experiment. Keep track of this.</li>" +
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

      let prob_change = [
        this_trl_prob[0] != prev_prob[0],
        this_trl_prob[1] != prev_prob[1]
      ];
      let neg_change = [
        this_trl_negs[0] != prev_negs[0],
        this_trl_negs[1] != prev_negs[1]
      ];

      if (prob_change.includes(true) && neg_change.includes(true)) {
        control_change = 3;
      } else if (prob_change.includes(true)) {
        control_change = 1;
      } else if (neg_change.includes(true)) {
        control_change = 2;
      } else {
        control_change = 0;
      }

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
        expected_value: test_structure.EV[trl],
        transition: test_structure.Forced[trl]-1,
        show_score: false,
        open_prob: this_trl_prob, // door opening probabilities
        trial_negs: this_trl_negs,
        is_practice: 0,
        stimuli_info: image_info,
        end_pause: 500,
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
          is_practice: 0,
          stimuli_info: image_info,
          stim_num: stim,
          trial_duration: parameters.timing.test_animation_dur, // in seconds
          trial_num: trl,
          this_val: test_structure.V[trl],
          negs: this_trl_negs,
          last_trial: trl == ntrials - 1 && stim == 3,
          isi: [0.5, 0.8],
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

      prev_prob = this_trl_prob.slice(0);
      prev_negs = this_trl_negs.slice(0);

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
