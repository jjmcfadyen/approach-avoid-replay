/* --- SET UP USER ID AND DATA STORAGE ON FIRESTORE --- */
/*jshint esversion: 6 */

// Sentry.init({
//   dsn: "https://fccb0dfb5ee34985a3297e8fd1fab8b1@sentry.io/1793765",
//   beforeSend(event, hint) {
//     // Check if it is an exception, and if so, show the report dialog
//     if (event.exception) {
//       Sentry.showReportDialog({
//         eventId: event.event_id,
//         title: "It looks like we're having issues.",
//         subtitle: "The researcher has been notified.",
//         user: {
//           name: "anonymous",
//           email: "anonymous@anonymous.com"
//         }
//       });
//     }
//     return event;
//   },
//   release: "risky-replay@2.2.8"
// });

import "/jspsych-6.1.0/jspsych.js";

import "/scripts/custom-survey.js";
import "/scripts/custom-instructions.js";
import "/scripts/exploration-choice.js";
import "/scripts/exploration-animation.js";
import "/scripts/exploration-test.js";
import "/scripts/value-test.js";
import "/scripts/test-choice.js";
import "/scripts/test-animation.js";
import "/scripts/end-experiment.js";
import "/scripts/localiser.js";

import "/parameters.js";
import "/timeline.js";

///////////////////////
// INITIATE FIREBASE //
///////////////////////
const firebaseConfig = {
  apiKey: "AIzaSyBrNXaD1XrN1XaflA0vMYzoPhmsTY_DBLE",
  authDomain: "risky-decision-2019.firebaseapp.com",
  databaseURL: "https://risky-decision-2019.firebaseio.com",
  projectId: "risky-decision-2019",
  storageBucket: "risky-decision-2019.appspot.com",
  messagingSenderId: "961436604846",
  appId: "1:961436604846:web:3aba4d0cfa4d030407503f"
};

firebase.initializeApp(firebaseConfig); // Initialize Firebase

// Enable persistence
firebase
  .firestore()
  .enablePersistence()
  .catch(function(err) {
    if (err.code == "failed-precondition") {
      console.log(
        "Multiple tabs open, persistence can only be enabled in one tab at a a time."
      );
    } else if (err.code == "unimplemented") {
      console.log(
        "The current browser does not support all of the features required to enable persistence"
      );
    }
  });

// Set up firestore
firebase
  .auth()
  .signInAnonymously()
  .catch(function(error) {
    // Sign in
    // Handle Errors here.
    var errorCode = error.code;
    var errorMessage = error.message;
  });

var uid;
firebase.auth().onAuthStateChanged(function(user) {
  if (user) {

    var isAnonymous = user.isAnonymous;
    uid = user.uid;

    // check for pre-existing file
    var docRef = db.collection("iterations")
      .doc("pilot_v2.5")
      .collection("subjects")
      .doc(uid);

    docRef.get().then(function(doc) {
        if (doc.exists) {
          alert("A document for this user ID (" + uid + ") already exists. Open an incognito browser window, " +
                "or delete this document from Firestore.");
        }
    }).catch(function(error) {
        console.log("Error getting document:", error);
    });
  }
});

var subjectID = "";
if (window.location.search.indexOf("PROLIFIC_PID") > -1) {
  subjectID = getQueryVariable("PROLIFIC_PID");
} else {
  subjectID = Math.random()
    .toString()
    .slice(2, 8); // if no prolific ID, generate random ID (for testing)
}

var db = firebase.firestore(); // link to firestore database
var imagelist;
var sesstype;
var structype;
var parameters;
var timeline;


//////////////////////
// START EXPERIMENT //
//////////////////////

window.startExp = function() {

  // Experimenter input
  if (document.getElementById("btn-behav").checked == true) {
    sesstype = "behav";
  } else if (document.getElementById("btn-meg").checked == true) {
    sesstype = "meg";
  } else {
    alert("Session type has not been defined!");
    return;
  }

  if (document.getElementById("btn-sessA").checked == true) {
    structype = "A";
  } else if (document.getElementById("btn-sessB").checked == true) {
    structype = "B";
  } else {
    alert("Structure type has not been defined!");
    return;
  }

  imagelist = document.getElementById("input-imagelist").value.split(' ');

  if (document.getElementById("input-subjectID").value != "") {
    subjectID = document.getElementById("input-subjectID").value;
  }

  // Subject info
  console.log("Firebase User ID: " + uid);
  console.log("Subject ID: " + subjectID);
  if (window.location.search.indexOf("PROLIFIC_PID") > -1) {
    console.log("Prolific ID: " + getQueryVariable("PROLIFIC_PID"));
  }

  // Load experiment info
  parameters = loadParameters(); // load experiment parameters (e.g. trial timing, debug mode, etc.)

  if (sesstype == "meg") {
    parameters.exp_variables.meg_mode = true;
    parameters.exp_variables.questionnaires = false;
  } else {
    parameters.exp_variables.meg_mode = false;
    parameters.exp_variables.questionnaires = true;
  }
  parameters.exp_variables.structure_type = structype;

  if (imagelist[0]!="") {
    parameters.state_images.all_states = imagelist;
  }

  timeline = loadTimeline(parameters, db, uid);

  // Set up entry on Firestore
  db.collection("iterations")
    .doc("pilot_v2.5")
    .collection("subjects")
    .doc(uid)
    .set({
      subjectID: subjectID,
      date: new Date().toLocaleDateString(),
      time: new Date().toLocaleTimeString(),
      // structure: parameters.exp_variables.structureID,
      meg: parameters.exp_variables.meg_mode,
      structure: parameters.exp_variables.structure_type,
      path1_imgs: timeline[1][0],
      path2_imgs: timeline[1][1]
    });

  // Begin experiment
  jsPsych.init({
    timeline: timeline[0].flat(Infinity),
    preload_images: [
      parameters.state_images.all_states.map(
        x => "assets/img/states/" + x + parameters.state_images.file_ext
      ),
    ].flat(Infinity)
  });
};
