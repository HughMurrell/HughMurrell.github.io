// -----------------------------------------------------------

latest_data = [1, 1, 1, 3, 3, 7, 13, 17, 24, 38, 51, 62, 62, 116, 150, 202, 240, 274, 402, 554, 709, 927, 1170, 1187, 1280, 1326, 1353, 1380, 1462, 1505, 1585]

var valid_data = 1;
var running = 0;
var p_SI_default = 0.46;
var p_IR_default = 0.1;
var p_SI_ld_default = 0.08;
var p_IR_ld_default = 0.1;

var p_SI = p_SI_default;
var p_IR = p_IR_default;
var p_SI_ld = p_SI_ld_default;
var p_IR_ld = p_IR_ld_default;

var time_interval = 200;
var count = 0;
var N = 60000000;
var timeseries;


var sir_color = {D: "#000000", S: "#00ffff", I: "#f00000", R: "#00f000", C: "#0000f0" }

var epi_state = { S: (N-1)/N, I: 1/N, R: 0 };

function reset_params () {
    
    p_SI = p_SI_default;
    p_IR = p_IR_default;
    p_SI_ld = p_SI_ld_default;
    p_IR_ld = p_IR_ld_default;
    
    $("#p_SI").val(p_SI)
    $("#p_SI").keyup(update_p_SI);

    $("#p_IR").val(p_IR)
    $("#p_IR").keyup(update_p_IR);

    $("#p_SI_ld").val(p_SI_ld)
    $("#p_SI_ld").keyup(update_p_SI_ld);

    $("#p_IR_ld").val(p_IR_ld)
    $("#p_IR_ld").keyup(update_p_IR_ld);
}

function reset_history () {

  timeseries = {D: {label: "D", color: sir_color.D, data: []},
                S: {label: "S", color: sir_color.S, data: []},
                I: {label: "I", color: sir_color.I, data: []},
                R: {label: "R", color: sir_color.R, data: []},
                C: {label: "C", color: sir_color.C, data: []}
  };
    
  for (i = 0; i < latest_data.length; i++) {
        timeseries.D.data.push([i, latest_data[i]]);
  }

  timeseries.S.data.push([count, N*epi_state.S]);
  timeseries.I.data.push([count, N*epi_state.I]);
  timeseries.R.data.push([count, N*epi_state.R]);
  timeseries.C.data.push([count, N*(epi_state.R+epi_state.I)]);
}

reset_params();
reset_history();

var plotOptions = {
        lines: { show: true },
	    points: { show: true },
        xaxis: {min: 0},
        series: { shadowSize: 0 }
    };

var plot = $.plot($("#epicurves"), [], plotOptions);

update_plot();
update_counters();

$("#reset-button").click(reset_all);
$("#start-button").click(start_all);
$("#stop-button").click(stop_all);
$("#continue-button").click(continue_all);

setInterval(run_SIR, time_interval);


function run_SIR() {
    
    // console.log(Math.round(N*epi_state.I))

    if (running == 0)
       return;

    
    var s = epi_state.S;
    var i = epi_state.I;
    var r = epi_state.R;
    
    var beta = p_SI;
    var gamma = p_IR;
    
    if (count > 21) {
       beta = p_SI_ld;
       gamma = p_IR_ld;
    }
    
    epi_state.S += ( - beta * s * i );
    epi_state.I += ( + beta * s * i - gamma * i );
    epi_state.R += (gamma * i);
                    
    count++;

    timeseries.S.data.push([count, Math.round(N*epi_state.S)]);
    timeseries.I.data.push([count, Math.round(N*epi_state.I)]);
    timeseries.R.data.push([count, Math.round(N*epi_state.R)]);
    timeseries.C.data.push([count, N*(epi_state.R+epi_state.I)]);

    update_plot();

    update_counters();

    if (Math.round(N*epi_state.I) == 0)
       running = 0;
}

function update_plot () {
 plot.setData([timeseries.D, timeseries.I, timeseries.R, timeseries.C ]);
 plot.setupGrid();
 plot.draw();
}


function update_counters () {
  $("#count_S").html(Math.round(N*epi_state.S));
  $("#count_I").html(Math.round(N*epi_state.I));
  $("#count_R").html(Math.round(N*epi_state.R));
}

function update_p_SI () {
  p = Number($("#p_SI").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_SI").css("background-color", "#f88");
  } else {
     p_SI = p;
     valid_data = 1;
     $("#p_SI").css("background-color", "#fff");
  }
}

function update_p_SI_ld () {
  p = Number($("#p_SI_ld").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_SI_ld").css("background-color", "#f88");
  } else {
     p_SI_ld = p;
     valid_data = 1;
     $("#p_SI_ld").css("background-color", "#fff");
  }
}

function update_p_IR () {
  p = Number($("#p_IR").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_IR").css("background-color", "#f88");
  } else {
     p_IR = p;
     valid_data = 1;
     $("#p_IR").css("background-color", "#fff");
  }
}

function update_p_IR_ld () {
  p = Number($("#p_IR_ld").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_IR_ld").css("background-color", "#f88");
  } else {
     p_IR_ld = p;
     valid_data = 1;
     $("#p_IR_ld").css("background-color", "#fff");
  }
}

function reset_all () {
	running = 0;

    count = 0;
    epi_state = { S: (N-1)/N, I: 1/N, R: 0 };

    reset_params();
    reset_history();
    update_plot();

  update_counters();

   $("#start-text").fadeIn();
}

function start_all () {
    running = 1;
    count = 0;
    epi_state = { S: (N-1)/N, I: 1/N, R: 0 };
    reset_history();
    update_plot();
    update_counters();
}

function stop_all () {
    running = 0;
    update_plot();
    update_counters();
}

function continue_all () {
    running = 1;
    update_plot();
    update_counters();
}
