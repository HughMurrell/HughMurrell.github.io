// -----------------------------------------------------------
// get the data
console.log("started")


var csv;
var latest_data = [];
var countries = [];
var country;
var data_read = false;
var break_point;
var series_point;


fetch('https://hughmurrell.github.io/CoVmodel/SIRou/data/confirmedcases-Table 1.csv')
.then(
  function(response) {
    if (response.status !== 200) {
      console.log('Looks like there was a problem. Status Code: ' +
        response.status);
      return;
    }

    // Examine the text in the response
    response.text().then(function(data) {
      // console.log(data);
      csv = data;
      jh_data = d3.csvParse(csv);
      console.log(jh_data[0]);

      for (var i = 0; i < jh_data.length; i++){
            countries.push(jh_data[i]["CountryName"]);
      }
      countries.sort();
      populateCountries("country2");
                         
    // set default country and default intervention break point
    var element = document.getElementById('country2');
    element.value = 'South Africa';
    var event = new Event('change');
    element.dispatchEvent(event);

    update_plot();

      data_read = true;
    });
  }
)
.catch(function(err) {
  console.log('Fetch Error :-S', err);
});


function populateCountries(countryElementId){
    // given the id of the <select> tag as function argument, it inserts <option> tags
    var countryElement = document.getElementById(countryElementId);
    countryElement.length=0;
    // countryElement.options[0] = new Option('Select Country','-1');
    countryElement.selectedIndex = 0;
    for (var i=0; i<countries.length; i++) {
        countryElement.options[countryElement.length] = new Option(countries[i],countries[i]);
    }
    // Assigned all countries. Now assign event listener.

    countryElement.onchange = function(){
        country = countries[countryElement.selectedIndex];
        console.log(country);
        for (var i=0; i<jh_data.length; i++){
            if (jh_data[i]["CountryName"] == country){
                values = Object.values(jh_data[i]);
                // the 9 below is to drop non-data fields
                int_values = values.slice(2,values.length).map(Number);
                index=int_values.findIndex(function(number) {
                  return number > 0;
                });
                latest_data = int_values.slice(index,int_values.length);
                // break_point = jh_data[i]["intervention"] - index;
                // p_SI = jh_data[i]["beta_before"];
                // p_SI_ld = jh_data[i]["beta_after"];
                // p_IR = jh_data[i]["gamma_before"];
                // p_IR_ld = jh_data[i]["gamma_after"];
                // N = jh_data[i]["population"];
                $("#p_N").prop('disabled', true);
                $('#p_N').val(N);
                p_R0 = (p_SI / p_IR).toFixed(1);
                p_R0_ld = (p_SI_ld / p_IR_ld).toFixed(1);
                $("#p_R0").prop('disabled', true);
                $('#p_R0').val(p_R0);
                $("#p_R0_ld").prop('disabled', true);
                $('#p_R0_ld').val(p_R0_ld);
            }
        }
        reset_all();
        // console.log(country)
    }
}


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
                                     
var p_R0;
var p_R0_ld;
        
var time_interval = 200;
var count = 0;
var N = 60000000;
var timeseries;


var sir_color = {D: "#000000", S: "#00ffff", I: "#f00000", R: "#00f000", C: "#0000f0" }

var epi_state = { S: (N-1)/N, I: 1/N, R: 0 };

function reset_params () {
    
    // p_SI = p_SI_default;
    // p_IR = p_IR_default;
    // p_SI_ld = p_SI_ld_default;
    // p_IR_ld = p_IR_ld_default;
    
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
    
   timeseries = {D: {label: "Data", color: sir_color.D, data: []},
                S: {label: "Susceptible", color: sir_color.S, data: []},
                I: {label: "Infective", color: sir_color.I, data: []},
                R: {label: "Recovered (or dead)", color: sir_color.R, data: []},
                C: {label: "Cumulative (Infective + Recovered)", color: sir_color.C, data: []},
                B: {label: "Start of Intervention", color: sir_color.B, data: [] }
  };
    
  for (i = 0; i < latest_data.length; i++) {
        timeseries.D.data.push([i, latest_data[i]]);
  }
    
  timeseries.B.data.push([break_point, 0]);
  timeseries.B.data.push([break_point, latest_data[latest_data.length-1]]);

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
        series: { shadowSize: 0 },
        grid: {hoverable: true, clickable:true},
         legend:{
                   // backgroundOpacity: 0.5,
                   // noColumns: 0,
                   // backgroundColor: "green",
                   position: "nw"
               }
    };

var placeholder = $("#epicurves");
var plot = $.plot($("#epicurves"), [], plotOptions);

$("#epicurves").on("plotclick",function(event,pos,item){
    if(item){
        break_point = item.series.data[item.dataIndex][0];
        console.log(break_point);
        reset_all();
    }
});


$("#reset-button").click(reset_all);
$("#start-button").click(start_all);
$("#stop-button").click(stop_all);
$("#continue-button").click(continue_all);

setInterval(run_SIR, time_interval);


function run_SIR() {
    
    if (running == 0)
       return;

    
    var s = epi_state.S;
    var i = epi_state.I;
    var r = epi_state.R;
    
    var beta = p_SI;
    var gamma = p_IR;
    
    if (count >= break_point) {
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

    if (Math.round(N*epi_state.I) == 0)
       running = 0;
}

function update_plot () {
 plot.setData([timeseries.D, timeseries.I, timeseries.R,
               timeseries.C, timeseries.B ]);
 plot.setupGrid();
 plot.draw();
}
                           
function update_p_R0 () {
    p_R0 = (p_SI / p_IR).toFixed(1);
    $('#p_R0').val(p_R0);
}
                                      
function update_p_R0_ld () {
   p_R0_ld = (p_SI_ld / p_IR_ld).toFixed(1);
   $('#p_R0_ld').val(p_R0_ld);
}

function update_p_SI () {
  p = Number($("#p_SI").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_SI").css("background-color", "#f88");
  } else {
     p_SI = p;
     valid_data = 1;
     update_p_R0();
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
     update_p_R0_ld();
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
     update_p_R0();
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
     update_p_R0_ld();
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
    
}

function start_all () {
    running = 1;
    count = 0;
    epi_state = { S: (N-1)/N, I: 1/N, R: 0 };
    reset_history();
    update_plot();
}

function stop_all () {
    running = 0;
    update_plot();
}

function continue_all () {
    running = 1;
    update_plot();
}
