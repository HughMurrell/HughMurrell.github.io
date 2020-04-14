// -----------------------------------------------------------
// get the data
console.log("started")


var countries = [];
var populations ={};
var country;
var case_data = [];
var cases = [];
var indicies_data = [];
var indicies = [];
var data_read = false;
var break_point;
var series_point;


fetch('https://hughmurrell.github.io/CoVmodel/SIRou/data/wb_population.csv')
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
      pop_data = d3.csvParse(data);

      for (var i = 0; i < pop_data.length; i++){
                         populations[ pop_data[i]["CountryName"] ] = pop_data[i]["2018"];
      }
      defaultStartPlot();
    });
  }
)
.catch(function(err) {
  console.log('Fetch Population Error', err);
});

fetch('https://hughmurrell.github.io/CoVmodel/SIRou/data/stringencyindex-Table 1.csv')
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
      indicies_data = d3.csvParse(csv);
      defaultStartPlot();
    });
  }
)
.catch(function(err) {
  console.log('Fetch Cases Error', err);
});


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
      case_data = d3.csvParse(csv);

      for (var i = 0; i < case_data.length; i++){
            countries.push(case_data[i]["CountryName"]);
      }
      countries.sort();
      populateCountries("country2");
                         
      defaultStartPlot();

    });
  }
)
.catch(function(err) {
  console.log('Fetch Cases Error', err);
});

function defaultStartPlot() {
    // set default country and default intervention break point
       var element = document.getElementById('country2');
       element.value = 'South Africa';
       var event = new Event('change');
       element.dispatchEvent(event);

       update_plot();
}





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
        for (var i=0; i<case_data.length; i++){
            if (case_data[i]["CountryName"] == country){
                values = Object.values(case_data[i]);
                keys = Object.keys(case_data[i]);
                // the 2 below is to drop non-data fields
                int_values = values.slice(2,values.length).map(Number);
                index=int_values.findIndex(function(number) {
                  return number > 0;
                });
                var first_case_date = keys[index+2];
                console.log(first_case_date);
                cases = int_values.slice(index,int_values.length);
                if (typeof indicies_data[i] != "undefined") {
                    values = Object.values(indicies_data[i]);
                    // replace missing values with previous value
                    for (var j=0; j<values.length-1; j++){
                        if ( values[j+1] == "" ) {
                            values[j+1]=values[j];
                        }
                    }
                    keys = Object.keys(indicies_data[i]);
                    int_values = values.slice(2,values.length).map(Number);
                    indicies = int_values.slice(index,int_values.length);
                    for (var j=0; j<indicies.length; j++){
                        indicies[j] = (indicies[j] ) / 100.0;
                    }
                }
                
                // break_point = case_data[i]["intervention"] - index;
                // p_SI = case_data[i]["beta_before"];
                // p_SI_ld = case_data[i]["beta_after"];
                // p_IR = case_data[i]["gamma_before"];
                // p_IR_ld = case_data[i]["gamma_after"];
                N = populations[country];
                $("#p_N").prop('disabled', true);
                $('#p_N').val(N);
                p_R0 = (p_SI / p_IR).toFixed(1);
                $("#p_R0").prop('disabled', true);
                $('#p_R0').val(p_R0);
            }
        }
        reset_all();
        // console.log(country)
    }
}


var valid_data = 1;
var running = 0;
var p_SI_default = 0.463;
var p_IR_default = 0.1;
var p_lambda_default = 0.87;
var p_mu_default = 0.95;
var p_eta_default = 1.0;

var p_SI = p_SI_default;
var p_IR = p_IR_default;
                           
var p_lambda = p_lambda_default;
var p_mu = p_mu_default;
var p_eta = p_eta_default;
                                     
var p_R0;
        
var time_interval = 200;
var count = 0;
var N = 60000000;
var timeseries;

var sir_color = {D: "#000000", S: "#00ffff", I: "#f00000", R: "#00f000", C: "#0000f0" ,L: "#ffff00"}

var epi_state = { S: (N-1)/N, I: 1/N, R: 0 };

function reset_params () {
    
    p_SI = p_SI_default;
    p_IR = p_IR_default;
    p_lambda = p_lambda_default;
    p_mu = p_mu_default;
    p_eta = p_eta_default;
    
    $("#p_SI").val(p_SI)
    $("#p_SI").keyup(update_p_SI);

    $("#p_IR").val(p_IR)
    $("#p_IR").keyup(update_p_IR);

    $("#p_lambda").val(p_lambda)
    $("#p_lambda").keyup(update_p_lambda);

    $("#p_mu").val(p_mu)
    $("#p_mu").keyup(update_p_mu);
                        
    $("#p_eta").val(p_eta)
    $("#p_eta").keyup(update_p_eta);
}

function reset_history () {
    
   timeseries = {D: {label: "Data", color: sir_color.D, data: []},
                S: {label: "Susceptible", color: sir_color.S, data: []},
                I: {label: "Infective", color: sir_color.I, data: []},
                R: {label: "Recovered (or dead)", color: sir_color.R, data: []},
                C: {label: "Cumulative (Infective + Recovered)", color: sir_color.C, data: []},
                L: {label: "Stringency Index * max(Data)", color: sir_color.L, data: [] }
  };
    
  var max_case = Math.max(...cases);
  console.log(max_case);
  for (i = 0; i < cases.length; i++) {
        timeseries.D.data.push([i, cases[i]]);
        timeseries.L.data.push([i, indicies[i] * max_case]);
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
    var lambda = p_lambda;
    var mu = p_mu;
    var eta = p_eta;
    
    // function sigmoid(x,x0,b) {
    //     return Math.pow(  ( 1/(1+Math.pow(Math.E, -b*(x-(1-x0))) ) * x ), 1-x ) ;
    // }
                           
    function sigmoid(x,lambda,mu,eta) {
        var k = 400*eta;
        return 1 - lambda / (1+Math.pow(Math.E, - k * (x-mu) )) ;
    }
    
                        /*
    function sigmoid(x,lambda,mu) {
                        if (x > mu) {
                          return 1 - lambda * x;
                        } else {
                          return 1
                        }
    }
                        */
                        
    var max_count = Math.abs(Math.min(count-1, indicies.length-1))
    var f = sigmoid( indicies[max_count], lambda, mu, eta)
                        
                        console.log(f);
    
    epi_state.S += ( - f * beta * s * i );
    epi_state.I += ( + f * beta * s * i - gamma * i );
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
               timeseries.C, timeseries.L ]);
 plot.setupGrid();
 max_y = plot.getOptions().yaxes[0].max;
 plot.draw();
}
                           
function update_p_R0 () {
    p_R0 = (p_SI / p_IR).toFixed(1);
    $('#p_R0').val(p_R0);
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

function update_p_lambda () {
    p = Number($("#p_lambda").val());
    if (isNaN(p) || p<0.0 || p>1.0) {
        valid_data = 0;
    $("#p_lambda").css("background-color", "#f88");
    } else {
        p_lambda = p;
        valid_data = 1;
        $("#p_lambda").css("background-color", "#fff");
    }
}

function update_p_mu () {
  p = Number($("#p_mu").val());
  if (isNaN(p) || p<0.0 || p>1.0) {
     valid_data = 0;
     $("#p_mu").css("background-color", "#f88");
  } else {
     p_mu = p;
     valid_data = 1;
     $("#p_mu").css("background-color", "#fff");
  }
}
            
function update_p_eta () {
    p = Number($("#p_eta").val());
    if (isNaN(p) || p<0.0 || p>1.0) {
        valid_data = 0;
        $("#p_eta").css("background-color", "#f88");
    } else {
        p_eta = p;
        valid_data = 1;
        $("#p_eta").css("background-color", "#fff");
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
