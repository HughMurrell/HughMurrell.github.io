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
var break_point = 20;
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
                         populations[ pop_data[i]["CountryCode"] ] = pop_data[i]["2018"];
      }
      defaultStartPlot();
    });
  }
)
.catch(function(err) {
  console.log('Fetch Population Error', err);
});

fetch('https://hughmurrell.github.io/CoVmodel/SIRou/data/stringencyindex.csv')
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


fetch('https://hughmurrell.github.io/CoVmodel/SIRou/data/confirmedcases.csv')
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
            countries.push(case_data[i]["CountryName"]+" "+case_data[i]["CountryCode"]);
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
       element.value = 'South Africa ZAF';
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
            if (case_data[i]["CountryName"] == country.substring(0,country.length-4)){
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
                N = populations[country.substring(country.length-3)];
                $("#p_N").prop('disabled', true);
                $('#p_N').val(N);
                p_R0 = (p_SI / p_IR).toFixed(2);
                $("#p_R0").prop('disabled', true);
                $('#p_R0').val(p_R0);
                p_Rn = (p_SI / p_IR).toFixed(2);
                $("#p_Rn").prop('disabled', true);
                $('#p_Rn').val(p_Rn);
            }
        }
                        timeseries.M.data = [];
        reset_all();
        console.log("Country = "+country+" with population "+N+" log2(N) = "+Math.log2(N))
    }
}


var valid_data = 1;
var running = 0;
var p_SI_default = 0.463;
var p_IR_default = 0.1;
var p_lambda_default = 0;
var p_eta_default = 1;
var p_mu_default = 0;
var p_eps_default = 1;
                        

var p_SI = p_SI_default;
var p_IR = p_IR_default;
                           
var p_lambda = p_lambda_default;
var p_eta = p_eta_default;
var p_mu = p_mu_default;
var p_eps = p_eps_default;

var p_R0;
var p_Rn;
        
var time_interval = 200;
var count = 0;
var sim_count = 0;
var N = 60000000;
var timeseries;

var sir_color = {D: "#000000", S: "#ffffff", I: "#ff0000", R: "#008800", C: "#000088" , L: "#f0f000", A: "#00ff00", M: "#00f0f0"}

var epi_state = { S: (N-1)/N, I: 1/N, R: 0 };
                        
var ameliorates = [];
                        


function reset_default_params () {
    
    p_SI = p_SI_default;
    p_IR = p_IR_default;
    p_lambda = p_lambda_default;
    p_eta = p_eta_default;
    p_mu = p_mu_default;
    
    reset_params();
}
                        
function reset_params () {
                            
    $("#p_SI").val(p_SI)
    $("#p_SI").keyup(update_p_SI);

    $("#p_IR").val(p_IR)
    $("#p_IR").keyup(update_p_IR);

    $("#p_lambda").val(p_lambda)
    $("#p_lambda").keyup(update_p_lambda);

    $("#p_eta").val(p_eta)
    $("#p_eta").keyup(update_p_eta);
    
    $("#p_mu").val(p_mu)
    $("#p_mu").keyup(update_p_mu);
                        
    p_R0 = (p_SI / p_IR).toFixed(2);
    $('#p_R0').val(p_R0);
    p_Rn = p_R0;
    $('#p_Rn').val(p_Rn);
                        
    recompute_ameliorates(p_mu);
                                                 
}

function reset_history () {
    
            var model_fit = [];
             
            if (timeseries != null) {
              model_fit = timeseries.M.data;
            }

   timeseries = {D: {label: "Data", color: sir_color.D, data: []},
                S: {label: "Susceptible", color: sir_color.S, data: []},
                I: {label: "Infective", color: sir_color.I, data: []},
                R: {label: "Recovered (or dead)", color: sir_color.R, data: []},
                C: {label: "Cumulative (Infective + Recovered)", color: sir_color.C, data: []},
                L: {label: "Stringency Index * max(Data)", color: sir_color.L, data: []},
                A: {label: "Amelioration Factor * max(Data)", color: sir_color.A, data: []},
                M: {label: "Last Optimal Model fit", color: sir_color.M, data: model_fit}
  };
    
  var max_case = Math.max(...cases);
  for (var i = 0; i < cases.length; i++) {
        timeseries.D.data.push([i, cases[i]]);
        timeseries.L.data.push([i, indicies[i] * max_case]);
        timeseries.A.data.push([i, ameliorates[i] * max_case]);
  }
                        
  count = 0;
  epi_state = { S: (N-cases[0])/N, I: cases[0]/N, R: 0 };
    
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
$("#opt-beta-button").click(optimize_beta);
$("#opt-lambda-button").click(optimize_lambda);
$("#opt-beta-lambda-button").click(optimize_beta_lambda);

setInterval(run_SIR, time_interval);

function run_SIR() {
    
    if (running == 0)
       return;

    count++;

    var s = epi_state.S;
    var i = epi_state.I;
    var r = epi_state.R;
    
    var beta = p_SI;
    var gamma = p_IR;
                        
    var prev_count = Math.abs(Math.min(count, indicies.length-1))
    var f = ameliorates[prev_count];
    
    epi_state.S += ( - f * beta * s * i );
    epi_state.I += ( + f * beta * s * i - gamma * i );
    epi_state.R += (gamma * i);
                    
    timeseries.S.data.push([count, Math.round(N*epi_state.S)]);
    timeseries.I.data.push([count, Math.round(N*epi_state.I)]);
    timeseries.R.data.push([count, Math.round(N*epi_state.R)]);
    timeseries.C.data.push([count, N*(epi_state.R+epi_state.I)]);
            
    update_p_Rn(f);

    update_plot();

    /*
            if (count >= cases.length - 1) {    // (Math.round(N*epi_state.I) == 0) {
              console.log("loss = "+sim_loss(p_SI, p_IR, p_lambda, p_eta, p_mu, p_eps));
              running = 0;
            }
     */
       
}

                           
function update_plot () {
 // console.log(timeseries.M);

 plot.setData([timeseries.D, timeseries.I, timeseries.R,
               timeseries.L, timeseries.A, timeseries.M, timeseries.C ]);
 plot.setupGrid();
 max_y = plot.getOptions().yaxes[0].max;
 plot.draw();
}
                           
function update_p_R0 () {
    p_R0 = (p_SI / p_IR).toFixed(2);
    $('#p_R0').val(p_R0);
}
            
function update_p_Rn (lambda) {
    p_Rn = (lambda * p_SI / p_IR).toFixed(2);
    $('#p_Rn').val(p_Rn);
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
        reset_params();
        reset_history();
        update_plot();
        valid_data = 1;
        $("#p_lambda").css("background-color", "#fff");
    }
}
            
            function update_p_eta () {
               p = Number($("#p_eta").val());
               if (isNaN(p) || p<1 || p>16) {
                  valid_data = 0;
                  $("#p_eta").css("background-color", "#f88");
               } else {
                  p_eta = p;
                  reset_params();
                  reset_history();
                  update_plot();
                  valid_data = 1;
                  $("#p_eta").css("background-color", "#fff");
               }
             }


function update_p_mu () {
  p = Number($("#p_mu").val());
  if (isNaN(p) || p<0.0 || p>10.0) {
     valid_data = 0;
     $("#p_mu").css("background-color", "#f88");
  } else {
     p_mu = p;
     reset_params();
     reset_history();
     update_plot();
     valid_data = 1;
     $("#p_mu").css("background-color", "#fff");
  }
}
            
 
function reset_all () {
	running = 0;
    count = 0;

    reset_default_params();
    reset_history();
    update_plot();
    
}

function start_all () {
    count = 0;
    epi_state = { S: (N-1)/N, I: 1/N, R: 0 };
    reset_history();
    running = 1;
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
            
function fractional_part(xv) {
            var x = Math.abs(xv)
            if (x>1) {
                return (x % 1);
            } else {
                return x;
            }
}
            
function sigmoid(t,lambda,eta) {
                var r = 1 - lambda * Math.pow( t, eta ); //  / (1+Math.pow(Math.E, - eta * (t-mu) )) ;
                if (isNaN(r)) {
                    r = 1;
                }
                return r;
}
                                    
function recompute_ameliorates(lag) {
                ameliorates = [];
                for (var i = 0; i < cases.length; i++) {
                    ind = Math.max(0,i-lag)
                    ameliorates.push(sigmoid( indicies[ind], p_lambda, p_eta));
                }
}
            
function sim_loss(beta, gamma, lambda, eta, mu, eps) {
    reset_history();
    timeseries.M.data = [];
    var state = { S: (N-1)/N, I: 1/N, R: 0, C: 1/N };
    timeseries.M.data.push([0, state.C]);
    var loss = Math.abs(cases[0] - N * state.C );
    var s;
    var i;
    var f;
    var j=1;
            //indicies[j] < mu   j<break_point &
    while (j<cases.length) {
        var ind = Math.max(0,j-mu)   // allowing for lag
        f = sigmoid( indicies[ind], lambda, eta);
        s = state.S;
        i = state.I;
        state.S += ( - f * beta * s * i );
        state.I += ( + f * beta * s * i - gamma * i );
        state.R += (gamma * i);
        state.C = state.I + state.R;
        // console.log("cases = "+cases[j]+" state = "+(N*state.C));
        // if ( Math.abs(cases[j] - N*state.C ) > (eps* Math.log2(N)) ) {
        //               loss += 1;
        // }
        loss += Math.pow((Math.log(1+cases[j]) - Math.log(1+N * state.C)),2 );
        if ( N*state.C < cases[cases.length-1] ) {
            timeseries.M.data.push([j, N*state.C]);
        }
        j++;
    }
    return loss
}
            
function simulate_loss_beta(x){
    sim_count += 1;
    var beta = x;
    var gamma = p_IR;
    var lambda = p_lambda;
    var eta = p_eta;
    var mu = p_mu;
    var eps = p_eps;
    var loss = sim_loss(beta, gamma, lambda, eta, mu, eps);
    // console.log("loss = "+loss+" at "+x+" after call "+sim_count);
    return loss;
}

function optimize_beta () {
    running = 0;
    console.log("optimizing beta ....");
    
    function loss(x) {
        return simulate_loss_beta(x);
    }
    sim_count = 0;
    var best = 0.01;
    var guess = 0.01;
    var best_loss = loss(guess);
    var curr_loss = best_loss;
    var delta = 0.01;
            
            while ( guess < 1.0 ) {
              guess += delta;
              curr_loss = loss(guess);
              if ( curr_loss < best_loss ) {
                best = guess;
                best_loss = curr_loss;
                update_plot();
              }
            }
            
            var crude_best = best;
            
            guess = best - delta;
            best = best - delta;
            var best_loss = loss(guess);
            var curr_loss = best_loss;
            var fine_delta = delta/10;
            
            while ( guess < crude_best + delta ) {
              guess += fine_delta;
              curr_loss = loss(guess);
              if ( curr_loss < best_loss ) {
                best = guess;
                best_loss = curr_loss;
              }
            }
    
    console.log("loss = "+loss(best)+" at "+best+" after "+sim_count+" function calls");
    p_SI=best; // p_lambda=x[1]; //  p_mu=x[2];  // p_IR=x[1];
    reset_params();
    reset_history();
    update_plot();
}
                                     
function simulate_loss_lambda(x){
    sim_count += 1;
    var beta = p_SI;
    var gamma = p_IR;
    var lambda = x;
    var eta = p_eta;
    var mu = p_mu;
    var eps = p_eps;
                                                     
    var loss = sim_loss(beta, gamma, lambda, eta, mu, eps)
    // console.log("loss = "+loss+" at "+x+" after call "+sim_count);

    return loss;
                                                     
}

                                     function optimize_lambda () {
                                         running = 0;
                                         console.log("optimizing lambda ....");
                                         
                                         function loss(x) {
                                             return simulate_loss_lambda(x);
                                         }
                                         sim_count = 0;
                                         var best = 0.01;
                                         var guess = 0.01;
                                         var best_loss = loss(guess);
                                         var curr_loss = best_loss;
                                         var delta = 0.001;
                                                 
                                                 while ( guess <= 1.0-delta ) {
                                                   guess += delta;
                                                   curr_loss = loss(guess);
                                                   if ( curr_loss < best_loss ) {
                                                     best = guess;
                                                     best_loss = curr_loss;
                                                   }
                                                 }
                                                 
                                                 var crude_best = best;
                                                 
                                                 guess = best - delta;
                                                 best = best - delta;
                                                 var best_loss = loss(guess);
                                                 var curr_loss = best_loss;
                                                 var fine_delta = delta/10;
                                                 
                                                 while ( guess < crude_best + delta - fine_delta) {
                                                   guess += fine_delta;
                                                   curr_loss = loss(guess);
                                                   if ( curr_loss < best_loss ) {
                                                     best = guess;
                                                     best_loss = curr_loss;
                                                   }
                                                 }
                                         
                                         console.log("loss = "+loss(best)+" at "+best+" after "+sim_count+" function calls");
                                         p_lambda=best; //  p_mu=x[2];  // p_IR=x[1]; // p_SI=best; //
                                         reset_params();
                                         reset_history();
                                         update_plot();
                                     }

function simulate_loss_beta_lambda(x_beta, x_lambda){
    sim_count += 1;
    var beta = x_beta;
    var gamma = p_IR;
    var lambda = x_lambda;
    var eta = p_eta;
    var mu = p_mu;
    var eps = p_eps;
                                                     
    var loss = sim_loss(beta, gamma, lambda, eta, mu, eps)
    // console.log("loss = "+loss+" at "+x+" after call "+sim_count);

    return loss;
                                                     
}

function optimize_beta_lambda () {
        running = 0;
        console.log("optimizing beta and lambda ....");
                                         
        function loss(x_beta,x_lambda) {
            return simulate_loss_beta_lambda(x_beta, x_lambda);
        }
            
        sim_count = 0;
        var delta = 0.001;
        var best_loss = loss(0,0);
        var curr_loss = best_loss;
        var best_beta = 0;
        var best_lambda = 0;
                                         

            
        for (var x_beta = 0; x_beta <= 1; x_beta += delta) {
            for (var x_lambda = 0; x_lambda <= 1; x_lambda+=delta ) {
                curr_loss = loss(x_beta, x_lambda);
                if (curr_loss < best_loss) {
                    best_loss = curr_loss;
                    best_beta = x_beta;
                    best_lambda = x_lambda;
                }
            }
        }
            
        console.log("loss = "+loss(best_beta,best_lambda)+" at "+best_beta+" , "+best_lambda+
                    " after "+sim_count+" function calls");
        p_SI=best_beta;
        p_lambda=best_lambda;
        reset_params();
        reset_history();
        update_plot();
}
