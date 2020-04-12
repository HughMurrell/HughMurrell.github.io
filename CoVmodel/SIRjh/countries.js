
var country_arr = new Array("Afghanistan", "Albania");


function populateCountries(countryElementId){
	// given the id of the <select> tag as function argument, it inserts <option> tags
	var countryElement = document.getElementById(countryElementId);
	countryElement.length=0;
	countryElement.options[0] = new Option('Select Country','-1');
	countryElement.selectedIndex = 0;
	for (var i=0; i<country_arr.length; i++) {
		countryElement.options[countryElement.length] = new Option(country_arr[i],country_arr[i]);
	}
    // Assigned all countries. Now assign event listener.

    countryElement.onchange = function(){
        console.log( country_arr[countryElement.selectedIndex-1] );
}


};

