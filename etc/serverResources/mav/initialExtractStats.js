$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	//Set Up Page
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		//change title to current name
		setHeadTitle("Initial Read Stats");
		$("#jumboTitle").html("Initial Read Extraction Info");
		addDiv("#mainContent", "sampStatsTable");
		addDiv("#mainContent", "sampNameMenu");
		//get sample names and the table with the sample names
		var sampNames = names["samples"];
		postJSON("/" + rName + "/getInitialReadStats", {samples : sampNames}).then(function(mainSampStatsInfoTab){
			//create sample table 
			var sampleTable =  new njhTable("#sampStatsTable", mainSampStatsInfoTab, rName + "_allSampExtractInfo", true);	
			function updateChartOnClick() { 
				var gifLoading = prsentDivGifLoading();
				//get all currently checked sample boxes and then update the current samples  
			    var allVals = [];
			    $('#sampNameMenu :checked').each(function() {
			      allVals.push($(this).val());
			    });
			    var currentSampNames = _.intersection(sampNames, allVals);
			    postJSON("/" + rName + "/getInitialReadStats", {samples : currentSampNames}).then(function(mainSampStatsInfoTab){
				 	sampleTable.updateWithData(mainSampStatsInfoTab);
				 	$("#sampNameMenu").scrollView();
			    }).catch(logRequestError).then(function(){
					//done loading 
					removeAllDivGifLoading();
				});
			};
			//create samp menu 
			var sampMenu = new njhCheckboxMenu("#sampNameMenu", sampNames, updateChartOnClick);
		}).catch(logRequestError).then(function(){
			//done loading samples
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});