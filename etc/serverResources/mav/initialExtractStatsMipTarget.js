$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var mipTar = locSplit.pop();
	//Set Up Page
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		//change title to current name
		setHeadTitle(mipTar + " Stats");
		$("#jumboTitle").html(mipTar + " Extraction Info");
		addDiv("#mainContent", "mipTarStatsTable");
		addDiv("#mainContent", "samplesNameMenu");
		//get sample names and the table with the sample names
		getJSON("/" + rName + "/samplesForExtractedMip/" + mipTar).then(function(allSamples){
			postJSON("/" + rName + "/getInitialReadStatsPerMipTar/" + mipTar,{samples:allSamples["samples"]} ).then(function(mainMipTarStatsInfoTab){
				//sample table
				var mipTarTable =  new njhTable("#mipTarStatsTable", mainMipTarStatsInfoTab,  mipTar + "_mipTarExtractInfo", true);	
				function updateChartOnClick() { 
					var gifLoading = prsentDivGifLoading();
					//get all currently checked sample boxes and then update the current samples  
				    var allVals = [];
				    $('#samplesNameMenu :checked').each(function() {
				      allVals.push($(this).val());
				    });
				    var currentSamples = _.intersection(allSamples["samples"], allVals);
				    postJSON("/" + rName + "/getInitialReadStatsPerMipTar/" + mipTar,{samples:currentSamples} ).then(function(mainMipTarStatsInfoTab){
				    	mipTarTable.updateWithData(mainMipTarStatsInfoTab);
					    $("#sampNameMenu").scrollView();
				    }).catch(logRequestError).then(function(){
						//done loading 
						removeAllDivGifLoading();
					});
				};
				var sampMenu = new njhCheckboxMenu("#samplesNameMenu", allSamples["samples"], updateChartOnClick);
			}).catch(logRequestError).then(function(){
				//done loading table
			});
		}).catch(logRequestError).then(function(){
			//done loading samples
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});

