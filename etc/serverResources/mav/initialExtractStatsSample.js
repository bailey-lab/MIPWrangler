$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var sample = locSplit.pop();
	//Set Up Page
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		$.extend(names, {sample:sample});
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		//change title to current name
		setHeadTitle(sample + " Stats");
		$("#jumboTitle").html(sample + " Extraction Info");
		addDiv("#mainContent", "sampleStatsTable");
		addDiv("#mainContent", "targetNameMenu");
		//get sample names and the table with the sample names
		getJSON("/" + rName + "/extractedMipsForSample/" + sample).then(function(allMipTargets){
			postJSON("/" + rName + "/getInitialReadStatsPerSample/" + sample,{mipTargets:allMipTargets["mipTargets"]} ).then(function(mainsampleStatsInfoTab){
				//sample table
				var sampleTable =  new njhTable("#sampleStatsTable", mainsampleStatsInfoTab,  sample + "_sampleExtractInfo", true);	
				function updateChartOnClick() { 
					var gifLoading = prsentDivGifLoading();
					//get all currently checked sample boxes and then update the current samples  
				    var allVals = [];
				    $('#targetNameMenu :checked').each(function() {
				      allVals.push($(this).val());
				    });
				    var currentTargets = _.intersection(allMipTargets["mipTargets"], allVals);
				    postJSON("/" + rName + "/getInitialReadStatsPerSample/" + sample,{mipTargets:currentTargets} ).then(function(mainsampleStatsInfoTab){
				    	sampleTable.updateWithData(mainsampleStatsInfoTab);
					    $("#sampNameMenu").scrollView();
				    }).catch(logRequestError).then(function(){
						//done loading 
						removeAllDivGifLoading();
					});
				};
				var targetMenu = new njhCheckboxMenu("#targetNameMenu", allMipTargets["mipTargets"], updateChartOnClick);
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