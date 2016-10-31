$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var sampName = locSplit.pop();
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		//Set up Page
		addDiv("body", "topNav");
		$.extend(names, {sample:sampName});
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		setHeadTitle(sampName);
		$("#jumboTitle").html(sampName);
		addDiv("#mainContent", "sampTable");
		addDiv("#mainContent", "sampNameMenu");
		addDiv("#mainContent", "sampleChartMaster");
		getJSON("/" + rName + "/mipFamNamesForSamp/" + sampName).then(function(mipNames){
			postJSON("/" + rName + "/getOneSampAllMipData/" + sampName, {mipFamilies:mipNames["mipFamilies"]}).then(function(mainSampInfoTab){
				//get sample names and the table with the sample names
				//sample table
				var sampleTable =  new njhTable("#sampTable", mainSampInfoTab, sampName + "_sampInfo", false);	
				// bar graph for the sample info
				var sampleChart = new njhSampleChart("#sampleChartMaster", mainSampInfoTab, sampName +  "_sampChart","p_targetName", "c_barcodeFrac","c_clusterID", ["p_targetName", "h_popUID", "c_clusterID", "c_barcodeCnt", "c_barcodeFrac"]);
				//update the chart and table on click of the sample checkboxes
				function updateChartOnClick(columns) { 
					var gifLoading = prsentDivGifLoading();
					//get all currently checked sample boxes and then update the current samples  
				    var allVals = [];
				    $('#sampNameMenu :checked').each(function() {
				      allVals.push($(this).val());
				    });
				    var currentMipNames = _.intersection(mipNames["mipFamilies"], allVals);
				    postJSON("/" + rName + "/getOneSampAllMipData/" + sampName, {mipFamilies:currentMipNames}).then(function(mainSampInfoTab){
					 	sampleTable.updateWithData(mainSampInfoTab);
					 	sampleChart.updateWithData(mainSampInfoTab);
					 	$("#sampNameMenu").scrollView();
					}).catch(logRequestError).then(function(){
						//done loading 
						removeAllDivGifLoading();
					});
				};
				//create samp menu 
				var sampMenu = new njhCheckboxMenuOrganized("#sampNameMenu", mipNames["mipFamilies"], updateChartOnClick);
			}).catch(logRequestError).then(function(){
				//done loading sample table
			});
		}).catch(logRequestError).then(function(){
			//done loading mip names
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});