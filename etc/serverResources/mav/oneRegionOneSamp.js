$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var sampName = locSplit.pop();
	var regionName = locSplit.pop();
	//set up page
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		$.extend(names, {sample:sampName, regionName:regionName});
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		addPanelWithDiv("#mainContent", "snpsDiv", "View Snps and Sequence Info");
		addDiv("#mainContent", "sampTable");
		addDiv("#mainContent", "sampNameMenu");
		addDiv("#mainContent", "sampleChartMaster");
		addDiv("#mainContent", "mipOverlapGraphChart");
		//change title to current name
		$("title", "head").html(sampName);
		$("#jumboTitle").html(regionName + ":" + sampName);
		//d3.select("#snpsDiv")
		//	.append("a")
		//		.attr("href", "/" + rName + "/showOneGeneOneSampSnps/" + regionName + "/" + sampName)
		//		.text("View Snps and Sequence Info");
		//get sample names and the table with the sample names
		getJSON("/" + rName + "/getNamesForRegion/" + regionName).then(function(allMipNames){
			postJSON("/" + rName + "/getOneSampAllMipData/" + sampName, {"mipFamilies":allMipNames["mipFamilies"]}).then(function(mainSampInfoTab){
				//sample table
				var sampleTable =  new njhTable("#sampTable", mainSampInfoTab, sampName + "_sampInfo", true);	
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
				    var currentMipNames = _.intersection(allMipNames["mipFamilies"], allVals);
				    postJSON("/" + rName + "/getOneSampAllMipData/" + sampName, {"mipFamilies":currentMipNames}).then(function(mainSampInfoTab){
					    sampleChart.updateWithData(mainSampInfoTab);
					 	setTimeout(function(){
					 		sampleTable.updateWithData(mainSampInfoTab);
					 		$("#sampNameMenu").scrollView(60);
					 	}, 750);
				    }).catch(logRequestError).then(function(){
						//done loading sample table
				    	removeAllDivGifLoading();
					});
				};
				//create samp menu 
				var sampMenu = new njhCheckboxMenuOrganized("#sampNameMenu", allMipNames["mipFamilies"], updateChartOnClick);
				getJSON("/" + rName + "/mipOverlapGraphData/" + regionName + "/" + sampName).then(function(response) {
				  var mOverLapDrawer = new MipOverlapper("#mipOverlapGraphChart", response, $(window).width() -200);
				  		$(window).bind("resize", function(){
							mOverLapDrawer.updateSize();
						});
				}).catch(logRequestError).then(function(){
					//done loading sample table
				});
			}).catch(logRequestError).then(function(){
				//done loading sample table
			});
		}).catch(logRequestError).then(function(){
			//done loading mip fam names
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});
    	
    	
    	