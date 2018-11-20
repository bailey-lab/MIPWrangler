$(document).ready(function(){
	//get current name from window location
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var mipName = locSplit.pop();
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		getJSON("/" + rName + "/getNamesForMipFam/" + mipName).then(function (namesForMipFam) {
			addDiv("body", "topNav");
			$.extend(names, {regionName:namesForMipFam["regionName"]})
			createNavBar("#topNav", names);
			addMainDiv("body", true);
			//change title to current name
			setHeadTitle(mipName);
			$("#jumboTitle").html(mipName);
			addDiv("#mainContent", "sampTable");
			addDiv("#mainContent", "sampNameMenu");
			addDiv("#mainContent", "sampleChartMaster");
			addDiv("#mainContent", "dnaViewer");
			addDiv("#mainContent", "popTable");
			postJSON("/" + rName + "/getOneMipAllSampsData/" + mipName, {"samples": namesForMipFam["samples"]}).then(function(mainPopInfoTab){
				//get sample names and the table with the sample names
				//sample table
				var sampleTable =  new njhTable("#sampTable", mainPopInfoTab, mipName + "_sampInfo", false);	
				// bar graph for the sample info
				var sampleChart = new njhSampleChart("#sampleChartMaster", mainPopInfoTab, mipName +  "_sampChart","s_Sample", "c_barcodeFrac","h_popUID", ["s_Sample", "h_popUID", "c_clusterID", "c_barcodeCnt", "c_barcodeFrac"]);
				var popUrls = ["/" + rName + "/getOneMipAllSampsPopData/" + mipName]
				popUrls.push("/" + rName + "/getOneMipPopSeqs/" + mipName);
				Promise.all(popUrls.map(function(popUrl){
					return postJSON(popUrl, {popUIDs:mainPopInfoTab["popUIDs"]});
				})).then(function(popData){
					//0 is popIno, 1 is popSeqs
					var popInfoTab = popData[0];
					var popSeqData = popData[1];
					//create seqViewer for the population final sequences 
					var sesUid = popSeqData["sessionUID"];
					var seqViewer = new njhSeqView("#dnaViewer", popSeqData);
					setUpCloseSession(sesUid);
					//create the population table and populate it 
					var popTable = new njhTable("#popTable", popInfoTab, mipName + "_popInfo", true);
					//update the chart and table on click of the sample checkboxes
					function updateChartOnClick() { 
						var gifLoading = prsentDivGifLoading();
						//get all currently checked sample boxes and then update the current samples  
					    var allVals = [];
					    $('#sampNameMenu :checked').each(function() {
					      allVals.push($(this).val());
					    });
					    var currentSampNames = _.intersection(namesForMipFam["samples"], allVals);
					    postJSON("/" + rName + "/getOneMipAllSampsData/" + mipName, {"samples": currentSampNames}).then(function(mainPopInfoTab){
						 	Promise.all(popUrls.map(function(popUrl){
								return postJSON(popUrl, {popUIDs:mainPopInfoTab["popUIDs"],  sessionUID:sesUid});
							})).then(function(popData){
								sampleChart.updateWithData(mainPopInfoTab);
							 	sampleTable.updateWithData(mainPopInfoTab);
								popTable.updateWithData(popData[0]);
							 	seqViewer.updateData(popData[1]);
						 		$("#sampNameMenu").scrollView(60, 0);
							}).catch(logRequestError).then(function(){
								//done loading pop  
							});
						}).catch(logRequestError).then(function(){
							//done loading update on click
							removeAllDivGifLoading();
						});
					};
					//create samp menu 
					var sampMenu = new njhCheckboxMenu("#sampNameMenu", namesForMipFam["samples"], updateChartOnClick);
				}).catch(logRequestError).then(function(){
					//done loading pop  
				});
			}).catch(logRequestError).then(function(){
				//done loading mip sample table 
			});
		}).catch(logRequestError).then(function(){
			//done loading names for mip fam
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});

