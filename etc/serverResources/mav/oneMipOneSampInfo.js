$(document).ready(function(){
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var sampName = locSplit.pop();
	var mipName = locSplit.pop();
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		getJSON("/" + rName + "/getNamesForMipFam/" + mipName).then(function (namesForMipFam) {
			var sampNames = namesForMipFam["samples"]
			//Set up Page
			addDiv("body", "topNav");
			$.extend(names, {regionName:namesForMipFam["regionName"], sample:sampName})
			createNavBar("#topNav", names);
			addMainDiv("body", true);
			setHeadTitle(mipName);
			$("#jumboTitle").html(mipName + ":" + sampName);
			addH1 ("#mainContent", "Final Clusters");
			addDiv("#mainContent", "finalSeqs");
			addDiv("#mainContent", "popTable");
			postJSON("/" + rName + "/getMipOneSampFinalSeqs/" + mipName + "/" + sampName).then(function(mainFinalData){
				getJSON("/" + rName + "/getOneMipOneSampsData/" + mipName + "/" + sampName).then(function(mainPopInfoTab){
					var seqViewerFinal = new njhSeqView("#finalSeqs", mainFinalData);
					var popTable = new njhTable("#popTable", mainPopInfoTab, mipName + "_" + sampName + "_popInfo", true);
				}).catch(logRequestError).then(function(){
					//done loading names for mip samp info
				});
			}).catch(logRequestError).then(function(){
				//done loading names for mip final seqs
			});
		}).catch(logRequestError).then(function(){
			//done loading names for mip fam
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});
    	
