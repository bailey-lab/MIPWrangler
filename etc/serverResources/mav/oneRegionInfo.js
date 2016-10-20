$(document).ready(function() {
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var regionName = locSplit.pop();
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		$.extend(names, {regionName:regionName})
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		setHeadTitle(regionName);
		$("#jumboTitle").html("Mip Targets for " + regionName);
		addPanelWithDiv("#mainContent","mipLinks", regionName);
		addPanelWithDiv("#mainContent","sampLinks", "Samples");
		getJSON("/" + rName + "/getNamesForRegion/" + regionName).then(function (names) {	
			var cols = 10;
			var mouseOverC = "#999";
			var mouseLeaveC = "#FFF";
			//add links for the mip families
			var linkPre = "/" + rName + "/showMipInfo/";
			var addTo = "#mipLinks";
			createLinksTable(addTo, linkPre, names["mipFamilies"], cols, mouseOverC, mouseLeaveC);
			//add links for the samples
			var linkPreSamp = "/" + rName + "/showRegionInfoForSamp/" + regionName + "/";
			var addToSamp = "#sampLinks";
			createLinksTable(addToSamp, linkPreSamp, names["samples"], cols, mouseOverC, mouseLeaveC);
		}).catch(logRequestError).then(function(){
			//done loading names for region
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
	
});