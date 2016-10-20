$(document).ready(function() {
	
	var rName = getRootName();
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		setHeadTitle("Mips!");
		$("#jumboTitle").html("Mips");
		addPanelWithDiv(".container.theme-showcase","regionLinks", "Mip Regions");
		addPanelWithDiv(".container.theme-showcase","sampLinks", "Samples");
		addPanelWithDiv(".container.theme-showcase","sampStatsLinks", "Initial Read Extraction Amounts Per Sample", "panel-info");
		addPanelWithDiv(".container.theme-showcase","mipStatsLinks", "Initial Read Extraction Amounts Per Mip Target", "panel-info");
		
		var cols = 10;
		var linkPre = "/" + rName + "/showRegionInfo/";
		var mouseOverC = "#999";
		var mouseLeaveC = "#FFF";
		var addTo = "#regionLinks";
		createLinksTable(addTo, linkPre, names["mipRegions"],cols, mouseOverC, mouseLeaveC);

		var linkPreSamp = "/" + rName + "/showOneSampAllMipData/";
		var addToSamp = "#sampLinks";
		createLinksTable(addToSamp, linkPreSamp, names["samples"],cols, mouseOverC, mouseLeaveC);
		
		$( "<a href=" + rName + "/showInitialReadStats>OverView of all Samples</a><p></p>" ).insertBefore( "#sampStatsLinks" );
		var linkPreSampStats = "/" + rName + "/showInitialReadStatsPerSample/";
		var addToSampStats = "#sampStatsLinks";
		createLinksTable(addToSampStats, linkPreSampStats, names["samples"],cols, mouseOverC, mouseLeaveC);

		var linkPreMipStats = "/" + rName + "/showInitialReadStatsPerMipTar/";
		var addToMipStats = "#mipStatsLinks";
		createLinksTable(addToMipStats, linkPreMipStats, names["mipTargets"],5, mouseOverC, mouseLeaveC);
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});

	

});