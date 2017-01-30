$(document).ready(function() {
	
	var rName = getRootName();
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		//createNavBar("#topNav", names);
		addMainDiv("body", true);
		setHeadTitle("Mips!");
		$("#jumboTitle").html("Mips on Genomes");
		addPanelWithDiv(".container.theme-showcase","mipsLinks", "Mip");
		//addPanelWithDiv(".container.theme-showcase","sampLinks", "Ge");

		var cols = 10;
		var linkPre = "/" + rName + "/showMipSeqs/";
		var mouseOverC = "#999";
		var mouseLeaveC = "#FFF";
		var addTo = "#mipsLinks";
		createLinksTable(addTo, linkPre, names["mips"], cols, mouseOverC, mouseLeaveC);
		d34.select("#mipsLinks").select("table").style("overflow", "scroll")
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});

	

});