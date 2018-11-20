$(document).ready(function(){
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = locSplit[2];
	var mipName = locSplit.pop();
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		getJSON("/" + rName + "/getNamesForMipFam/" + mipName).then(function (namesForMipFam) {
			var sampNames = namesForMipFam["samples"]
			//Set up Page
			addDiv("body", "topNav");
			$.extend(names, {regionName:namesForMipFam["regionName"]})
			createNavBar("#topNav", names);
			addMainDiv("body", true);
			setHeadTitle(mipName);
			$("#jumboTitle").html(mipName);
			addPanelWithDiv("#mainContent","linkPanelBody", "Data Across all Samples");
			addPanelWithDiv("#mainContent","sampLinks", "Per Sample");
		
			$('<a>',{
			    text: 'See All Sample Data',
			    title: 'Sample Data',
			    href: "/" + rName + "/showOneMipAllSampsData/" + mipName,
			    id: "samp"
			}).appendTo('#linkPanelBody');
		
			var cols = 10;
			var linkPre = "/" + rName + "/showOneMipOneSampData/" + mipName + "/" ;
			var mouseOverC = "#999";
			var mouseLeaveC = "#FFF";
			var addTo = "#sampLinks";
			createLinksTable(addTo, linkPre, sampNames,cols, mouseOverC, mouseLeaveC);
		}).catch(logRequestError).then(function(){
			//done loading names for mip fam
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});