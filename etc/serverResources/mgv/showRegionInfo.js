/**
 * 
 */


$(document).ready(function() {
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = getRootName();
	var mipRregionName = locSplit.pop();
	setHeadTitle(mipRregionName);
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		getJSON("/" + rName + "/getMipsForRegion/" + mipRregionName).then(function (mipNames) {
			names["regionName"] = mipRregionName;
			addDiv("body", "topNav");
			createNavBar("#topNav", names);
			addMainDiv("body", true);
			setHeadTitle(mipRregionName);
			$("#jumboTitle").html(mipRregionName);
			addPanelWithDiv(".container.theme-showcase","mipsLinks", "Mip");
			//addPanelWithDiv(".container.theme-showcase","sampLinks", "Ge");
			var cols = 10;
			var linkPre = "/" + rName + "/showMipSeqs/";
			var mouseOverC = "#999";
			var mouseLeaveC = "#FFF";
			var addTo = "#mipsLinks";
			createLinksTable(addTo, linkPre, mipNames["mips"], cols, mouseOverC, mouseLeaveC);		
			addH1 ("#mainContent", "General Infoformation Per Target");
			addDiv("#mainContent", "tarInfoTab");
			postJSON("/" + rName + "/getMipTarInfoForGenome", {"genome":names["primaryGenome"], "mipTars": mipNames["mips"]}).then(function (tarInfo) {
				var tarsTable =  new njhTable("#tarInfoTab", tarInfo, names["primaryGenome"] + "_" + mipRregionName + "_tarInfos", false);	
			}).catch(logRequestError).then(function(){
				//done loading names
			});
		}).catch(logRequestError).then(function(){
			//done loading names
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});


