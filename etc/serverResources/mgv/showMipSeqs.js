$(document).ready(function(){
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = getRootName();
	var mipName = locSplit.pop();
	setHeadTitle(mipName);
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		$("#jumboTitle").html(mipName);
		addDiv("#mainContent", "dnaViewer");
		addDiv("#mainContent", "locTable");
		var gifLoading = prsentDivGifLoading();
		postJSON('/' + rName + '/getMipSeqs/' + mipName, {mipTar:mipName}).then(function (mainData) {
			getJSON('/' + rName + '/getMipGenomeLocs/' + mipName).then(function (locData) {
				var sesUid = mainData["sessionUID"];
				var SeqViewer = new njhSeqView("#dnaViewer", mainData);
				var locTable =  new njhTable("#locTable", locData, mipName + "_genLocs", false);	
				gifLoading.remove();
				setUpCloseSession(sesUid);
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
