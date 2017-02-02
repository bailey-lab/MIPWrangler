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
		var gifLoading = prsentDivGifLoading();
		postJSON('/' + rName + '/getMipSeqs/' + mipName, {mipTar:mipName}).then(function (mainData) {
			var sesUid = mainData["sessionUID"];
			var SeqViewer = new njhSeqView("#dnaViewer", mainData);
			gifLoading.remove();
			setUpCloseSession(sesUid);
		}).catch(logRequestError).then(function(){
			//done loading names
		});
	}).catch(logRequestError).then(function(){
		//done loading 
		removeAllDivGifLoading();
	});
});
