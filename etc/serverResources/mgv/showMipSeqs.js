$(document).ready(function(){
	var locSplit = window.location.toString().split(/[\/]+/);
	var rName = getRootName();
	var mipName = locSplit.pop();
	setHeadTitle(mipName);
	addDiv("body", "viewer");
	var gifLoading = prsentDivGifLoading();
	postJSON('/' + rName + '/getMipSeqs/' + mipName, {mipTar:mipName}).then(function (mainData) {
		var sesUid = mainData["sessionUID"];
		var SeqViewer = new njhSeqView("#viewer", mainData);
		gifLoading.remove();
		setUpCloseSession(sesUid);
	}).catch(function(err){
  		removeAllDivGifLoading();
  		logRequestError(err);
  	});
});
