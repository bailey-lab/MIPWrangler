$(document).ready(function(){
	//get current name from window location
	var rName = getRootName();
	$("title", "head").html(rName);
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		addDiv("body", "topNav");
		createNavBar("#topNav", names);
		addMainDiv("body", true);
		$("#jumboTitle").html("Ruh Roh");
		$(".jumbotron", "#mainContent").append("<h3>You attempted to go to a page that doesn't exist, or it could be there is no data for where you attempted to go.</h3>");
	}).catch(logRequestError).then(function(){
		//done loading 
	});	
});


