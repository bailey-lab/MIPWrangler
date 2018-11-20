$(document).ready(function() {
	
	var rName = getRootName();
	var gifLoading = prsentDivGifLoading();
	getJSON("/" + rName + "/getAllNames").then(function (names) {
		postJSON("/" + rName + "/getGenomeLens", {"genomes":[names["primaryGenome"]]}).then(function (genomeLens) {
			addDiv("body", "topNav");
			createNavBar("#topNav", names);
			addMainDiv("body", true);
			setHeadTitle("Mips!");
			$("#jumboTitle").html("Mips on Genomes");
			addPanelWithDiv(".container.theme-showcase","mipsLinks", "Regions");
			var cols = 10;
			var linkPre = "/" + rName + "/showRegionInfo/";
			var mouseOverC = "#999";
			var mouseLeaveC = "#FFF";
			var addTo = "#mipsLinks";
			createLinksTable(addTo, linkPre, names["mipRegions"], cols, mouseOverC, mouseLeaveC);		
			addH1 ("#mainContent", "Regions General Infoformation");
			addDiv("#mainContent", "regionInfoTab");
			addPanelOnlyHead("#mainContent", "General Infoformation Per Target For " + names["primaryGenome"]);
			addDiv("#mainContent", "tarInfoTab");
			addPanelOnlyHead("#mainContent", "General Infoformation Per Target For All Genomes");
			addDiv("#mainContent", "tarInfoTabAll");
			addPanelOnlyHead("#mainContent", "Input Mip Info Table");
			addDiv("#mainContent", "inputTarInfo");
			
			postJSON("/" + rName + "/getMipRegionInfoForGenome", {"genome":names["primaryGenome"]}).then(function (regionInfo) {
				var regionTable =  new njhTable("#regionInfoTab", regionInfo, names["primaryGenome"] + "_regionInfo", false);	
				var tarUrls = ["/" + rName + "/getMipTarInfoForGenome" ]
				tarUrls.push("/" + rName + "/getAllMipTarInfoAllGenomes" );
				tarUrls.push("/" + rName + "/getInfoMipArmsInfo" );
				Promise.all(tarUrls.map(function(tarUrl){
					return postJSON(tarUrl, {"genome": names["primaryGenome"], "mipTars": names["mips"]});
				})).then(function(allTarInfos){
					
					var tarInfo = allTarInfos[0]
					var allGenomeTarInfo = allTarInfos[1]
					var inputTarInfoTable = allTarInfos[2]
					var tarsTable =  new njhTable("#tarInfoTab", tarInfo, names["primaryGenome"] + "_tarInfos", false);	
					var tarInfoTabAllTable =  new njhTable("#tarInfoTabAll", allGenomeTarInfo, "allGenomes" + "_tarInfos", false);	
					var inputTarInfoTable =  new njhTable("#inputTarInfo", inputTarInfoTable, "mip_arms.tab.txt", false);	
					
					addH1 ("#mainContent", names["primaryGenome"]);
					addDiv("#mainContent", "genomeDiv");
					d34.select("#genomeDiv").style("margin", "10px");
					
					var margin = {top: 50, right: 20, bottom: 30, left: 100},
					    width = window.innerWidth *.90 - margin.left - margin.right,
					    height = window.innerHeight *.80 - margin.top - margin.bottom;
					
					var yScale = d34.scaleBand()
					    .range([0, height])
					    .domain(genomeLens[names["primaryGenome"]].map(function(d){ return d.name;}) )
					    .padding(.2)
					    .round(true);

					var xScale = d34.scaleLinear()
					    .rangeRound([0, width])
					    .domain([0, d34.max(genomeLens[names["primaryGenome"]], function(d) { 
						 return d.len;   })]);

					var	xAxis = d34.axisBottom(xScale).tickFormat(function(e){
		        if(Math.floor(e) != e){
		            return;
		        }
		        return e;
						});
		     var yAxis = d34.axisLeft(yScale);

					var svg = d34.select("#genomeDiv").append("svg")
					    .attr("width", width + margin.left + margin.right)
					    .attr("height", height + margin.top + margin.bottom)
					    .attr("id", "njhGenome")
					  .append("g")
					    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
					svg.append('g').classed('x-axis-group axis', true);
					svg.append('g').classed('y-axis-group axis', true);
					svg.select('.x-axis-group.axis')
			      .attr('transform', 'translate(0,' + height + ')')
			      .call(xAxis);

				  svg.select('.y-axis-group.axis')
				      .call(yAxis);
				  
				  var tooltip = d3.select("body")
						.append("div")
							.style("position", "absolute")
							.style("visibility", "hidden")
							.style("background-color", "#88aaaa")
							.style("width", "300")
							.attr("id", "popHover");
				  
				  var bars = svg.selectAll("rect")
		      		.data(genomeLens[names["primaryGenome"]])
		      		.enter()
		      			.append("g").attr("class", "subbar")
								  .append("rect")
						      .attr("height", yScale.bandwidth())
						      .attr("y", function(d) { return yScale(d.name)})
						      .attr("x", function(d) { return xScale(0); })
						      .attr("width", function(d) { return xScale(d.len); })
						      .attr("id", function(d) { return d.name; })
						      .style("fill", function(d) { return "#005AC8"; });
					
				}).catch(logRequestError).then(function(){
					//done loading names
				});
				
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

