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
			
			addDiv("body", "genome");
			d34.select("#genome").style("margin", "10px");
			
			var margin = {top: 50, right: 20, bottom: 30, left: 100},
			    width = window.innerWidth *.90 - margin.left - margin.right,
			    height = window.innerHeight *.80 - margin.top - margin.bottom;
			
			 console.log(genomeLens[names["primaryGenome"]].map(function(d){ return d.name;}))
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

			var svg = d34.select("#genome").append("svg")
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
		//done loading 
		removeAllDivGifLoading();
	});

	

});

