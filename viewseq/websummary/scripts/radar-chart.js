/* Practically all this code comes from https://github.com/alangrafu/
 * radar-chart-d3 I only made some additions and aesthetic adjustments to make 
 * the chart look better (of course, that is only my point of view), such as a 
 * better placement of the titles at each line end, adding numbers that reflect
 *  what each circular level stands for, not placing the last level and slight
 *  differences in color. For a bit of extra information check the blog about 
 *  it: http://nbremer.blogspot.nl/2013/09/making-d3-radar-chart-look-bit-
 *  better.html
 * 
 * This is an update of the above update to radar-chart.js library. I've 
 * improved the way it constructs the plot by drawing the axes vertically, 
 * adding labels, tickmarks and arbirary scale values, and then animating 
 * rotating them into position and then make the ticks, scale values and labels 
 * I've also changed it to one series per plot and multiple plots (at the 
 * moment). I've also added improved tooltips.
 */

	
function randomColors(total)
{
    var i = 360/total; // distribute the colors evenly on the hue range
    var r = []; // hold the generated colors
    for (var x=0; x<total; x++)
    {
        r.push(d3.hsl(i * x, 0.8, 0.4)); // you can also alternate the saturation and value for even more contrast between the colors
    }
    return r;
}

var RadarChart = {
draw: function(id, d, options){
var cfg = {
	 radius: 5,
	 w: 600,
	 h: 600,
	 factor: 1,
	 factorLegend: .85,
	 levels: 3,
	 maxValue: 0,
	 opacityArea: 0.5,
	 ToRight: 5,
	 TranslateX: 80,
	 TranslateY: 30,
	 ExtraWidthX: 100,
	 ExtraWidthY: 100,
	 color: d3.scale.category10(),
	 labelFontSize: 14,
	 axisFontSize: 10,
	 levelStrokeWidth: 0.3,
	 axesStrokeWidth: 1,
	 showLastLevel:false,
	 axisIsPercent:false,
	 tAnim:3000,
	 title:"",
	 titlePad:50,
	 pointSize:5
	};
	
	if('undefined' !== typeof options){
	  for(var i in options){
		if('undefined' !== typeof options[i]){
		  cfg[i] = options[i];
		}
	  }
	}
	
	// get data
	var allAxis = (d[0].map(function(i, j){return i.axis}));
	var total = allAxis.length;
	var radius = Math.min(cfg.w/2, cfg.h/2);
	var radians = 2*Math.PI;
	
	// remove existing svg from container and add a new one
	d3.select(id).select("svg").remove();
	var svg = d3.select(id)
			.append("svg")
			.attr("width", cfg.w+cfg.ExtraWidthX)
			.attr("height", cfg.h+cfg.ExtraWidthY)
			.attr("transform", "translate(" + cfg.TranslateX + "," + cfg.TranslateY + ")");
	
	// plot title
	if (cfg.title!="") {
		var title = svg.append("text")
		.attr("class", "title")
		.text(cfg.title)
		.style("font-family", "sans-serif")
		.style("font-size", cfg.labelFontSize+"px")
		.attr("text-anchor", "middle")
		.attr("x", cfg.w/2)
		.attr("y", 0)
		.attr("dy", "1em");
	}
	
	var g = svg.append("g")
	.attr("id", "outerg")
	.attr("transform", "translate(0,"+cfg.titlePad+")");
	
	// setup animations
	var t0 = g.transition().duration(cfg.tAnim/3);
	
	// axis lines	
	var axis = g.selectAll(".axis")
	.data(allAxis)
	.enter()
	.append("g")
	.attr("class", "axis");

	axis.append("line")
		.attr("x1", cfg.w/2)
		.attr("y1", cfg.h/2)
		.attr("x2", cfg.w/2)
		.attr("y2", 0)
		.attr("class", "line")
		.style("stroke", "grey")
		.style("stroke-width", cfg.axesStrokeWidth+"px");
	
	// axis labels
	axis.append("text")
	.attr("class", "axisLabel")
	.attr("id", function(d, i) {return "axisLabel_"+i;})
	.text(function(d){return d})
	.style("font-family", "sans-serif")
	.style("font-size", cfg.labelFontSize+"px")
	.attr("text-anchor", function(d,i) {
		if ((i*(360/total))>180) {
			return "start";
		} else {
			return "end";
		}
	})
	.attr("x", cfg.w/2)
	.attr("y", 0)
	.attr("dy", "-0.5em")
	.attr("transform", function(d,i) {
		if ((i*(360/total))>180) {
			return "rotate(90,"+cfg.w/2+",0)"; 
		} else {
			return "rotate(-90,"+cfg.w/2+",0)";
		}
	})
	.style("opacity", "0");
	
	// choose wether to draw the outermost ticks
	var subi=1;
	if (cfg.showLastLevel) {
		subi=0;
	}

	// axis ticks
	var Format = d3.format("s")
	for(var j=subi; j<cfg.levels; j++){
		axis.append("text")
		.attr("class", "axisTick")
		.attr("class", "axisTick"+j)
		.style("font-family", "sans-serif")
		.style("font-size", cfg.axisFontSize+"px")
		.attr("text-anchor", "middle")
		.attr("x", cfg.w/2)
		.attr("y", j*((cfg.h/2)/cfg.levels))
		.attr("dx", function(d, i) {
			if ((i*(360/total))>180) {
				return "-1.0em"; 
			} else {
				return "1.0em";
			}
		})
		.style("opacity", "0")
		.text(function(d, i){
			return Format(((cfg.levels-j)*(cfg.maxValue[d]/cfg.levels)).toPrecision(2));
		});
	}
	
	// animate them positioning themselves
	t0.selectAll(".axis")
    .attrTween("transform", function(d, i) {
    	var start = "rotate(0,"+cfg.w/2+","+cfg.h/2+")";
    	var end = "rotate("+i*(360/total)+","+cfg.w/2+","+cfg.h/2+")";
    	var rtrans = d3.interpolateString(start,end);
       	return rtrans;
    });
	t1 = t0.transition().duration(cfg.tAnim/3);
	t1.selectAll(".axisLabel").style("opacity", "1");
	tloop = t1.transition().duration(cfg.tAnim/20);
	for(var j=(cfg.levels-subi); j>-1; j--){
		tloop.selectAll(".axisTick"+j).style("opacity", "1");
		tloop = tloop.transition().duration(cfg.tAnim/20);
	}

	// draw the connecting lines for each specified level	
	//Circular segments
	for(var j=0; j<cfg.levels-subi; j++){
	  var levelFactor = radius*((j+1)/cfg.levels);
	  g.selectAll(".levels")
		.data(allAxis)
		.enter()
		.append("svg:line")
		.attr("x1", function(d, i){return levelFactor*(1-Math.sin(i*radians/total));})
		.attr("y1", function(d, i){return levelFactor*(1-Math.cos(i*radians/total));})
		.attr("x2", function(d, i){return levelFactor*(1-Math.sin((i+1)*radians/total));})
		.attr("y2", function(d, i){return levelFactor*(1-Math.cos((i+1)*radians/total));})
		.attr("class", "webline")
		.style("stroke", "grey")
		.style("stroke-opacity", "0")
		.style("stroke-width", cfg.levelStrokeWidth+"px")
		.attr("transform", "translate("+(cfg.w/2-levelFactor)+","+(cfg.h/2-levelFactor)+")");
	}

	// animate them positioning themselves
	//var t1 = t0.transition().duration(1000);
	t0.selectAll(".webline").style("stroke-opacity","0.75");
	
	g.selectAll(".area")
	 .data(d)
	 .enter()
	 .append("polygon")
	 .attr("class", "radar-chart-area")
	 .style("stroke-width", "2px")
	 .style("stroke-opacity", "0")
	 .style("stroke", "grey")
	 .attr("points",function(d) {
		 var str="";
		 for(var pti=0;pti<total;pti++){
			 init_x = cfg.w/2;
			 len_on_axis = ((d[pti]["value"]/cfg.maxValue[d[pti]["axis"]])*cfg.h/2);
			 new_x = init_x+(Math.sin(pti*((2*Math.PI)/total))*len_on_axis);
			 new_y = (cfg.h/2)-Math.cos(pti*((2*Math.PI)/total))*len_on_axis;
			 str=str+new_x+","+new_y+" ";
		 }
		 return str;
	  })
	 .style("fill", cfg.color)
	 .style("fill-opacity", 0)
	 .on('mouseover', function (d){
						z = "polygon."+d3.select(this).attr("class");
						g.selectAll("polygon")
						 .transition(200)
						 .style("fill-opacity", 0.1); 
						g.selectAll(z)
						 .transition(200)
						 .style("fill-opacity", .7);
					  })
	 .on('mouseout', function(){
						g.selectAll("polygon")
						 .transition(200)
						 .style("fill-opacity", cfg.opacityArea);
	 })
	 ;

	// value points
	for(var j=0; j<total; j++){
		g.append("circle")
		 .data(d)
		 .attr("class", "radar-chart-points")
		 .attr("title", function(d, i) {
			 return d[j]["value"];
		 })
		 .attr("cx", function(d,i) {
			 init_x = cfg.w/2;
			 len_on_axis = ((d[j]["value"]/cfg.maxValue[d[j]["axis"]])*cfg.h/2);
			 new_x = init_x+(Math.sin(j*((2*Math.PI)/total))*len_on_axis);
			 return new_x;
		 })
		 .attr("cy", function(d,i) {
			 len_on_axis = ((d[j]["value"]/cfg.maxValue[d[j]["axis"]])*cfg.h/2);
			 new_y = (cfg.h/2)-Math.cos(j*((2*Math.PI)/total))*len_on_axis;
			 return new_y;
		 })
		 .attr("r",cfg.pointSize)
		 .style("fill", cfg.color)
		 .style("fill-opacity", 0)
		 .style("stroke-width", "2px")
		 .style("stroke-opacity", "0")
		 .style("stroke", "black")
		 .on('mouseover', function (d, j){
			 newX =  parseFloat(d3.select(this).attr('cx'))+15;
			 newY =  parseFloat(d3.select(this).attr('cy'))+10;
			 d3.select(this).transition(200).attr('r', cfg.pointSize*2)

			var tooltip_box = g.append("rect")
						.attr("class", id.replace("#","")+"_tt")
						.attr("x",0)
						.attr("y",0)
						.attr("width",0)
						.attr("height",0)
						.attr("fill","Ivory")
						.attr("fill-opacity",0)
						.style("stroke", "grey")
						.style("stroke-opacity", "0")
						.style("stroke-width", "2px");

			 var tooltip = g.append('text')
						.attr("class", id.replace("#","")+"_tt")
					   	.style('opacity', 0)
					   	.style('font-family', 'sans-serif')
					   	.style('font-size', cfg.axisFontSize+'px')
					 	.attr('x', newX)
						.attr('y', newY)
						.text(d3.select(this).attr('title'))
						.transition(200)
						.style('opacity', 1);
			 			 
			 tooltip.each(function(){
				 var thisbox = this.getBBox();
				 tooltip_box.attr('x', thisbox.x-5)
					.attr('y', thisbox.y-3)
					.attr("width", thisbox.width+10)
					.attr("height", thisbox.height+6)
					.transition(200)
					.style('stroke-opacity', 1)
					.style('fill-opacity', 1);
			 });
		  })
		  .on('mouseout', function(){
			  d3.select(this).transition(200).attr('r', cfg.pointSize);
			  d3.selectAll("."+id.replace("#","")+"_tt")
			  	.transition(200)
				.remove();
		  });
	}
	
	// add final animation
	t2 = tloop.transition().duration(cfg.tAnim/3);
	t2.selectAll(".radar-chart-area")
		.style("fill-opacity",cfg.opacityArea)
		.style("stroke-opacity", "1");
	t2.selectAll(".radar-chart-points")
		.style("fill-opacity","1")
		.style("stroke-opacity", "1");
}
};
