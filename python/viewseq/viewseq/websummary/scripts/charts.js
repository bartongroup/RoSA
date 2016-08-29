var charts = (function () {

var module = {};

var plot_col = 'steelblue';
var selection_col = 'red';
var highlight_col = 'lightgreen';
var point_radius = 3; //radius of plotted points

module.selected_point = null;
//module.selected_row = null;

var WRAPPER = "wrapper"
var CANVAS = "canvas"
var SELECT_CANVAS = "select_canvas"
var HIGHLIGHT_CANVAS = "highlight_canvas"

var tooltip_enabled = false;
var highlight_enabled = false;


var chart_data_list = {}; //dictionary of chart objects, to be indexed by id

function init_chart_data(id)
{
    // svg : reference to svg, initially null
    // xscale : reference to xscale, initially null
    // yscale : reference to yscale, initially null
    // xaxis : reference to x axis, initially null
    // yaxis : reference to y axis, initially null
    // highlight_points : points currently highlighted on this chart
    var chart_data =
    {
        svg : null,
        quadTree : null,
        xscale : null,
        yscale : null,
        xaxis : null,
        yaxis : null,
        highlight_points : [],
        highlight_context : null,
        select_context : null
    }
    chart_data_list[id] = chart_data;
}

// Initialise document for a scatterplot
// id: id of the element in which the scatterplot should be placed
function init_scatter_plot(id)
{
    // set up holder for chart data
    init_chart_data(id);

    // wrapper div
    init_wrapper(id);

    //canvas
    var body_canvas = document.createElement(CANVAS);
    body_canvas.id = CANVAS + id;

    //selection canvas
    var select_canvas = document.createElement(CANVAS);
    select_canvas.id = SELECT_CANVAS + id;
    select_canvas.className = SELECT_CANVAS;

    document.getElementById(WRAPPER + id).appendChild(body_canvas);

    //check if tooltip is enabled
    if (typeof(tooltip_callback) == "function")
    {
        tooltip_enabled = true;
    }

    //check if highlight is enabled
    if (typeof(highlight_callback) == "function")
    {
        highlight_enabled = true;

        //selection canvas
        var highlight_canvas = document.createElement(CANVAS);
        highlight_canvas.id = HIGHLIGHT_CANVAS + id;
        highlight_canvas.className = HIGHLIGHT_CANVAS;

        document.getElementById(WRAPPER + id).appendChild(highlight_canvas);
    }

    document.getElementById(WRAPPER + id).appendChild(select_canvas); //must be appended last or we lose mouse events?!
};

function init_wrapper(id)
{
    // wrapper div
    var wrapper_div = document.createElement('div');
    wrapper_div.className = WRAPPER;

    var width = $('#'+id).width(),
    aspect = 1;
    var height = width * aspect;

    wrapper_div.id = WRAPPER + id;
    wrapper_div.style.width = width + 'px';
    wrapper_div.style.height = height + 'px';

    document.getElementById(id).appendChild(wrapper_div);
}


// Draw a scatterplot
// scatterplot_data: the json formatted data for the scatterplot
// plot: configuration of the plot
module.draw_scatter_plot = function (plot) {

    if (chart_data_list[plot.id] == null)
    {
        init_scatter_plot(plot.id);
    }
    var c = chart_data_list[plot.id];

    //remove any old plot, in case we are redrawing
    if (c.svg != null)
    {
        c.svg.selectAll("*").selectAll("*").remove();
    }

    c.quadTree = d3.quadtree(plot.data,
                             function(d) {return +d[plot.xindex]; },
                             function(d) {return +d[plot.yindex]; })

    //set up width and height
    var width = $('#'+plot.id).width(),
    aspect = 1;
    var height = width * aspect;

    // Set up plot margins
    var margin = {top: 10, right: 10, bottom: 40, left: 70},
    plot_width = width - margin.left - margin.right,
    plot_height = height - margin.top - margin.bottom;

    // add the tooltip area to the webpage
    var tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0); //needs to be 0, as this plonks a box on the web page ready to use as tooltip

    /* define x and y scales: these translate data into x.y positions on the graph
     as they're often not the same. */

    var scatterplot_data = plot.data;

    // Define x scale
    if (plot.scale_type == "log")
    {
        // filter 0 values out of the data
        scatterplot_data = scatterplot_data.filter(function (d)
                            {
                                if(d[plot.xindex] == 0)
                                {
                                    return false;
                                }
                                else if (d[plot.yindex] == 0)
                                {
                                    return false;
                                }
                                else
                                {
                                    return true;
                                }
                            });

        var maxx = d3.max(scatterplot_data, function(d) { return d[plot.xindex]; })
        var maxy = d3.max(scatterplot_data, function(d) { return d[plot.yindex]; })
        var minx = d3.min(scatterplot_data, function(d) { return d[plot.xindex]; })
        var miny = d3.min(scatterplot_data, function(d) { return d[plot.yindex]; })

        c.xscale = d3.scaleLog()
        .domain([minx, maxx])
        .rangeRound([0, plot_width]).nice();

        c.yscale = d3.scaleLog()
        .domain([miny, maxy])
        .rangeRound([plot_height, 0]).nice();
    }
    else
    {
        var maxx = d3.max(scatterplot_data, function(d) { return d[plot.xindex]; })
        var maxy = d3.max(scatterplot_data, function(d) { return d[plot.yindex]; })

        //default to linear
        c.xscale = d3.scaleLinear()
        .domain(d3.extent(scatterplot_data, function(d) { return d[plot.xindex]; }))
        .rangeRound([0, plot_width]);

        c.yscale = d3.scaleLinear()
        .domain(d3.extent(scatterplot_data, function(d) { return d[plot.yindex]; }))
        .rangeRound([plot_height, 0]);
    }

    /* define the scatterplot: we want to put some svg in the scatterplot div,
     with width and height as defined earlier */
    c.svg = d3.select('#' + WRAPPER + plot.id)
    .append("svg")
    .attr("class", "absolute_svg")
    .attr("id", "svg" + plot.id)
    .attr("width", width)
    .attr("height", height)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // plot x and y axes
    draw_axes(plot, plot_height, plot_width, c, margin.bottom, margin.left);

    // draw the points on the canvas
    // set up the canvas
    var canvas = d3.select("#" + CANVAS + plot.id)
    .attr("width", plot_width + point_radius)
    .attr("height", plot_height + point_radius)
    .style("transform", "translate(" + margin.left + "px ," + margin.top + "px )");

    //make a selection canvas the same size as the plot canvas
    var selectcanvas = d3.select("#" + SELECT_CANVAS + plot.id)
    .attr("width", plot_width + point_radius)
    .attr("height", plot_height + point_radius)
    .style("transform", "translate(" + margin.left + "px ," + margin.top + "px )");

    if (highlight_enabled)
    {
        //make a highlight canvas too
        var highlightcanvas = d3.select("#" + HIGHLIGHT_CANVAS + plot.id)
        .attr("width", plot_width + point_radius)
        .attr("height", plot_height + point_radius)
        .style("transform", "translate(" + margin.left + "px ," + margin.top + "px )");

        var highlight_context = highlightcanvas.node().getContext('2d');
        c.highlight_context = highlight_context;
    }

    var context = canvas.node().getContext('2d');
    var select_context = selectcanvas.node().getContext('2d');
    c.select_context = select_context;

    c.txscale = c.xscale;
    c.tyscale = c.yscale;
    draw(scatterplot_data, context, plot, c);

    var zoom = d3.zoom()
               .scaleExtent([1, 30])
               .on('zoom', function(event) { onZoom(event, scatterplot_data, context, plot, c); });

    // zoom handler
    selectcanvas.call(zoom)

    // on onclick handler
    selectcanvas.on("click", function(event) { onClick(this, plot, c); } );

    // mouseover handler
    selectcanvas.on("mousemove", function(event) { onMouseOver(this, plot, tooltip, c); });

    // mouseout handler
    selectcanvas.on("mouseout", function(event) { onMouseOut(this, plot, tooltip); });

    //selectcanvas.call(zoomBehaviour);

    module.update_scatter_selection(plot);
};

function onZoom(event, data, context, plot, chart_data) {

    transform = d3.event.transform;

    // transform the scale
    chart_data.txscale = transform.rescaleX(chart_data.xscale);
    chart_data.tyscale = transform.rescaleY(chart_data.yscale);

    draw(data, context, plot, chart_data, transform);

    // transform and translate the axes
    chart_data.svg.select(".x.axis").call(chart_data.xaxis.scale(transform.rescaleX(chart_data.xscale)));
    chart_data.svg.select(".y.axis").call(chart_data.yaxis.scale(transform.rescaleY(chart_data.yscale)));

    //event.stopPropagation(); //try to avoid mad scrolling when we get to the full zoom extent
}

// Event handler for clicking on a scatterplot point corresponding to a gene
function onClick(canvas, plot, chart_data) {

    if (d3.event.defaultPrevented) return; // ignore drag

    if (typeof(update_plots_callback) == "function") //check update_plots has been defined
    {
        var closest = getMouseOverPoint(canvas, plot, chart_data);

        if (closest)
        {
            //something new was selected
            module.selected_point = closest[plot.iindex];
            update_plots_callback(true);
        }
        else if (module.selected_point)
        {
            //something was selected but now it's not
            module.selected_point = null;
            update_plots_callback(true);
        }
    }
}

function onMouseOver(canvas, plot, tooltip, chart_data)
{
    var d = getMouseOverPoint(canvas, plot, chart_data);

    if (tooltip_enabled)
    {
        if (d)
        {
            //show tooltip
            var view_width = Math.max(document.documentElement.clientWidth, window.innerWidth || 0);
            var leftpos = d3.event.pageX + 15;
            if (d3.event.pageX + 15 + 280 > view_width) //tooltip will hang off page
            {
                leftpos = d3.event.pageX - 290; //switch tooltip to display on left instead of right
            }
            tooltip.transition()
            .duration(200)
            .style("opacity", 100);

            tooltip.html(tooltip_callback(d,plot))
            .style("left", leftpos + "px")
            .style("top", (d3.event.pageY - 28) + "px");
        }
        else
        {
            tooltip.transition()
            .duration(300)
            .style("opacity", 0);
        }
    }

    if (highlight_enabled)
    {
        if (d)
        {
            //highlight related points
            points = highlight_callback(d, plot); //returns ids of points to highlight

            var related = [];
            for (point of points)
            {
                // don't highlight current point
                if (point != d)
                {
                    var next = plot.data.filter(function(v) { return v[plot.iindex] == point[plot.iindex]; });
                    if ((next[0][plot.xindex] != null) && (next[0][plot.yindex] != null))
                    {
                        related.push(next[0]);
                    }
                }
            }
            chart_data.highlight_points = related;
            drawHighlightedPoints(plot, chart_data);
        }
        else
        {
            //turn off old selection
            chart_data.highlight_points = [];
            chart_data.highlight_context.clearRect(0, 0, chart_data.highlight_context.canvas.width,
                                                   chart_data.highlight_context.canvas.height);
        }
    }
}

function onMouseOut(canvas, plot, tooltip)
{
    tooltip.transition()
        .duration(300)
        .style("opacity", 0);
}

function getMouseOverPoint(canvas, plot, chart_data)
{
    var mouse = d3.mouse(canvas);

    // map the point to the data space
    var xClicked = chart_data.txscale.invert(mouse[0]);
    var yClicked = chart_data.tyscale.invert(mouse[1]);

    // find the closest point in the dataset to the clicked point
    var closest = chart_data.quadTree.find(xClicked, yClicked);

    // yes, find has a radius option, which is all very well
    // but it seems to be in chart units, not pixels, so it varies depending on the scale
    // and it's not clear how to scale a radius by different x and y scales
    // and frankly this works a lot better
    var dX = chart_data.txscale(closest[plot.xindex]);
    var dY = chart_data.tyscale(closest[plot.yindex]);

    // register the click if the clicked point is in the radius of the point
    var distance = euclideanDistance(mouse[0], mouse[1], dX, dY);

    if (distance < point_radius)
    {
        return closest;
    }
    else
    {
        return null;
    }
}


function euclideanDistance(x1, y1, x2, y2) {
    return Math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

function draw(data, context, plot, chart_data, transform)
{
    // set the point details
    context.clearRect(0, 0, context.canvas.width, context.canvas.height);
    context.fillStyle = plot_col;
    context.strokeWidth = 1;
    context.strokeStyle = 'black';
    context.globalAlpha = 0.5; //opacity

    // draw the points
    var i = -1, n = data.length, d, active;
    while (++i < n) {
        d = data[i];
        drawPoint(d, context, plot, chart_data);
        if (d[plot.iindex] == module.selected_point)
        {
            active = d;
        }
    }

    // sort out any points which were selected before we redrew
    if (active)
    {
        drawSelectedPoint(active, plot, chart_data);
    }

    //clear any highlighted points
    if (highlight_enabled)
    {
        drawHighlightedPoints(plot, chart_data);
    }
}

function drawSelectedPoint(point, plot, chart_data)
{
    var context = chart_data.select_context;
    context.clearRect(0, 0, context.canvas.width, context.canvas.height);

    context.fillStyle = selection_col;
    context.globalAlpha = 1;
    drawPoint(point, context, plot, chart_data);
}

function drawHighlightedPoints(plot, chart_data)
{
    var context = chart_data.highlight_context;
    context.clearRect(0, 0, context.canvas.width, context.canvas.height);

    context.fillStyle = highlight_col;
    context.globalAlpha = 1;
    for (point of chart_data.highlight_points)
    {
        drawPoint(point, context, plot, chart_data);
    }
}

function drawPoint(point, context, plot, chart_data)
{
    var cx, cy;

    // transform and translate points (for zoom and pan)
    cx = chart_data.txscale(point[plot.xindex]);
    cy = chart_data.tyscale(point[plot.yindex]);

    context.beginPath();
    context.arc(cx, cy, point_radius, 0, 2 * Math.PI);
    context.closePath();
    context.fill();
    //context.stroke(); too slow on redraw
}

// Initialise document for a barchart
// id: id of the element in which the barchart should be placed
init_barchart = function (id)
{
    //check if already initialised, if so return
    if (document.getElementById(WRAPPER + id) != null)
    {
        return;
    }
    init_chart_data(id);
    init_wrapper(id);
};

module.draw_barchart = function (plot)
{
    init_barchart(plot.id);
    var c = chart_data_list[plot.id];

    //set up width and height
    //var width = $('#'+plot.id).width(),
    var width = $('#tufu-scatter').width(), //hack because plots on other tabs don't have width/height TODO fix!
    aspect = 1;
    var height = width * aspect;

    // Set up plot margins
    // higher top margin to allow text values to be written above bars
    var margin = {top: 30, right: 40, bottom: 40, left: 70};
    var plot_width = width - margin.left - margin.right;
    var plot_height = height - margin.top - margin.bottom;

    // scale to ordinal because x axis is not numerical
    //c.xscale = d3.scaleOrdinal().rangeRoundBands([0, plot_width], 0.1);
    c.xscale = d3.scaleBand().range([0, plot_width]);

    //scale to numerical value by height
    c.yscale = d3.scaleLinear().range([plot_height, 0]);

    c.svg = d3.select('#' + WRAPPER + plot.id)
                  .append("svg")  //append svg element inside #chart
                  .attr("class", "absolute_svg")
                  .attr("width", width) //set width
                  .attr("height", height) //set height
                  .append("g")
                  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    c.xscale.domain(plot.data.map(function(d){ return d[plot.xname]; }));
    c.yscale.domain([0, d3.max(plot.data, function(d){return d[plot.yname]})]);

    var bar = c.svg.selectAll(".bar")
        .data(plot.data)
        .enter().append("g")
        .attr("class", "bar");

    bar.append("rect")
          .attr("class", "bar")
          .attr("x", function(d) { return c.xscale(d[plot.xname]); })
          .attr("width", c.xscale.step())
          .attr("y", function(d) { return c.yscale(d[plot.yname]); })
          .attr("height", function(d) { return plot_height - c.yscale(d[plot.yname]); });

    bar.append("text")
      .attr("x", function(d) { return c.xscale(d[plot.xname]) + c.xscale.step()/2; }) //c.xscale.rangeBand()/2; })
      .attr("y", function(d) { return c.yscale(d[plot.yname]) - 15; })
      .attr("dy", ".75em")
      .style("text-anchor", "middle")
      .text(function(d) { return d[plot.ydisplay]; });

    draw_axes(plot, plot_height, plot_width, c, margin.bottom, margin.left);
};


// Initialise document for a histogram
// id: id of the element in which the histogram should be placed
init_histogram = function (id)
{
    //check if already initialised, if so return
    if (document.getElementById(WRAPPER + id) != null)
    {
        return;
    }

    init_chart_data(id);
    init_wrapper(id);
};

// Draw a histogram
// hist_data: json formatted data for the histogram
// plot: configuration for the plot
module.draw_histogram = function (plot)
{
    init_histogram(plot.id);
    var c = chart_data_list[plot.id];

    //remove any old plot, in case we are redrawing
    if (c.svg != null)
    {
        c.svg.selectAll("*").selectAll("*").remove();
    }

    // set up width and height
    var width = $('#'+plot.id).width(),
    aspect = 1;
    var height = width * aspect;

    // Set up plot margins
    var margin = {top: 10, right: 10, bottom: 40, left: 70},
    plot_width = width - margin.left - margin.right,
    plot_height = height - margin.top - margin.bottom;

    var maxval = d3.max(plot.data, function(d) { return d[plot.xindex]; }) + 1; //correction for small values
    // Define x scale
    if (plot.scale_type == "log")
    {
        c.xscale = d3.scaleLog()
        .domain([1, maxval])
        .range([0, plot_width]);
    }
    else
    {
        //default to linear
        c.xscale = d3.scaleLinear()
        .domain([0, maxval])
        .rangeRound([0, plot_width]);
    }

    // Generate a histogram using num_bins uniformly-spaced bins.
    if  (maxval < 16)
    {
        plot.num_bins = maxval; // 0 bin as well as 1-maxval
    }

    var data = d3.histogram()
        .thresholds(c.xscale.ticks(plot.num_bins))
        .domain(c.xscale.domain())  //need this to ensure first bin extends lower limit to 0
        .value(function(d) { return +d[plot.xindex]; })
        (plot.data);

    // Define y scale now we have the bins set up
    c.yscale = d3.scaleLinear()
        .domain([0, d3.max(data, function(d) { return d.length; })])
        .rangeRound([plot_height, 0]);

    var num_ticks = 10;

    //ordinalscale places ticks in centre of bar instead of between bars
    var ordinalscale = d3.scaleBand()
    .domain(c.xscale.ticks(plot.num_bins)) //map from bin to range
    .range([0, plot_width]);

    // Set up an svg for the histogram plot
    c.svg = d3.select('#' + WRAPPER + plot.id)
        .append("svg")
        .attr("class", "absolute_svg")
        .attr("id", "svg" + plot.id)
        .attr("preserveAspectRatio", "xMidYMid")
        .attr("viewBox", "0 0 " + width + " " + height)
        .attr("width", width)
        .attr("height", height)
        .append("g")
        .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

    // Draw the bars
    var bar = c.svg.selectAll(".bar")
        .data(data)
        .enter().append("g")
        .attr("class", "bar")
        .attr("transform", function(d)
              {
                return "translate(" + ordinalscale(d.x0) + "," + c.yscale(d.length) + ")";
              });

    bar.append("rect")
        .attr("x", 1)
        .attr("width", ordinalscale.step()) //c.xscale(data[0].x1) - c.xscale(data[0].x0) - 1) //rangeBand())
        .attr("height", function(d) { return plot_height - c.yscale(d.length); });

    module.update_hist_selection(plot);

    c.xscale = ordinalscale;
    draw_axes(plot, plot_height, plot_width, c, margin.bottom, margin.left);
};

module.draw_radar_plot = function (chartdata)
{
    var classes = {};
    $.each(chartdata.conditions, function(i, val)
    {
        classes[val] = i
    });

    var cols = randomColors(Object.keys(classes).length); //randomColors fn contained in radar chart library
    for(var j=0; j<cols.length; j++)
    {
        classes[Object.keys(classes)[j]] = cols[j];
    };

    $.each(chartdata.data, function(i, val)
    {
        var classstr = chartdata.conditions[i];
        var idstr = String(i).replace(/\//g,"").replace(/\./g,"");
        $('<div id="quality_'+idstr+'" class="'+classstr+'" style="float:left"></div>').appendTo($("#qualitydata"));

        var b = [];
        $.each(val, function(j,val2)
        {
            b.push({axis:j,value:val2});
        });

        var config = {
            w:400,
            h:400,
            maxValue:chartdata.limits,
            labelFontSize:14,
            axisFontSize:12,
            levels:3,
            showLastLevel:true,
            title:i,
            color:classes[classstr]
        };

        RadarChart.draw("#quality_"+idstr, [b], config);
    });

};



module.update_scatter_selection = function (scatter)
{
    //turn off old selection
    var chart_data = chart_data_list[scatter.id];
    chart_data.select_context.clearRect(0, 0, chart_data.select_context.canvas.width, chart_data.select_context.canvas.height);

    if (!module.selected_point)
    {
        return;
    }

    var datum = scatter.data.filter(function(v) { return v[scatter.iindex] == module.selected_point; });

    drawSelectedPoint(datum[0], scatter, chart_data);
}

module.update_hist_selection = function (histogram)
{
    // Update call for individual histograms so we can update after a scale change
    hist_index = histogram.iindex; //column name of gene index in this histogram
    chart_data_list[histogram.id].svg.selectAll(".bar")
    .style("fill", function(f)
           {
           if (f.length > 0) //if f is empty this data point can't be in the f bin
           {
           var result = $.grep(f, function(e){ return e[hist_index] == module.selected_point; });
           var update = result.length > 0;
           return update ? selection_col : plot_col;
           }
           return plot_col
           });
}

function draw_axes(plot, plot_height, plot_width, chart_data, bottom_margin, left_margin)
{
    var numticks = 5;

    max_xval = d3.max(plot.data, function(d) { return d[plot.xindex]; });

    if (plot.chart_type == "histogram")
    {
        //TODO ? use ordinal scale to get tickmarks in centre of bar instead of at start
        /*var xscale = d3.scaleBand()
        .domain(chart_data.xscale.ticks(plot.num_bins))
        .range([0, plot_width]);*/

        /* Define axes: use corresponding scales and attach to bottom and left of graph */
        chart_data.xaxis = d3.axisBottom(chart_data.xscale)
            .tickValues(set_tickvalues())
            .tickFormat(set_tickformat());
            //.ticks(5)
        // ordinal scale does not have ticks property so must set tick values explicitly
        //.tickValues(function(d) { return set_tickvalues(d);} )
            //.tickFormat(d3.format(".1")); //set_tickformat());
            //.tickSubdivide(0);
    }
    else if (plot.chart_type == "barchart")
    {
        chart_data.xaxis = d3.axisBottom(chart_data.xscale)
            .tickValues(set_tickvalues())
            .tickFormat(set_tickformat());
    }
    else
    {
        //var xscale = chart_data.xscale;
        chart_data.xaxis = d3.axisBottom(chart_data.xscale)
            .tickFormat(logFormatx); //d3.format("e"));
    }



    function set_tickvalues()
    {
        if (plot.scale_type == "log")
        {
            vals = chart_data.xscale.domain().filter(function(d) { return powerOfTen(d); });
        }
        else
        {
            // make array of multiples of 10^n where maxval is of the order of 10^(n+1)
            maxval = d3.max(chart_data.xscale.domain());

            if (maxval < 16)
            {
                vals = d3.range(0,maxval)
            }
            else
            {
				step = Math.pow(10, Math.round(Math.log10(maxval / 5)));

                vals = d3.range(0, maxval, step);

                while (vals.length > 6)
                {
                    step = step * 2;
                    vals = d3.range(0, maxval, step);
                }
            }
        }
        return vals;
    }

    function logFormatx(d)
    {
        return logFormat(d, chart_data.xscale.domain());
    }

    function logFormaty(d)
    {
        return logFormat(d, chart_data.yscale.domain());
    }

    function logFormat(d, axisrange) {

        var numberFormat = d3.format(".1e");

        // want to make sure we're always displaying some tickmarks
        if (Math.log10(axisrange[1]) - Math.log10(axisrange[0]) < 2)
        {
            // powers of ten won't work here
            // go for default unless there's going to be more than 6 labels,
            // then take every second label
            var step = Math.pow(10, Math.floor(Math.log10(axisrange[1])) - 1);
            if ((axisrange[1] - axisrange[0])/step > 6 )
            {
                return ((d / step) % 2 === 0) ? numberFormat(d) : "";
            }
            else
            {
                return(numberFormat(d));
            }
        }
        else
        {
            return powerOfTen(d) ? numberFormat(d) : "";
        }
    }

    function set_tickformat()
    {
        if (max_xval < 10001)
        {
            formatValue = d3.format("d");
        }
        else
        {
            formatValue = d3.format(".1e");
        }

        return formatValue;
    }

    /* Now draw x axis*/
    chart_data.svg.append("g")
    .attr("class", "x axis")
    .attr("transform", "translate(0, "+ plot_height +")")
    .call(chart_data.xaxis);


    if ((plot.scale_type == "log") && (plot.chart_type == "scatter"))
    {
        chart_data.yaxis = d3.axisLeft(chart_data.yscale)
        .tickFormat(logFormaty);
    }
    else
    {
        chart_data.yaxis = d3.axisLeft(chart_data.yscale)
        .ticks(numticks)
        .tickFormat(d3.format(".1e"));
    }



    /* Draw y axis */
    chart_data.svg.append("g")
    .attr("class", "y axis")
    .attr("transform", "translate("+ 0 +", 0)")
    .call(chart_data.yaxis);


    function powerOfTen(d) {
     return d / Math.pow(10, Math.ceil(Math.log10(d))) === 1;
    }

    xlabel = plot.xlabel;
    ylabel = plot.ylabel;
    
    // add the axes labels
    chart_data.svg.append("g")
    .attr("class", "axis-label")
    .append("text")
    .attr("text-anchor", "middle")
    .attr("x", plot_width/2)    // changes x location
    .attr("y", plot_height + bottom_margin/2) // changes y location
    .attr("dy", "1em")
    .text(xlabel);
    
    // y-axis label requires different handling
    chart_data.svg
    .append("g")
    .attr("class", "axis-label")
    .attr("transform", "translate(" + -left_margin*4/5 + ", " + plot_height/2 + ")") //position label
    .append("text")
    .attr("text-anchor", "middle")
    .attr("transform", "rotate(-90)") //rotate about mid-point
    .text(ylabel);
    
}

    return module;
}());



