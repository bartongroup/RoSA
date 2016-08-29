/***********************************************************************************/
// Definitions of scatter plots
/***********************************************************************************/

// chart_type : "scatter" or "histogram"
// scale_type : "linear" or "log"
// xlabel : label for x axis on plot
// ylabel : label for y axis on plots
// suffix : reference for suffix (e.g. FPKM) to add to data column name to get alternative values
// xname : name of data column for x values
// yname : name of data column for y values
// xdisplay : name of data column for tooltip x values (not always same as data column!)
// ydisplay : name of data column for tooltip y values (not always same as data column!)
// title : title of plot, goes in panel header
// index : name of data column containing gene indices
// id : document element id for this plot
// parent_id : document element id of parent container (as id won't exist until created on load)
// json : name of json file containing data for this plot
// buttons : list of button definitions for buttons to go in panel header

var svg_sense_scatter =
{
    config :
        {
            chart_type : "scatter",
            scale_type : "linear",
            xlabel : "Number of sense counts",
            ylabel : "Number of antisense counts",
            suffix : "",
            xname : "exon",
            yname : "antisense",
            xdisplay : "exon",
            ydisplay : "antisense",
            index : "gene_id",
            id : "scatterplot",
            data : null
        },

    title : "Sense versus anti-sense counts",
    parent_id : "anti-parent",
    json : "genes.json",
    buttons : [ {button_label: "Raw counts", buttonset: "Norm",
                  button_list: [{id: "anti_counts_raw", label: "Raw counts", norm: "raw"},
                                {id: "anti_counts_fpkm", label: "FPKM", norm: "fpkm"}]},
               {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "anti_counts_linear", label: "Linear scale", scale: "linear"},
                                {id: "anti_counts_log", label: "Log scale", scale: "log"}]}]
};
var svg_intron_exon_scatter =
{
    config :
    {
        chart_type : "scatter",
        scale_type : "linear",
        xlabel : "Number of counts in exons",
        ylabel : "Number of counts in introns",
        suffix : "",
        xname : "exon",
        yname : "intron",
        xdisplay : "exon",
        ydisplay : "intron",
        index : "gene_id",
        id : "intron-exon-scatter",
        data : null
    },
    title : "Intron counts versus exon counts",
    parent_id : "intron-exon-parent",
    json : "genes.json",
    buttons : [ {button_label: "Raw counts", buttonset: "Norm", button_list: [{id: "intron_counts_raw", label: "Raw counts", norm: "raw"},
                                                           {id: "intron_counts_fpkm", label: "FPKM", norm: "fpkm"}]},
               {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "intron_counts_linear", label: "Linear scale", scale: "linear"},
                                {id: "intron_counts_log", label: "Log scale", scale: "log"}]}]
};

var svg_tutr_exon_scatter =
{
    config :
    {
        chart_type : "scatter",
        scale_type : "linear",
        xlabel : "Number of counts in exons",
        ylabel : "Number of counts in 3'UTRs",
        suffix : "",
        xname : "exon",
        yname : "three_prime_UTR",
        xdisplay : "exon",
        ydisplay : "three_prime_UTR",
        index : "gene_id",
        id : "tutr-exon-scatter",
        data : null
    },
    title : "3'UTR counts versus exon counts",
    parent_id : "tutr-parent",
    json : "genes.json",
    buttons : [ {button_label: "Raw counts", buttonset: "Norm", button_list: [{id: "tutr_counts_raw", label: "Raw counts", norm: "raw"},
                                                           {id: "tutr_counts_fpkm", label: "FPKM", norm: "fpkm"}]},
               {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "tutr_counts_linear", label: "Linear scale", scale: "linear"},
                                {id: "tutr_counts_log", label: "Log scale", scale: "log"}]}]
};
var svg_futr_exon_scatter =
{
    config :
    {
        chart_type : "scatter",
        scale_type : "linear",
        xlabel : "Number of counts in exons",
        ylabel : "Number of counts in 5'UTRs",
        suffix : "",
        xname : "exon",
        yname : "five_prime_UTR",
        xdisplay : "exon",
        ydisplay : "five_prime_UTR",
        index : "gene_id",
        id : "futr-exon-scatter",
        data : null
    },
    title : "5'UTR counts versus exon counts",
    parent_id : "futr-parent",
    json : "genes.json",
    buttons : [ {button_label: "Raw counts", buttonset: "Norm", button_list: [{id: "futr_counts_raw", label: "Raw counts", norm: "raw"},
                                                           {id: "futr_counts_fpkm", label: "FPKM", norm: "fpkm"}]},
               {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "futr_counts_linear", label: "Linear scale", scale: "linear"},
                                {id: "futr_counts_log", label: "Log scale", scale: "log"}]}]
};

var svg_tufu_scatter =
{
    config :
    {
        chart_type : "scatter",
        scale_type : "linear",
        xlabel : "Number of counts in 3'UTRs",
        ylabel : "Number of counts in 5'UTRs",
        suffix : "",
        xname : "three_prime_UTR",
        yname : "five_prime_UTR",
        xdisplay : "three_prime_UTR",
        ydisplay : "five_prime_UTR",
        index : "gene_id",
        id : "tufu-scatter",
        data : null
    },
    title : "5' UTR versus 3' UTR counts",
    parent_id : "tufu-parent",
    json : "genes.json",
    buttons : [ {button_label: "Raw counts", buttonset: "Norm", button_list: [{id: "tufu_counts_raw", label: "Raw counts", norm: "raw"},
                                                           {id: "tufu_counts_fpkm", label: "FPKM", norm: "fpkm"}]},
               {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "tufu_counts_linear", label: "Linear scale", scale: "linear"},
                                {id: "tufu_counts_log", label: "Log scale", scale: "log"}]}]
    
};

/***********************************************************************************/
// Definitions of histograms
/***********************************************************************************/
var svg_hist_exon =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 100,
        ylabel : "Number of exons",
        xlabel : "Exon length",
        suffix : "",
        xname : "length",
        index : "gene_id",
        id : "exon-lengths-hist",
        data : null
    },
    title : "Exon lengths distribution",
    parent_id : "exon_lengths",
    json : "exon-data.json",
    buttons : [ {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "exon_lengths_linear", label: "Linear scale", scale: "linear"},
                                {id: "exon_lengths_log", label: "Log scale", scale: "log"}]}]
};
var svg_hist_intron =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 100,
        ylabel : "Number of introns",
        xlabel : "Intron length",
        suffix : "",
        xname : "length",
        index : "gene_id",
        id : "intron-lengths-hist",
        data : null
    },
    title : "Intron lengths distribution",
    parent_id : "intron_lengths",
    json : "intron-data.json",
    buttons : [ {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "intron_lengths_linear", label: "Linear scale", scale: "linear"},
                                {id: "intron_lengths_log", label: "Log scale", scale: "log"}]}]
};
var svg_exon_count =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 100,
        ylabel : "Number of genes",
        xlabel : "Number of exons/gene",
        suffix : "",
        xname : "exon_counts",
        index : "gene_id",
        id : "exon-counts-hist",
        data : null
    },
    title : "Exon counts distribution",
    parent_id : "exon_counts",
    json : "genes.json",
    buttons : [ {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "exon_counts_linear", label: "Linear scale", scale: "linear"},
                                {id: "exon_counts_log", label: "Log scale", scale: "log"}]}]
};

var svg_tr_count =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 10,
        ylabel : "Number of genes",
        xlabel : "Number of transcripts/gene",
        suffix : "",
        xname : "transcript_counts",
        index : "gene_id",
        id : "tr-counts-hist",
        data : null
    },
    title : "Transcript counts distribution",
    parent_id : "tr_counts",
    json : "genes.json",
    buttons : []
};
var svg_gene_exons =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 10,
        ylabel : "Number of exons",
        xlabel : "Exon length",
        suffix : "",
        xname : "length",
        index : "gene_id",
        id : "gene_exon_lengths",
        data : null
    },
    title : "Exon lengths for selected gene",
    parent_id : "subchart1",
    json : "exon-data.json",
    buttons : []
};
var svg_gene_introns =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 10,
        ylabel : "Number of introns",
        xlabel : "Intron length",
        suffix : "",
        xname : "length",
        index : "gene_id",
        id : "gene_intron_lengths",
        data : null
    },
    title : "Intron lengths for selected gene",
    parent_id : "subchart2",
    json : "intron-data.json",
    buttons : []
};
var svg_tutr_lengths =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 100,
        ylabel : "Number of 3' UTRs",
        xlabel : "3' UTR length",
        suffix : "",
        xname : "length",
        index : "gene_id",
        id : "tutr_lengths",
        data : null
    },
    title : "3' UTR lengths distribution",
    parent_id : "tutr_lengths_parent",
    json : "three_prime_UTR-data.json",
    buttons : [ {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "tutr_lengths_linear", label: "Linear scale", scale: "linear"},
                                {id: "tutr_lengths_log", label: "Log scale", scale: "log"}]}]
};

var svg_futr_lengths =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 100,
        ylabel : "Number of 5' UTRs",
        xlabel : "5' UTR length",
        suffix : "",
        xname : "length",
        index : "gene_id",
        id : "futr_lengths",
        data : null
    },
    title : "5' UTR lengths distribution",
    parent_id : "futr_lengths_parent",
    json : "five_prime_UTR-data.json",
    buttons : [ {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "futr_lengths_linear", label: "Linear scale", scale: "linear"},
                                {id: "futr_lengths_log", label: "Log scale", scale: "log"}]}]
};

var svg_tr_lengths =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 100,
        ylabel : "Number of genes",
        xlabel : "Mean transcript length",
        suffix : "",
        xname : "length",
        index : "gene_id",
        id : "tr_lengths",
        data : null
    },
    title : "Gene lengths distribution",
    parent_id : "tr_lengths_parent",
    json : "genes.json",
    buttons : [ {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "tr_lengths_linear", label: "Linear scale", scale: "linear"},
                                {id: "tr_lengths_log", label: "Log scale", scale: "log"}]}]
};

var svg_first_exons =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 100,
        ylabel : "Number of genes",
        xlabel : "Counts in first exon",
        suffix : "",
        xname : "counts",
        index : "gene_id",
        id : "first_exons",
        data : null
    },
    title : "Distribution of counts in first exons",
    parent_id : "first_exons_parent",
    json : "firstexons.json",
    buttons : [ {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "first_exons_linear", label: "Linear scale", scale: "linear"},
                                {id: "first_exons_log", label: "Log scale", scale: "log"}]}]
};

var svg_last_exons =
{
    config :
    {
        chart_type : "histogram",
        scale_type : "linear",
        num_bins : 100,
        ylabel : "Number of genes",
        xlabel : "Counts in last exon",
        suffix : "",
        xname : "counts",
        index : "gene_id",
        id : "last_exons",
        data : null
    },
    title : "Distribution of counts in last exons",
    parent_id : "last_exons_parent",
    json : "lastexons.json",
    buttons : [ {button_label: "Linear scale", buttonset: "Scale",
                  button_list: [{id: "last_exons_linear", label: "Linear scale", scale: "linear"},
                                {id: "last_exons_log", label: "Log scale", scale: "log"}]}]
};
/***********************************************************************************/
// Definitions of barcharts
/***********************************************************************************/
var svg_gene_len_by_exon =
{
    config :
    {
        chart_type : "barchart",
        scale_type : "linear",
        num_bins : 10,
        ylabel : "Length",
        xlabel : "Exon",
        suffix : "",
        yname : "length",
        index : "gene_id",
        id : "gene_len_by-exon",
        data : null
    },
    title : "Exon lengths for selected gene",
    parent_id : "subchart3",
    json : "exon-data.json",
    buttons : []
};

var svg_annot_totals =
{
    config :
    {
        chart_type : "barchart",
        xname : "name",
        yname : "value",
        ydisplay : "value_text",
        xdisplay : "name",
        ylabel : "Total number of features in annotation",
        xlabel : "",
        id: "annot_totals",
        data : null
    },
    plot: null,
    parent_id: "annot_totals_parent",
    title: "Annotation totals",
    buttons: []
};

var svg_read_totals_fwd =
{
    config :
    {
        svg: null,
        xscale : null,
        yscale : null,
        xaxis : null,
        yaxis : null,
        chart_type : "barchart",
        xname : "name",
        yname : "value",
        ydisplay : "value_text",
        xdisplay : "name",
        ylabel : "Total number of reads",
        xlabel : "",
        id: "read_totals_fwd",
        data : null,
    },
    plot: null,
    parent_id: "read_totals_fwd_parent",
    title: "Forward strand read totals",
    buttons: []
};

var svg_read_totals_rev =
{
    config :
    {
        chart_type : "barchart",
        xname : "name",
        yname : "value",
        ydisplay : "value_text",
        xdisplay : "name",
        ylabel : "Total number of reads",
        xlabel : "",
        id: "read_totals_rev",
        data : null
    },
    plot: null,
    parent_id: "read_totals_rev_parent",
    title: "Reverse strand read totals",
    buttons: []
};