jQuery.extend( jQuery.fn.dataTableExt.oSort, {
              
    "excludeNaN-asc": function ( a, b ) {
              
        var x = Number.MAX_SAFE_INTEGER;
        var y = Number.MAX_SAFE_INTEGER;
              
        if (!(isNaN(a) || (a==Infinity)))
        {
            x = parseFloat(a);
        }
        if (!(isNaN(b) || (b==Infinity)))
        {
            y = parseFloat(b);
        }
        return ((x < y) ? -1 : ((x > y) ?  1 : 0));
              
    },
              
    "excludeNaN-desc": function ( a, b ) {
              
        var x = -Number.MAX_SAFE_INTEGER;
        var y = -Number.MAX_SAFE_INTEGER;
              
        if (!(isNaN(a) || (a==Infinity)))
        {
            x = parseFloat(a);
        }
        if (!(isNaN(b) || (b==Infinity)))
        {
            y = parseFloat(b);
        }
        return ((x < y) ? 1 : ((x > y) ?  -1 : 0));
              
    }
});
