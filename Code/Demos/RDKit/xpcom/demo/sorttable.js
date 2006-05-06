ts_addEvent(window, "load", ts_sortables_init);

var SORT_COLUMN_INDEX;
var GLOBAL_ROWS;
var GLOBAL_TABLE_SIZE=3;
var GLOBAL_START_POINT=0;

function ts_sortables_init() {
  // Find all tables with class sortable and make them sortable
  if (!document.getElementsByTagName) return;
  tbls = document.getElementsByTagName("table");
  for (ti=0;ti<tbls.length;ti++) {
    thisTbl = tbls[ti];
    if (((' '+thisTbl.className+' ').indexOf("sortable") != -1) && (thisTbl.id)) {
      //initTable(thisTbl.id);
      ts_makeSortable(thisTbl);
    }
  }
}

function ts_getInnerText(el) {
  if (typeof el == "string") return el;
  if (typeof el == "undefined") { return el };
  if (el.innerText) return el.innerText;	//Not needed but it is faster
  var str = "";
	
  var cs = el.childNodes;
  var l = cs.length;
  for (var i = 0; i < l; i++) {
    switch (cs[i].nodeType) {
    case 1: //ELEMENT_NODE
      str += ts_getInnerText(cs[i]);
      break;
    case 3:	//TEXT_NODE
      str += cs[i].nodeValue;
      break;
    }
  }
  return str;
}


function ts_makeSortable(table) {
  var firstRow;	
  if (table.rows && table.rows.length > 0) {
    firstRow = table.rows[0];
  }
  if (!firstRow) return;
    
  // We have a first row: assume it's the header, and make its contents clickable links
  for (var i=0;i<firstRow.cells.length;i++) {
    var cell = firstRow.cells[i];
    var txt = ts_getInnerText(cell);
    cell.innerHTML = '<a href="#" class="sortheader" onclick="ts_sortAndRefreshTable(this);return false;">'+txt+'<span class="sortarrow">&nbsp;&nbsp;&nbsp;</span></a>';
  }

  var newRows = new Array();
  var j;
  for (j=1;j<table.rows.length;j++) { newRows[j-1] = table.rows[j]; }
  GLOBAL_ROWS=newRows;
  ts_refreshTable(table,GLOBAL_START_POINT,GLOBAL_TABLE_SIZE);
}

function ts_clearArrows(table){
  // Delete any other arrows there may be showing
  var allspans = document.getElementsByTagName("span");
  var ci;
  for (ci=0;ci<allspans.length;ci++) {
    if (allspans[ci].className == 'sortarrow') {
      if (ts_getParent(allspans[ci],"table") == table) { // in the same table as us?
	allspans[ci].innerHTML = '&nbsp;&nbsp;&nbsp;';
      }
    }
  }
}

function ts_resortTable(lnk) {
  // get the span
  var span;
  var ci;
  for (ci=0;ci<lnk.childNodes.length;ci++) {
    if (lnk.childNodes[ci].tagName &&
	lnk.childNodes[ci].tagName.toLowerCase() == 'span'){
      span = lnk.childNodes[ci];
    }
  }
  var spantext = ts_getInnerText(span);
  var td = lnk.parentNode;
  var column = td.cellIndex;
  var table = ts_getParent(td,'TABLE');
    
  // Work out a type for the column
  if (!GLOBAL_ROWS || GLOBAL_ROWS.length <= 1) return;
  var itm = String(GLOBAL_ROWS[0][column]);
  var sortfn = ts_sort_caseinsensitive;
  if (itm.match(/^\d\d[\/-]\d\d[\/-]\d\d\d\d$/)) sortfn = ts_sort_date;
  if (itm.match(/^\d\d[\/-]\d\d[\/-]\d\d$/)) sortfn = ts_sort_date;
  if (itm.match(/^[£$]/)) sortfn = ts_sort_currency;
  if (itm.match(/^[-]?[\d\.]+$/)) sortfn = ts_sort_numeric;
  SORT_COLUMN_INDEX = column;

  var i,j;
  GLOBAL_ROWS.sort(sortfn);
  var ARROW;
  if (span.getAttribute("sortdir") == 'down') {
    ARROW = '&nbsp;&nbsp;&uarr;';
    GLOBAL_ROWS.reverse();
    span.setAttribute('sortdir','up');
  } else {
    ARROW = '&nbsp;&nbsp;&darr;';
    span.setAttribute('sortdir','down');
  }

  // Delete any other arrows there may be showing
  ts_clearArrows(table);
  span.innerHTML = ARROW;
}


function ts_refreshTable(table,startAt,numberToUse){
  var newRows;
  var i,j;
  if(!GLOBAL_ROWS){
    return;
  }
  if(startAt>GLOBAL_ROWS.length)
    return;
  rows = GLOBAL_ROWS;

  if(numberToUse <=0) numberToUse=rows.length;
  var lastEntry;
  lastEntry=startAt+numberToUse;
  if(lastEntry>=rows.length)
    lastEntry=rows.length;
    
  var offset=1;
  for (i=startAt;i<lastEntry;i++) {
    while(table.rows.length<=offset){
      var tmpRow=table.insertRow(table.rows.length);
      for(j=0;j<rows[i].length;j++){
	tmpRow.insertCell(j);
      }
    }
    var row=table.rows[offset];
    offset++;
    for(j=0;j<rows[i].length;j++){
      row.cells[j].innerHTML=rows[i][j];
    }
  }
  while(offset<table.rows.length){
    row = table.rows[offset++];
    for(j=0;j<rows[0].length;j++){
      row.cells[j].innerHTML="";
    }
  }
  button = document.getElementById('prevButton');
  if(startAt<=0){
    button.disabled=1;
  } else {
    button.disabled=0;
  }
  button = document.getElementById('nextButton');
  if(lastEntry>=rows.length){
    button.disabled=1;
  } else {
    button.disabled=0;
  }

}

function ts_refresh(lnk,startAt,numberToUse){
  var td = lnk.parentNode;
  var table = ts_getParent(td,'TABLE');
  ts_refreshTable(table,startAt,numberToUse);
}


function ts_sortAndRefreshTable(lnk){
  ts_resortTable(lnk);
  ts_refresh(lnk,GLOBAL_START_POINT,GLOBAL_TABLE_SIZE);
}

function ts_incrementTable(id){
  var table = document.getElementById(id);
  if(!table) return;
  GLOBAL_START_POINT += GLOBAL_TABLE_SIZE;
  ts_refreshTable(table,GLOBAL_START_POINT,GLOBAL_TABLE_SIZE);
}

function ts_decrementTable(id){
  var table = document.getElementById(id);
  if(!table) return;
  GLOBAL_START_POINT -= GLOBAL_TABLE_SIZE;
  if(GLOBAL_START_POINT<0) GLOBAL_START_POINT=0;
  ts_refreshTable(table,GLOBAL_START_POINT,GLOBAL_TABLE_SIZE);
}

function ts_getParent(el, pTagName) {
  if (el == null) return null;
  else if (el.nodeType == 1 && el.tagName.toLowerCase() == pTagName.toLowerCase()) // Gecko bug, supposed to be uppercase
    return el;
  else
    return ts_getParent(el.parentNode, pTagName);
}
function ts_sort_date(a,b) {
  // y2k notes: two digit years less than 50 are treated as 20XX,
  // greater than 50 are treated as 19XX
  aa = String(a[SORT_COLUMN_INDEX]);
  bb = String(b[SORT_COLUMN_INDEX]);
  if (aa.length == 10) {
    dt1 = aa.substr(6,4)+aa.substr(3,2)+aa.substr(0,2);
  } else {
    yr = aa.substr(6,2);
    if (parseInt(yr) < 50) { yr = '20'+yr; } else { yr = '19'+yr; }
    dt1 = yr+aa.substr(3,2)+aa.substr(0,2);
  }
  if (bb.length == 10) {
    dt2 = bb.substr(6,4)+bb.substr(3,2)+bb.substr(0,2);
  } else {
    yr = bb.substr(6,2);
    if (parseInt(yr) < 50) { yr = '20'+yr; } else { yr = '19'+yr; }
    dt2 = yr+bb.substr(3,2)+bb.substr(0,2);
  }
  if (dt1==dt2) return 0;
  if (dt1<dt2) return -1;
  return 1;
}

function ts_sort_currency(a,b) { 
  aa = String(a[SORT_COLUMN_INDEX]).replace(/[^0-9.]/g,'');
  bb = String(b[SORT_COLUMN_INDEX]).replace(/[^0-9.]/g,'');
  return parseFloat(aa) - parseFloat(bb);
}

function ts_sort_numeric(a,b) { 
  var aa = a[SORT_COLUMN_INDEX];
  if (isNaN(aa)) aa = 0;
  var bb = b[SORT_COLUMN_INDEX]; 
  if (isNaN(bb)) bb = 0;
  return aa-bb;
}

function ts_sort_caseinsensitive(a,b) {
  var aa = String(a[SORT_COLUMN_INDEX]).toLowerCase();
  var bb = String(b[SORT_COLUMN_INDEX]).toLowerCase();
  if (aa==bb) return 0;
  if (aa<bb) return -1;
  return 1;
}

function ts_sort_default(a,b) {
  var aa = a[SORT_COLUMN_INDEX];
  var bb = b[SORT_COLUMN_INDEX];
  if (aa==bb) return 0;
  if (aa<bb) return -1;
  return 1;
}


function ts_addEvent(elm, evType, fn, useCapture)
     // addEvent and removeEvent
     // cross-browser event handling for IE5+,  NS6 and Mozilla
     // By Scott Andrew
{
  if (elm.addEventListener){
    elm.addEventListener(evType, fn, useCapture);
    return true;
  } else if (elm.attachEvent){
    var r = elm.attachEvent("on"+evType, fn);
    return r;
  } else {
    alert("Handler could not be removed");
  }
  return false;
};


