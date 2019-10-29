// $Id$
//
//  Copyright (C) 2009 Greg Landrum
//
//   @@ All Rights Reserved  @@
//
var LINE=1;
var WEDGE=2;
var ATOM=3;
var BOUNDS=4;
var RESOLUTION=5;
var END=-1;

var elemDict={
  6:"rgb(0,0,0)",
  7:"rgb(0,0,255)",
  8:"rgb(255,0,0)",
  9:"rgb(51,204,204)",
  15:"rgb(255,128,0)",
  16:"rgb(204,204,0)",
  17:"rgb(0,204,0)",
  35:"rgb(128,77,26)",
  0:"rgb(128,128,128)",
}

function DrawMol(canvasId,molData){
  if(molData.length<7){
    alert("data too short");
    return;
  }
  var canvas = document.getElementById(canvasId);
  var ctx = canvas.getContext('2d');
  ctx.restore();
  ctx.clearRect(0,0,canvas.width,canvas.height); // clear canvas

  if(molData[0]!=RESOLUTION || molData[1]<=0) {
    alert("no resolution data!")
    return;
  }
  var dotsPerAngstrom=molData[1];

  if(molData[2]!=BOUNDS ){
    alert("no bounds data!")
    return;
  }
  var minx=molData[3];
  var miny=molData[4];
  var maxx=molData[5];
  var maxy=molData[6];

  // we use the canvas to do our scaling and centering:
  var xscale = canvas.width/(maxx-minx);
  var yscale = canvas.height/(maxy-miny);
  var cx,cy;
  if(xscale>1){
    // the mol is too small in the x direction to fill the canvas
    var fraction = 1./xscale;
    var unused = canvas.width*(1-fraction);
    cx=unused/2-minx
    xscale=1;
  } else {
    cx=-minx;
  }
  if(yscale>1){
    // the mol is too small in the x direction to fill the canvas
    var fraction = 1./yscale;
    var unused = canvas.height*(1-fraction);
    cy=unused/2-miny
    yscale=1;
  } else {
    cy=-miny;
  }

  var scale;
  if(xscale<yscale){
    scale=xscale;
  } else {
    scale=yscale;
  }

  
  ctx.save();
  ctx.translate(cx,cy);
  ctx.scale(scale,scale);
  ctx.font = "2em 'arial'";

  var idx=7;
  while(idx<molData.length){
    var tag = molData[idx++];
    if(tag==LINE){
      var width=molData[idx++];
      ctx.lineWidth=width;
      var dashed=molData[idx++];
      var n1=molData[idx++];
      var c1 = elemDict[n1];
      if(!c1) c1=elemDict[6];
      var n2=molData[idx++];
      var c2 = elemDict[n2];
      if(!c2) c2=elemDict[6];
      var x1=molData[idx++];
      var y1=molData[idx++];
      var x2=molData[idx++];
      var y2=molData[idx++];
      if(c1==c2){
        ctx.beginPath();
        ctx.strokeStyle=c1;
        ctx.moveTo(x1,y1);
        ctx.lineTo(x2,y2);
        ctx.closePath();
        ctx.stroke()
      } else {
        var mpx=x1+(x2-x1)/2;
        var mpy=y1+(y2-y1)/2;

        ctx.beginPath();
        ctx.strokeStyle=c1;
        ctx.moveTo(x1,y1);
        ctx.lineTo(mpx,mpy);
        ctx.closePath();
        ctx.stroke()

        ctx.beginPath();
        ctx.strokeStyle=c2;
        ctx.moveTo(mpx,mpy);
        ctx.lineTo(x2,y2);
        ctx.closePath();
        ctx.stroke()

      }
    } else if(tag==WEDGE){
      
    } else if(tag==ATOM){
      var n1=molData[idx++];
      var x1=molData[idx++];
      var y1=molData[idx++];
      var nChars=molData[idx++];
      var symb="";
      for(j=0;j<nChars;++j){
          symb+=String.fromCharCode(molData[idx++]);
      }
      
      // not all canvas implementations have fillText
      if(ctx.fillText){
          if(symb!="C"){
              var c1 = elemDict[n1];
              if(!c1) c1=elemDict[6];
              ctx.textAlign='center';
              ctx.textBaseline='middle';
              measure = ctx.measureText(symb);
              ctx.beginPath();
              ctx.fillStyle="rgb(255,255,255)";
              ctx.arc(x1,y1,measure.width,0,6.3,0);
              ctx.closePath();
              ctx.fill();
              ctx.beginPath();
              ctx.fillStyle=c1;
              ctx.fillText(symb,x1,y1);
              ctx.closePath();
              ctx.fill();
          }
      }
    } else if(tag==END){
        break;
    } else {
      alert("bad tag!");
      return;
    }
  }

}
