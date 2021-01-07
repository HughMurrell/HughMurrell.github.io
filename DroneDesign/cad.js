//render a simple house using ConnectTheDots and Square

function draw(){
var makerjs = require('makerjs');

    var R = parseFloat(document.getElementById("R").value);
    var D = parseFloat(document.getElementById("D").value);
    var Z = parseFloat(document.getElementById("Z").value);
    var X = parseFloat(document.getElementById("X").value);
    var Y = parseFloat(document.getElementById("Y").value);
    var P = parseFloat(document.getElementById("P").value);
    var Q = parseFloat(document.getElementById("Q").value);
    
    P = P*Math.PI/180;
    Q = Q*Math.PI/180;
    
var cp=R*Math.cos(P);
var sp=R*Math.sin(P);
var cq=R*Math.cos(Q);
var sq=R*Math.sin(Q);
    
console.log(R, " ", Z, " ", X, " ", Y, " ", P, " " , Q);
console.log(cp, " ", sp, " ", cq, " ", sq);
    
var top_dots = [
            [0, X + sp],
            [cp,X],
            [cp,0],
            [cp+Z,0],
            [cp+Z,Y],
            [cp+Z+cq,Y+sq]
                ];

console.log(top_dots);

var top_drone = new makerjs.models.ConnectTheDots(true, top_dots);
                
var bot_dots = [
        [0, -(X + sp)],
        [cp,-X],
        [cp,0],
        [cp+Z,0],
        [cp+Z,-Y],
        [cp+Z+cq,-(Y+sq)]
                ];

var bot_drone = new makerjs.models.ConnectTheDots(true, bot_dots);
    
    
var rotors = {
    paths: {
        "ft": new makerjs.paths.Arc([0, X + sp], D/2, 0, 360),
        "fb": new makerjs.paths.Arc([0, -(X + sp)], D/2, 0, 360),
        "bt": new makerjs.paths.Arc([cp+Z+cq,Y+sq], D/2, 0, 360),
        "bb": new makerjs.paths.Arc([cp+Z+cq,-(Y+sq)], D/2, 0, 360)
    }
};

var window1 = new makerjs.models.Square(40);
window1.origin = [40, 40];

var window2 = new makerjs.models.Square(40);
window2.origin = [-80, 40];

var rot=(Math.PI-P)*180/Math.PI;
var arm_ft = new makerjs.models.Oval(R+20, 20);
makerjs.model.move(arm_ft,[cp-10,X-10]);
makerjs.model.rotate(arm_ft,rot,[cp,X]);
                                    
var rot=(Q)*180/Math.PI;
var arm_bt = new makerjs.models.Oval(R+20, 20);
makerjs.model.move(arm_bt,[cp+Z-10,Y-10]);
makerjs.model.rotate(arm_bt,rot,[cp+Z,Y]);
                                    
var rot=(Math.PI+P)*180/Math.PI;
var arm_fb = new makerjs.models.Oval(R+20, 20);
makerjs.model.move(arm_fb,[cp-10,-X-10]);
makerjs.model.rotate(arm_fb,rot,[cp,-X]);

var rot=(-Q)*180/Math.PI;
var arm_bb = new makerjs.models.Oval(R+20, 20);
makerjs.model.move(arm_bb,[cp+Z-10,-Y-10]);
makerjs.model.rotate(arm_bb,rot,[cp+Z,-Y]);
                                    
var drone = {
  models: {
      "top": top_drone,
      "bottom": bot_drone,
      "rotor blades": rotors,
      "arm_ft": arm_ft,
      "arm_bt": arm_bt,
      "arm_fb": arm_fb,
      "arm_bb": arm_bb
                                    }
};
    
var svg = makerjs.exporter.toSVG(drone);
    
var frog = window.open("","wildebeast","width=900,height=900,scrollbars=1,resizable=1")

    frog.document.open()
    frog.document.write(svg)
    frog.document.close()

}
