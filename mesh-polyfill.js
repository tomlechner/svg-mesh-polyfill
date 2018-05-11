
// Polyfill for svg meshes.
// Copyright Tom Lechner, 2018
// Distributed under GNU General Public License version 3 or later. See <http://fsf.org/>.


(function() {

    // Namespaces -----------------------------------
    var svgNS    = "http://www.w3.org/2000/svg";
    var xlinkNS  = "http://www.w3.org/1999/xlink"

    // Test if mesh gradients are supported.
    var m = document.createElementNS( svgNS, "meshgradient" );
    if (m.x) {
        return;
    }




    //-------------------------- math functions

    /*! The bezier matrix.
     *  See <a href="lotsamath.html">here</a> for what this is.
     */
    var B = [	-1,  3, -3,  1,
                 3, -6,  3,  0,
                -3,  3,  0,  0,
                 1,  0,  0,  0
            ];

    /*! The inverse of the bezier matrix.
     *  See <a href="lotsamath.html">here</a> for what this is.
     */
    var Binv = [ 0,    0,    0,   1,
                 0,    0,  1./3,  1,
                 0,  1./3, 2./3,  1,
                 1,    1,    1,   1
               ];

    /*! v=[ t^3, t^2, t, 1 ]
     */
    function getT(t) 
    {
        return [ t*t*t, t*t, t, 1 ];
    }

    /*! Returns a*b, a and b are double[4]..
     */
    return dot(a, b)
    {
        var f=0;
        for (var c=0; c<4; c++) f += a[c]*b[c];
        return f;
    }

    /*! return v = m b (b is the vector), m is 4x4 matrix, b is double[4]
     */
    function m_times_v(m, b)
    {
        v = [];

        for (var r=0; r<4; r++) {
            v[r]=0;
            for (var c=0; c<4; c++) v[r] += m[r*4+c]*b[c];
        }

        return v;
    }

    /*! return m = a x b, a and b are 4x4 matrices
     */
    function m_times_m(a, b)
    {
        var m = [];
        for (var r=0; r<4; r++) {
            for (var c=0; c<4; c++) {
                m[r*4 + c] = 0;
                for (var c2=0; c2<4; c2++) m[r*4 + c] += a[r*4 + c2] * b[c2*4 + c];
            }
        }

        return m;
    }

    /*! Fill I with 4x4 Identity matrix
     */
    function getI()
    {
        var I = [];
        for (var c=0; c<16; c++) I[c]=0;
        I[0]=I[5]=I[10]=I[15]=1;
        return I;
    }

    /*! Fill 4x4 I with x,y,z,w scaled by a,b,c,d
     */
    function getScaledI(a, b, c, d)
    {
        for (var c1=0; c1<16; c1++) I[c1]=0;
        I[0]=a;
        I[5]=b;
        I[10]=c;
        I[15]=d;
        return I;
    }

    /*! For debugging: cout a 4x4 matrix G[16].
     */
    function printG(const char *ch,double *G)
    {
        for (r=0; r<4; r++) {
            console.log(G[r*4+0], G[r*4+1], G[r*4+2], G[r*4+3]);
        }
    }

    /*! Transpose the 4x4 matrix in place.
     */
    function m_transpose(M)
    {
        var t;
        for (var r=0; r<4; r++) {
            for (var c=0; c<4; c++) {
                if (r <= c) continue;
                t = M[r*4 + c];
                M[r*4 + c] = M[c*4 + r];
                M[c*4 + r] = t;
            }
        }
    }


    /*! This makes a polynomial column vector T:
     *  <pre>
     *    [ t^3 ]
     *  T=[ t^2 ], v=n*(t-t0), t=v/n + t0, 
     *    [  t  ]
     *    [  1  ]
     * 
     *    V = N T
     * 
     *     [ 1/n^3  3*t0/n^2  3*t0^2/n   t0^3 ] 
     *   N=[   0     1/n^2     2*t0/n    t0^2 ]
     *     [   0       0         1/n      t0  ]
     *     [   0       0          0        1  ]
     * </pre>  
     *
     */
    function getPolyT(n, t0)
    {
        return  [ 1./n/n/n, 3*t0/n/n, 3*t0*t0/n, t0*t0*t0,
                    0,       1./n/n,    2*t0/n,    t0*t0,
                    0,       0,         1./n,      t0,
                    0        0,         0,         1
                ];
    }


    // Flatvector class -----------------------------------
	function Flatvector(xx,yy) {
		if (xx === undefined) this.x = 0; else this.x = xx;
		if (yy === undefined) this.y = 0; else this.y = yy;
	};

	Flatvector.prototype.Invert    = function() { this.x = -this.x; this.y = -this.y; }
	Flatvector.prototype.Normalize = function() {
			var d = this.x*this.x + this.y*this.y;
			if (d == 0) return;
			d = Math.sqrt(d);
			this.x /= d;
			this.y /= d;
		}
	Flatvector.prototype.Set        = function(xx,yy) { this.x = xx; this.y=yy; }
	Flatvector.prototype.isZero     = function() { return this.x==0 && this.y==0; }
	Flatvector.prototype.angle      = function() { return Math.atan2(this.y,this.x); }
	Flatvector.prototype.norm       = function() { return Math.sqrt(this.x*this.x+this.y*this.y); }
	Flatvector.prototype.norm2      = function() { return this.x*this.x+this.y*this.y; }
	Flatvector.prototype.cross      = function(v) { return this.x*v.y-this.y*v.x; } /*magnitude and sign of cross product, which points in z direction */
	Flatvector.prototype.distanceTo = function(v) { return Math.sqrt((v.x-this.x)*(v.x-this.x)+(v.y-this.y)*(v.y-this.y)); }
	Flatvector.prototype.add        = function(v) { return new Flatvector(this.x + v.x, this.y + v.y); }
	Flatvector.prototype.subtract   = function(v) { return new Flatvector(this.x - v.x, this.y - v.y); }
	Flatvector.prototype.dot        = function(v) { return this.x*v.x + this.y*v.y; }
	Flatvector.prototype.scale      = function(s) { return new Flatvector(this.x*s, this.y*s); }


	 //parse a list of numbers into an array of Flatvectors
    function parse_points(s) {
        var points = [],
            values = s.split(/[ ,]+/),
            i;
        for (i = 0; i < values.length-1; i += 2) {
            points.push( new Flatvector( parseFloat( values[i]), parseFloat( values[i+1] )));
        }

        return points;
    }


    // Affine class -----------------------------------

    // As defined in the SVG spec
    // | a  c  e |
    // | b  d  f |
    // | 0  0  1 |
    function Affine(a, b, c, d, e, f) {
        if (a === undefined) {
            this.a = 1;
            this.b = 0;
            this.c = 0;
            this.d = 1;
            this.e = 0;
            this.f = 0;
        } else {
            this.a = a;
            this.b = b;
            this.c = c;
            this.d = d;
            this.e = e;
            this.f = f;
        }
    }

    Affine.prototype.a = null;
    Affine.prototype.b = null;
    Affine.prototype.c = null;
    Affine.prototype.d = null;
    Affine.prototype.e = null;
    Affine.prototype.f = null;

    Affine.prototype.append = function(v) {
        if (!(v instanceof Affine)) {
            console.log ( "mesh.js: argument to Affine.append is not affine!");
        }

        var a = this.a * v.a + this.c * v.b,
            b = this.b * v.a + this.d * v.b,
            c = this.a * v.c + this.c * v.d,
            d = this.b * v.c + this.d * v.d,
            e = this.a * v.e + this.c * v.f + this.e,
            f = this.b * v.e + this.d * v.f + this.f;

        return new Affine(a, b, c, d, e, f);
    };

    Affine.prototype.toString = function() {
        return ("affine: "   + this.a + " " + this.c + " " + this.e +
            "\n        " + this.b + " " + this.d + " " + this.f);
    };



    // Utility functions ---------------------------------

    // Browsers return a string rather than a transform list for gradientTransform!
    function parseTransform(t) {
        var affine = new Affine(),
            radian,
            tan,
            cos;

        for (var i in t = t.match(/(\w+\(\s*(\-?\d+\.?\d*e?\-?\d*\s*,?\s*)+\))+/g)) {
            var c = t[i].match(/[\w\.\-]+/g),
                type = c.shift();

            switch (type) {
                case "translate":
                    var trans;
                    if (c.length == 2) {
                        trans = new Affine( 1, 0, 0, 1, c[0], c[1] );
                    } else {
                        console.log( "mesh.js: translate does not have 2 arguments!" );
                        trans = new Affine( 1, 0, 0, 1, 0, 0 );
                    }
                    affine = affine.append( trans );
                    break;

                case "scale":
                    var scale;

                    if (c.length == 1) {
                        scale = new Affine( c[0], 0, 0, c[0], 0, 0 );
                    } else if (c.length == 2) {
                        scale = new Affine( c[0], 0, 0, c[1], 0, 0 );
                    } else {
                        console.log( "mesh.js: scale does not have 1 or 2 arguments!" );
                        scale = new Affine( 1, 0, 0, 1, 0, 0 );
                    }

                    affine = affine.append(scale);

                    break;

                case "rotate":
                    if (c.length == 3 ) {
                        trans = new Affine( 1, 0, 0, 1, c[1], c[2]);
                        affine = affine.append(trans);
                    }

                    if (c[0]) {
                        radian = c[0] * Math.PI/180.0;
                        cos = Math.cos(radian);
                        sin = Math.sin(radian);

                        if (Math.abs(cos) < 1e-16) { // I hate rounding errors...
                            cos = 0;
                        }
                        if (Math.abs(sin) < 1e-16) { // I hate rounding errors...
                            sin = 0;
                        }
                        var rotate = new Affine(cos, sin, -sin, cos, 0, 0);
                        affine = affine.append(rotate);
                    } else {
                        console.log( "math.js: No argument to rotate transform!" );
                    }

                    if (c.length == 3) {
                        trans = new Affine(1, 0, 0, 1, -c[1], -c[2]);
                        affine = affine.append(trans);
                    }

                    break;

                case "skewX":
                    if (c[0]) {
                        radian = c[0] * Math.PI/180.0;
                        tan = Math.tan(radian);
                        skewx = new Affine( 1, 0, tan, 1, 0, 0 );
                        affine = affine.append(skewx);
                    } else {
                        console.log("math.js: No argument to skewX transform!");
                    }

                    break;

                case "skewY":
                    if (c[0]) {
                        radian = c[0] * Math.PI/180.0;
                        tan = Math.tan(radian);
                        skewy = new Affine( 1, tan, 0, 1, 0, 0 );
                        affine = affine.append(skewy);
                    } else {
                        console.log("math.js: No argument to skewY transform!");
                    }

                    break;

                case "matrix":
                    if (c.length == 6) {
                        var matrix = new Affine( c[0], c[1], c[2], c[3], c[4], c[5] );
                        affine = affine.append(matrix);
                    } else {
                        console.log("math.js: Incorrect number of arguments for matrix!");
                    }

                    break;

                default:
                    console.log("mesh.js: Unhandled transform type: " + type);

                    break;
            }
        }

        return affine;
    }


    //------------------------------------- Mesh Object ------------------------


     //Object to simplify rendering by caching some stuff
    function Mesh(id) {

        this.id = id;
		var themesh = document.getElementById(id);

		this.xsize = 0;
		this.ysize = 0;
		this.points = null;
		this.colors = null;


		 //---------state used by the renderer:
        var Cx, Cy;

        this.V  = [0,0,0,0]; //temp space
        this.TT = [0,0,0,0];
        this.SS = [0,0,0,0]; //temp space for getPoint(double,double)

         //used primarily for PatchData::WhatColor() lookup
        this.s0 = -1;
        this.ds =  0;
        this.t0 = -1;
        this.dt =  0;

         //render context
        var buffer; //a temp, non-local buffer
        var bufferwidth = 0, bufferheight = 0;
        var numchannels = 4;    //usually 4 (argb) and 8bits
        var bitsperchannel = 8; //8 or 16
        var stride = 0; //usually bufferwidth * numchannels * bitsperchannel/8

        var recurse_depth = 0;

        var colUL, colUR, colLL, colLR;


        /*! Return the point (S Cx T,S Cy T).
         *  assumes Cx,Cy already set right.
         * 
         * Called from rpatchpoint().
         */
        flatpoint getPointST(double *S,double *T)
        {
            flatpoint p;
            m_times_v(Cx,T,V); 
            p.x=dot(S,V);
            m_times_v(Cy,T,V);
            p.y=dot(S,V);

            return p;
        }

        /*! Update SS and TT. s and t must be in range [0..1].
         * Recomputes SS and TT, when s!=SS[2] or t!=TT[2]. (see getT())
         */
        flatpoint getPoint(double s,double t)
        {
            if (s != SS[2]) getT(SS,s);
            if (t != TT[2]) getT(TT,t);

            return getPointST();
        }

		/*! Grab either x or y coordinates from a particular mesh square at roffset and coffset.
		 *
		 * This Gt refers only to the one 4x4 coordinate section starting at (roffset,coffset).
		 *
		 * roffset and coffset are point indices (in range [0..ysize) and [0..xsize) respectively),
		 * not subpatch indices.
		 *
		 */
		function getGt(roffset, coffset, isfory) 
		{
			var Gt = [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0];

			if (isfory) {
				for (var r=0; r<4; r++)
					for (var c=0; c<4; c++) 
						Gt[r*4+c] = points[(c+roffset)*xsize+(r+coffset)].y;
			} else {
				for (var r=0; r<4; r++)
					for (var c=0; c<4; c++) 
						Gt[r*4+c] = points[(c+roffset)*xsize+(r+coffset)].x;
			}

			return Gt;
		}

        var MESH_Full_Bezier = "Full_Bezier";
        var MESH_Coons       = "Coons"      ;
        var MESH_Border_Only = "Border_Only";
        var MESH_Linear      = "Linear"     ;

        /*! Interpolate control points according to whichcontrols.
         *  This will not change the controls that exist for each type,
         * but will change all the others as best as it knows how.
         */
        function InterpolateControls(whichcontrols)
        {
            if (whichcontrols == MESH_Full_Bezier || xsize==0 || ysize==0) return;
            
            if (whichcontrols == MESH_Linear) {
                 //redo all but the outermost corners
                var p00,p30,p03,p33;
                p00 = points[0];
                p03 = points[xsize-1];
                p30 = points[(ysize-1)*xsize];
                p33 = points[(ysize-1)*xsize+xsize-1];

                var i, s,t;
                for (var r=0; r<ysize; r++) {
                    for (var c=0; c<xsize; c++) {
                         //if point is a corner point, then skip
                        if ((c==0 && r==0)
                              || (c==0 && r==ysize-1)
                              || (c==xsize-1 && r==0 )
                              || (c==xsize-1 && r==ysize-1))
                            continue;

                        i = r*xsize+c;
                        s = c/(xsize-1);
                        t = r/(ysize-1);
                        //points[i] = s*(t*p33+(1-t)*p03) + (1-s)*(t*p30+(1-t)*p00);
                        points[i] = p33.scale(s*t).add(p03.scale(s*(1-t))).add(p30.scale((1-s)*t).add(p00.scale((1-s)*(1-t)));
                    }
                }

            } else if (whichcontrols == MESH_Border_Only) {
                 //redo all interior points
                var i, s,t;
                var pt,pb,pl,pr;
                for (var r=1; r<ysize-1; r++) {
                    for (var c=1; c<xsize-1; c++) {
                        i = r*xsize+c;
                         
                         //point is weight average of the following 4 points
                        pt = points[c];
                        pb = points[c+(ysize-1)*xsize];
                        pl = points[r*xsize];
                        pr = points[xsize-1+r*xsize];

                        s = c/(xsize-1);
                        t = r/(ysize-1);
                        //points[i] = ((1-s)*pl + s*pr + (1-t)*pt + t*pb)/2;
                        points[i] = pl.scale(.5*(1-s)).add(pr.scale(.5*s).add(pt.scale(.5*(1-t)).add(pb.scale(.5*t));
                    }
                }

            } else if (whichcontrols == MESH_Coons) {
                 //redo all 4 interior points per subpatch:
                 //  p11=1./9*(-4*p00+6*(p01+p10)-2*(p03+p30)+3*(p31+p13)-p33)
                 //  p12=1./9*(-4*p03+6*(p02+p13)-2*(p00+p33)+3*(p32+p10)-p30)
                 //  p21=1./9*(-4*p30+6*(p31+p20)-2*(p33+p00)+3*(p01+p23)-p03)
                 //  p22=1./9*(-4*p33+6*(p32+p23)-2*(p30+p03)+3*(p02+p20)-p00)
                var i;
                var p00,p30,p03,p33, p11,p12,p21,p22, p01,p02,p31,p32,p10,p13,p20,p23;
                for (var r=0; r<ysize-1; r+=3) {
                    for (var c=0; c<xsize-1; c+=3) {
                        i = r*xsize+c; //points p00 corner of subpatch
                        
                        p00 = points[i];
                        p01 = points[i+1];
                        p02 = points[i+2];
                        p03 = points[i+3];
                        
                        p10 = points[i+    xsize];
                        p13 = points[i+3+  xsize];
                        
                        p20 = points[i+  2*xsize];
                        p23 = points[i+3+2*xsize];
                        
                        p30 = points[i+  3*xsize];
                        p31 = points[i+1+3*xsize];
                        p32 = points[i+2+3*xsize];
                        p33 = points[i+3+3*xsize];

						// p11=(-4*p00+6*(p01+p10)-2*(p03+p30)+3*(p31+p13)-p33)/9; 
						// p12=(-4*p03+6*(p02+p13)-2*(p00+p33)+3*(p32+p10)-p30)/9; 
						// p21=(-4*p30+6*(p31+p20)-2*(p33+p00)+3*(p01+p23)-p03)/9; 
						// p22=(-4*p33+6*(p32+p23)-2*(p30+p03)+3*(p02+p20)-p00)/9; 
                        p11 = (p00.scale(-4).add((p01.add(p10).scale(6)).subtract((p03.add(p30).scale(2)).add((p31.add(p13).scale(3)).subtract(p33)) .scale(1/9); 
                        p12 = (p03.scale(-4).add((p02.add(p13).scale(6)).subtract((p00.add(p33).scale(2)).add((p32.add(p10).scale(3)).subtract(p30)) .scale(1/9); 
                        p21 = (p30.scale(-4).add((p31.add(p20).scale(6)).subtract((p33.add(p00).scale(2)).add((p01.add(p23).scale(3)).subtract(p03)) .scale(1/9); 
                        p22 = (p33.scale(-4).add((p32.add(p23).scale(6)).subtract((p30.add(p03).scale(2)).add((p02.add(p20).scale(3)).subtract(p00)) .scale(1/9); 
                        
                        points[i+  xsize +1] = p11;
                        points[i+  xsize +2] = p12;
                        points[i+2*xsize +1] = p21;
                        points[i+2*xsize +2] = p22;
                    }
                }
            }

            //NeedToUpdateCache(0,-1,0,-1);
        }

        // Weighted average to find Bezier points for linear sides.
        function lerp_third(p0, p1) {
			 //p0 + (p1-p0)/3  ==  p0*2/3 + p1*1/3
			return new Flatpoint(p0.x*2.0/3.0 + p1.x/3.0, p0.y*2.0/3.0 + p1.y/3.0)
        }

		function ReadInPoints() {
			if (!themesh) return;

			 //initial point is an attribute of main meshgradient element
			var x = Number(themesh.getAttribute("x")),
				y = Number(themesh.getAttribute("y"));
			var initialPoint = new Flatpoint(x,y);

			var rows = themesh.children;
			if (!rows || rows.length == 0) {
				 //bad file!
				return;
			}

			this.ysize = 3*rows.length + 1; //note: assumes children are ONLY rows!

			try {
				var i = 0, j = 0;

				colors = [];
				points = [];
				points[0] = initialPoint;

				for (int ri=0; ri<rows.length; ri++) {
					i = ri;

					var patches = rows[ri].children;
					if (ri == 0) {
						this.ysize = 3*rows   .length + 1;
						this.xsize = 3*patches.length + 1; //note: assumes children are ONLY patches!
					}

					colors[ri] = [];


					for (int pi=0; pi<patches.length; pi++) {
						j = pi;
						var ii = ri*3*this.xsize + pi*3;
						var stops = patches[pi].children;

						for (int si=0; si<stops.length; si++) {
							var stop = stops[si]; //note: assumes children are ONLY stops!

                            var edge = si;
                            if (ri !== 0) {
                                ++edge; // There is no top if row isn't first row.
                            }

							var path = stop.getAttribute("path");
							if (path === null) throw "Warning! bad mesh definition: missing path";

							var parts = path.match(/\s*([lLcC])\s*(.*)/);
							if (parts === null) throw "Warning! bad mesh definition: path has to be one of lLcC");
							var command = parts[1];
                            var stop_points = parse_points( parts[2] );

							switch(command) {
							  case 'l': //relative
                                  if (edge === 0) { // Top
                                      points[(3*i)*xsize + 3*j+3] = stop_points[0].add(points[(3*i)*xsize + 3*j]);
                                      points[(3*i)*xsize + 3*j+1] = lerp_third( points[(3*i)*xsize + 3*j], points[(3*i)*xsize + 3*j+3] );
                                      points[(3*i)*xsize + 3*j+2] = lerp_third( points[(3*i)*xsize + 3*j+3], points[(3*i)*xsize + 3*j] );

                                  } else if (edge == 1) { // Right
                                      points[(3*i+3)*xsize + 3*j+3] = stop_points[0].add(points[(3*i)*xsize + 3*j+3]);
                                      points[(3*i+1)*xsize + 3*j+3] = lerp_third( points[(3*i)*xsize + 3*j+3], points[(3*i+3)*xsize + 3*j+3] );
                                      points[(3*i+2)*xsize + 3*j+3] = lerp_third( points[(3*i+3)*xsize + 3*j+3], points[(3*i)*xsize + 3*j+3] );

                                  } else if (edge == 2) { // Bottom
                                      if(j===0) {
                                          points[(3*i+3)*xsize + 3*j+0] = stop_points[0].add(points[(3*i+3)*xsize + 3*j+3]);
                                      }
                                      points[(3*i+3)*xsize + 3*j+1] = lerp_third( points[(3*i+3)*xsize + 3*j], points[(3*i+3)*xsize + 3*j+3] );
                                      points[(3*i+3)*xsize + 3*j+2] = lerp_third( points[(3*i+3)*xsize + 3*j+3], points[(3*i+3)*xsize + 3*j] );

                                  } else { // Left
                                      points[(3*i+1)*xsize + 3*j] = lerp_third( points[(3*i)*xsize + 3*j], points[(3*i+3)*xsize + 3*j] );
                                      points[(3*i+2)*xsize + 3*j] = lerp_third( points[(3*i+3)*xsize + 3*j], points[(3*i)*xsize + 3*j] );
                                  }

                                  break;

							  case 'L': //absolute
                                  if (edge === 0) { // Top
                                      points[(3*i)*xsize + 3*j+3] = stop_nodes[0];
                                      points[(3*i)*xsize + 3*j+1] = lerp_third( points[(3*i)*xsize + 3*j], points[(3*i)*xsize + 3*j+3] );
                                      points[(3*i)*xsize + 3*j+2] = lerp_third( points[(3*i)*xsize + 3*j+3], points[(3*i)*xsize + 3*j] );

                                  } else if (edge === 1) { // Right
                                      points[(3*i+3)*xsize + 3*j+3] = stop_nodes[0];
                                      points[(3*i+1)*xsize + 3*j+3] = lerp_third( points[(3*i)*xsize + 3*j+3], points[(3*i+3)*xsize + 3*j+3] );
                                      points[(3*i+2)*xsize + 3*j+3] = lerp_third( points[(3*i+3)*xsize + 3*j+3], points[(3*i)*xsize + 3*j+3] );

                                  } else if (edge === 2) { // Bottom
                                      if(j === 0) {
                                          points[(3*i+3)*xsize + 3*j+0] = stop_nodes[0];
                                      }
                                      points[(3*i+3)*xsize + 3*j+1] = lerp_third( points[(3*i+3)*xsize + 3*j], points[(3*i+3)*xsize + 3*j+3] );
                                      points[(3*i+3)*xsize + 3*j+2] = lerp_third( points[(3*i+3)*xsize + 3*j+3], points[(3*i+3)*xsize + 3*j] );

                                  } else { // Left
                                      points[(3*i+1)*xsize + 3*j] = lerp_third( points[(3*i)*xsize + 3*j], points[(3*i+3)*xsize + 3*j] );
                                      points[(3*i+2)*xsize + 3*j] = lerp_third( points[(3*i+3)*xsize + 3*j], points[(3*i)*xsize + 3*j] );
                                  }

                                  break;

                              case "c":
                                  if (edge === 0) { // Top
                                      points[(3*i)*xsize + 3*j+1] = stop_nodes[0].add(points[(3*i)*xsize + 3*j]);
                                      points[(3*i)*xsize + 3*j+2] = stop_nodes[1].add(points[(3*i)*xsize + 3*j]);
                                      points[(3*i)*xsize + 3*j+3] = stop_nodes[2].add(points[(3*i)*xsize + 3*j]);

                                  } else if (edge === 1) { // Right
                                      points[(3*i+1)*xsize + 3*j+3] = stop_nodes[0].add(points[(3*i)*xsize + 3*j+3]);
                                      points[(3*i+2)*xsize + 3*j+3] = stop_nodes[1].add(points[(3*i)*xsize + 3*j+3]);
                                      points[(3*i+3)*xsize + 3*j+3] = stop_nodes[2].add(points[(3*i)*xsize + 3*j+3]);

                                  } else if (edge === 2) { // Bottom
                                      points[(3*i+3)*xsize + 3*j+2] = stop_nodes[0].add(points[(3*i+3)*xsize + 3*j+3]);
                                      points[(3*i+3)*xsize + 3*j+1] = stop_nodes[1].add(points[(3*i+3)*xsize + 3*j+3]);
                                      if(j === 0) {
                                          points[(3*i+3)*xsize + 3*j+0] = stop_nodes[2].add(points[(3*i+3)*xsize + 3*j+3]);
                                      }
                                  } else { // Left
                                      points[(3*i+2)*xsize + 3*j] = stop_nodes[0].add(points[(3*i+3)*xsize + 3*j]);
                                      points[(3*i+1)*xsize + 3*j] = stop_nodes[1].add(points[(3*i+3)*xsize + 3*j]);
                                  }

                                  break;

                              case "C":
                                  if (edge === 0) { // Top
                                      points[(3*i)*xsize + 3*j+1] = stop_nodes[0];
                                      points[(3*i)*xsize + 3*j+2] = stop_nodes[1];
                                      points[(3*i)*xsize + 3*j+3] = stop_nodes[2];

                                  } else if (edge == 1) { // Right
                                      points[(3*i+1)*xsize + 3*j+3] = stop_nodes[0];
                                      points[(3*i+2)*xsize + 3*j+3] = stop_nodes[1];
                                      points[(3*i+3)*xsize + 3*j+3] = stop_nodes[2];

                                  } else if (edge == 2) { // Bottom
                                      points[(3*i+3)*xsize + 3*j+2] = stop_nodes[0];
                                      points[(3*i+3)*xsize + 3*j+1] = stop_nodes[1];
                                      if (j === 0) {
                                          points[(3*i+3)*xsize + 3*j+0] = stop_nodes[2];
                                      }
                                  } else { // Left
                                      points[(3*i+2)*xsize + 3*j] = stop_nodes[0];
                                      points[(3*i+1)*xsize + 3*j] = stop_nodes[1];
                                  }

                                  break;

							  default:
								throw "Bad path in mesh!";
							} //switch(command)

                            if ((i === 0 && j === 0) || si > 0) {
                                var color_raw = getComputedStyle(stops[si]).stopColor.match(/^rgb\s*\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*\)$/i),
                                    alpha_raw = getComputedStyle(stops[si]).stopOpacity,
                                    alpha = 255;

                                if (alpha_raw) {
                                    alpha = parseInt(alpha_raw * 255);
                                }

                                if (color_raw) {
                                    if (edge === 0) { // upper left corner
                                        colors[i][j] = [];
                                        colors[i][j][0] = parseInt(color_raw[1]);
                                        colors[i][j][1] = parseInt(color_raw[2]);
                                        colors[i][j][2] = parseInt(color_raw[3]);
                                        colors[i][j][3] = alpha; // Alpha

                                    } else if (edge === 1) { // upper right corner
                                        colors[i][j+1] = [];
                                        colors[i][j+1][0] = parseInt(color_raw[1]);
                                        colors[i][j+1][1] = parseInt(color_raw[2]);
                                        colors[i][j+1][2] = parseInt(color_raw[3]);
                                        colors[i][j+1][3] = alpha; // Alpha

                                    } else if (edge === 2) { // lower right corner
                                        colors[i+1][j+1] = [];
                                        colors[i+1][j+1][0] = parseInt(color_raw[1]);
                                        colors[i+1][j+1][1] = parseInt(color_raw[2]);
                                        colors[i+1][j+1][2] = parseInt(color_raw[3]);
                                        colors[i+1][j+1][3] = alpha; // Alpha

                                    } else if (edge === 3) { // lower left corner
                                        colors[i+1][j] = [];
                                        colors[i+1][j][0] = parseInt(color_raw[1]);
                                        colors[i+1][j][1] = parseInt(color_raw[2]);
                                        colors[i+1][j][2] = parseInt(color_raw[3]);
                                        colors[i+1][j][3] = alpha; // Alpha
                                    }
                                }
                            }

						}
					}
				}

			} catch (error) {
				console.log(error);
				return;
			}
		}

		 //now do the actual setup processing
		ReadInPoints();
		InterpolateControls(MESH_Coons); //fill in the center 4 points of a 16 point patch

    }; //end Mesh object


//----------- mesh rendering



     //render the s,t rectangular area s1..s2, t1..t2
    Mesh.prototype.patchpoint = function(context, s1, t1,  s2, t2) {
        d++;
        if (d>20) {
            console.log("recurse at 20, probably an error!");
            return;
        }
        
        flatpoint c1,c2,c3,c4;

        var T = getT(t1);
        var S = getT(s1);
        c1 = context->getPoint(S,T); // computes (S Cx T,S Cy T), is already in screen coords

        //T = getT(t1);
        S = getT(s2);
        c2 = context->getPoint(S,T);

        T = getT(t2);
        S = getT(s1);
        c3 = context->getPoint(S,T);

        //T = getT(t2);
        S = getT(s2);
        c4 = context->getPoint(S,T);

        unsigned long color;
        if ( d>8 ||
                 ((int)c1.x==(int)c2.x && (int)c2.x==(int)c3.x && (int)c3.x==(int)c4.x &&
               (int)c1.y==(int)c2.y && (int)c2.y==(int)c3.y && (int)c3.y==(int)c4.y)) {
             //render the point

            color = pixelfromcolor(coloravg(&col,coloravg(&cola,colUL,colUR,s1), coloravg(&colb,colLL,colLR,s1), t1));

            dp->NewFG(color);
            dp->drawpoint(c1,1,1);//***

        } else {
            //divide into smaller squares:
            // s1,t1         (s1+s2)/2,t1          s2,t1
            // s1,(t1+t2)/2  (s1+s2)/2,(t1+t2)/2   s2,(t1+t2)/2
            // s1,t2         (s1+s2)/2,t2          s2,t2

    //		if (abs((int)c1.x-(int)c4.x)>1 && abs((int)c1.y-(int)c4.y)>1) {
                patchpoint(context, s1,t1, (s1+s2)/2,(t1+t2)/2);
                patchpoint(context, (s1+s2)/2,t1, s2,(t1+t2)/2);
                patchpoint(context, s1,(t1+t2)/2, (s1+s2)/2,t2);
                patchpoint(context, (s1+s2)/2,(t1+t2)/2, s2,t2);
    //		}
        }

        d--;
    }

    /*! Draws one patch.
     * called by drawpatches(). 
     * The whole patch is made of potentially a whole lot of adjacent patches.
     *
     * This function prepares up colUL,colUR,colLL,colLR and Cx,Cy matrices for patchpoint2().
     *
     * roff,coff is which patch, point start is == xoff*3
     */
    Mesh.prototype.drawpatch = function(roff, coff)
    {
        //DBG cerr <<"Draw Color Patch: roff:"<<roff<<"  coff:"<<coff<<"   mode:"<<rendermode<<endl;

        flatpoint fp;

        var C;
        var Gtx = getGt(roff*3, coff*3, 0);
        var Gty = getGt(roff*3, coff*3, 1);

        for (var r=0; r<4; r++) {
            for (var c=0; c<4; c++) {
                fp = flatpoint(Gtx[c*4+r],Gty[c*4+r]);
                fp = dp->realtoscreen(fp);
                Gtx[c*4+r] = fp.x;
                Gty[c*4+r] = fp.y;
            }
        }


        colUL = colors [roff  ][coff  ];
        colUR = colors [roff  ][coff+1];
        colLL = colors [roff+1][coff  ];
        colLR = colors [roff+1][coff+1];
        
        C  = m_times_m(B,Gty);
        Cy = m_times_m(C,B);
        C  = m_times_m(B,Gtx);
        Cx = m_times_m(C,B);  //Cx = B Gtx B

    	patchpoint(0,0,1,1); //draw all points
    }

     //step over each square defined in the mesh
    Mesh.prototype.drawpatches = function()
    {
        for (var r=0; r < ysize/3; r++) {
            for (var c=0; c < xsize/3; c++) {
                drawpatch(r,c);
            }
        }
    }



    //---------------------- Check over all fillables in document

    var fillables = document.querySelectorAll('rect,circle,ellipse,path,text');

    for (var i = 0; i < fillables.length; ++i) {
        var fillable = fillables[i];
		if (!fillable.id) {
			fillable.id = "meshjs" + i;
		}
        console.log( "id: " + fillable.id );


		 //get mesh element from fill property like "url(#myMesh)"
        var fill = fillable.style.fill;
        var url_value = fill.match(/^url\(\s*\"?\s*#([^\s\"]+)\"?\s*\)/);

        if ( !(url_value && url_value[1]) ) continue;

        var mesh = document.getElementById(url_value[1]);
        if (mesh.nodeName !== "meshgradient" ) continue;


        //console.log( "We now have a confirmed mesh" );

         //Make a canvas to hold mesh render
        //var my_canvas = document.createElementNS( xhtmlNS, "canvas" );
        var my_canvas = document.createElement( "canvas" );  // Both work for HTML...
        var bbox = fillable.getBBox();
        my_canvas.width  = bbox.width;
        my_canvas.height = bbox.height;

        // console.log ( "Canvas: " + my_canvas.width + "x" + my_canvas.height );
        var my_context = my_canvas.getContext("2d");

        var my_canvas_image = my_context.getImageData( 0, 0, my_canvas.width, my_canvas.height);
        var my_data = my_canvas_image.data;


         // Create mesh object
        var my_mesh = new Mesh( url_value[1] );

        // Adjust for bounding box if necessary.
        if (mesh.getAttribute( "gradientUnits" ) === "objectBoundingBox") {
            my_mesh.scale( new Point( bbox.width, bbox.height ) );
        }

        // Apply gradient transform.
        var gradientTransform = mesh.getAttribute("gradientTransform");
        // console.log( typeof gradientTransform );
        if ( gradientTransform != null ) {
            var affine = parseTransform( gradientTransform );
            my_mesh.transform( affine );
        }

        // Position to Canvas coordinate.
        var t = new Point( -bbox.x, -bbox.y );
        if (mesh.getAttribute( "gradientUnits" ) === "userSpaceOnUse") {
            my_mesh.transform(t);
        }

        // Do the actual rendering
        my_mesh.paint(my_data, my_canvas.width);


		// all done rendering! put back buffer to the image element
        my_context.putImageData(my_canvas_image, 0, 0);

        // Create image element of correct size
        var my_image = document.createElementNS( svgNS, "image" );
        my_image.setAttribute("width", my_canvas.width);
        my_image.setAttribute("height",my_canvas.height);
        my_image.setAttribute("x", bbox.x);
        my_image.setAttribute("y", bbox.y);

        // Set image to data url
        var my_png = my_canvas.toDataURL();
        my_image.setAttributeNS(xlinkNS, "xlink:href", my_png);

        // Insert image into document
        fillable.parentNode.insertBefore( my_image, fillable );
        fillable.style.fill = "none";

        // Create clip referencing shape and insert into document
        var clip = document.createElementNS( svgNS, "clipPath");
        var clip_id = "patchjs_clip" + i;
        clip.setAttribute("id", clip_id);
        var use = document.createElementNS( svgNS, "use");
        use.setAttributeNS( xlinkNS, "xlink:href", "#" + fillable.id);
        clip.appendChild(use);
        fillable.parentElement.insertBefore(clip, fillable);
        my_image.setAttribute("clip-path", "url(#" + clip_id + ")");
    }

})();

