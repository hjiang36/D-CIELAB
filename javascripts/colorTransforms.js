function colorTransformMatrix(name){
	// Color Transformation matrix
	// Multiply the matrix with the column vector to get output color
	var transM;
	switch (name) {
		case "srgb2xyz":
			transM = [[0.4124564,  0.3575761,  0.1804375], 
				  [0.2126729,  0.7151522,  0.0721750],
 				  [0.0193339,  0.1191920,  0.9503041]];
			break;
		case "xyz2lms":
			transM = [[ 0.2689,    0.8518,   -0.0358],
   				  [-0.3962,    1.1770,    0.1055],
    				  [ 0.0214,   -0.0247,    0.5404]];
			break;
		case "xyz2srgb":
			transM = [[ 3.2404542, -1.5371385, -0.4985314],
				  [-0.9692660,  1.8760108,  0.0415560],
 				  [ 0.0556434, -0.2040259,  1.0572252]];
			break;
		case "lms2xyz":
			transM = [[ 1.7896,   -1.2865,    0.3645],
    				  [ 0.6079,    0.4072,   -0.0379],
   				  [-0.0499,    0.0808,    1.8040]];
			break;
		default:
			transM = null;
	}
	return transM;
}

function applyColorTransform(c, transM) {
	// Transform color c with matrix transM
	if (c.length != transM[0].length) {
		return null;
	}
	var cOut = new Array(transM.length);
	for (row = 0; row < transM.length; row++) {
		cOut[row] = 0;
		for (col = 0; col < c.length; col++) {
			cOut[row] = cOut[row] + c[col] * transM[row][col];
		}
	}
	return cOut;
}

function srgb2xyz(srgb) {
	if (srgb.length != 3) {
		return null;
	}
	var transM = colorTransformMatrix("srgb2xyz");
	var XYZ = applyColorTransform(srgb, transM);
	return XYZ;
}

function xyz2srgb(xyz) {
	if (xyz.length != 3) {
		return null;
	}
	var transM = colorTransformMatrix("xyz2srgb");
	var srgb = applyColorTransform(xyz, transM);
	return srgb;
}

function lms2xyz(lms) {
	if (lms.length != 3) {
		return null;
	}
	var transM = colorTransformMatrix("lms2xyz");
	var XYZ = applyColorTransform(lms, transM);
	return XYZ;
}

function xyz2lms(xyz) {
	if (xyz.length != 3) {
		return null;
	}
	var transM = colorTransformMatrix("xyz2lms");
	var LMS = applyColorTransform(xyz, transM);
	return LMS;
}

function xyy2xyz(xyy) {
	if (xyy.length != 3) {
		return null;
	}
	var XYZ = new Array(3);
	XYZ[1] = xyy[2];
	XYZ[0] = xyy[0] * xyy[2] / xyy[1];
	XYZ[2] = (1 - xyy[0] - xyy[1]) * xyy[2] / xyy[1];
	return XYZ;	
}

function xyz2xyy(xyz) {
	if (xyz.length != 3) {
		return null;
	}
	var xyY = new Array(3);
	xyY[2] = xyz[1];
	xyY[0] = xyz[0] / (xyz[0] + xyz[1] + xyz[2]);
	xyY[1] = xyz[1] / (xyz[0] + xyz[1] + xyz[2]);
	return xyY;
}

function brettelTransform(c, cbType){
	// Compute brettel equivalent color for LMS color c
	// The output will also be in LMS space
	if (c.length != 3) {
		return null;
	}
	var cOut = c;
	var e = [79.1286, 64.6775, 39.067];
	var anchor = [0.0801, 0.1579, 0.5897,
		      0.1284, 0.2237, 0.3636,
		      0.9856, 0.7325, 0.0011,
		      0.0914, 0.0070, 0];
	var coef = [[e[1]*anchor[8] - e[2]*anchor[7], e[2]*anchor[6]-e[0]*anchor[8], e[0]*anchor[7] - e[1]*anchor[6]],
		    [e[1]*anchor[2] - e[2]*anchor[1], e[2]*anchor[0]-e[0]*anchor[2], e[0]*anchor[1] - e[1]*anchor[0]]];
	switch (cbType){
		case 1: // Protanopia
			if (c[2] / c[1] < e[2] / e[1]) {
				cOut[0] = - (coef[0][1]*c[1] + coef[0][2]*c[2])/coef[0][0];
			}
			else {
				cOut[0] = - (coef[1][1]*c[1] + coef[1][2]*c[2])/coef[1][0];
			}
			break;
		case 2: // Deutanopia
			if (c[2] / c[0] < e[2] / e[0]) {
				cOut[1] = - (coef[0][0] * c[0] + coef[0][2]*c[2])/coef[0][1];
			}
			else {
				cOut[1] = - (coef[1][0] * c[0] + coef[1][2]*c[2])/coef[1][1];
			}
			break;
		case 3: // Tritanopia
			coef = [[e[1]*anchor[11] - e[2]*anchor[10], e[2]*anchor[9]-e[0]*anchor[11], e[0]*anchor[10]-e[1]*anchor[9]],
				[e[1]*anchor[5] - e[2]*anchor[4], e[2]*anchor[3] - e[0]*anchor[5], e[0]*anchor[4]-e[1]*anchor[3]]];
			if (c[1] / c[0] < e[1] / e[0]) {
				cOut[2] = - (coef[0][0] * c[0] + coef[0][1]*c[1])/coef[0][2];
			}
			else {
				cOut[2] = - (coef[1][0] * c[0] + coef[1][1]*c[1])/coef[1][2];
			}
			break;
		default:
			break;

	}
	return cOut;
}

function xyz2lab(xyz, whiteXYZ) {
	// xyz to lab
	// Reference: 
	//   http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_Lab.html
	if (xyz.length != 3 || whiteXYZ.length != 3)
		return null;
	var Lab = new Array(3);
	var f = new Array(3);
	var e = 0.008856;
	var k = 903.3;
	for (ii = 0; ii < 3; ii++){
		xyz[ii] = xyz[ii] / whiteXYZ[ii];
		if (xyz[ii] > 1) return null;
		if (xyz[ii] > e)
			f[ii] = Math.pow(xyz[ii], 1/3);
		else
			f[ii] = (k*xyz[ii] + 16) / 116;
	}
	Lab[0] = 116 * f[1] - 16;
	Lab[1] = 500 * (f[0] - f[1]);
	Lab[2] = 200 * (f[1] - f[2]);
	return Lab
}

function deltaEab(xyz1, xyz2, whiteXYZ) {
	// deltaEab 1976
	var dLab = new Array(3);
	var Lab1 = xyz2lab(xyz1, whiteXYZ);
	var Lab2 = xyz2lab(xyz2, whiteXYZ);
	for (ii = 0; ii < 3; ii++)
		dLab[ii] = Lab1[ii] - Lab2[ii];
	return Math.sqrt(dLab[0]*dLab[0] + dLab[1]*dLab[1] + dLab[2]*dLab[2]);
}
