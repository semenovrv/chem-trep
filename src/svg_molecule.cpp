#include <string>

#include <cmath>
#include "svg_molecule.h"
#include "simple_molecule.h"
#include "molutils.h"

/*
QFont myFont(fontName, fontSize);;
QString str("I wonder how wide this is?");

QFontMetrics fm(myFont);
int width=fm.horizontalAdvance(str);
*/

namespace MolStruct {

	void crossLines(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double & x, double & y) {
		double a, b, c, d, det;
	
		x = 1000000;
		y = 1000000;
		a = y2 - y1; 
		b = x2 - x1; 
		c = y4 - y3; 
		d = x4 - x3;
		det = c*b - a*d;
		if (abs(det) < 0.00001) return;
		y = (a*c*(x3 - x1) + y1*c*b - y3*a*d) / det;
		x = -(d*b*(y3 - y1) + x1*d*a - x3*c*b) / det;
	};

	void drawArrowSVG(std::string & pathList, std::vector<std::string> & polygonList, double x1, double y1, double x2, double y2, double lgth, double widthDIV2) {

		//	x1, y1 - from
		//		x2, y2 - to
		//		width - arrow width(pixels);
		//	lgth - length of nakonechnik(opixel)
		//		widthDIV2 - half of width of nakonechnik

		double k1, a, b, c, d;
		double xp, yp, x3, y3, x4, y4, xx;
		std::string ss;
		
			//if Width>=0 then Cv.Pen.Width := round(width);
			//  Cv.Pen.Color := clBlack;
			//  Cv.Brush.Color := clBlack;
			//Cv.PenPos := Point(round(x1), round(y1));
			//Cv.LineTo(round(x2), round(y2));
		ss = " M " + formatPrecision(x1,1) + " " + formatPrecision(y1,1) + " L " + formatPrecision(x2,1) + " " + formatPrecision(y2,1);
		pathList = pathList + ss;
		if (abs(x2 - x1) > 8) {
			k1 = (y2 - y1) / (x2 - x1);
			a = k1*k1 + 1;
			b = 2 * k1*(x2*y1 - x1*y2) / (x2 - x1) - 2 * y2*k1 - 2 * x2;
			xx = (x2*y1 - x1*y2) / (x2 - x1);
			c = x2*x2 + y2*y2 - lgth*lgth + xx*xx - 2 * y2*(x2*y1 - x1*y2) / (x2 - x1);
			d = b*b - 4 * a*c;
			xp = round((-b + sqrt(d)) / (2 * a));
			if (((xp > x1) && (xp > x2)) || ((xp < x1) && (xp < x2))) xp = round((-b - sqrt(d)) / (2 * a));
			yp = round(xp*k1 + (x2*y1 - x1*y2) / (x2 - x1));
			if (y2 != y1) {
				
				x3 = round(xp + widthDIV2*sin(atan(k1)));
				y3 = round(yp - widthDIV2*cos(atan(k1)));
				x4 = round(xp - widthDIV2*sin(atan(k1)));
				y4 = round(yp + widthDIV2*cos(atan(k1)));
				
			} else {
				x3 = xp;
				y3 = yp - widthDIV2;
				x4 = xp;
				y4 = yp + widthDIV2;
			};
		}
		else {
			xp = x2;
			yp = y2 - lgth;
			if (((yp > y1) && (yp > y2)) || ((yp < y1) && (yp < y2))) yp = y2 + lgth;
			x3 = xp - widthDIV2;
			y3 = yp;
			x4 = xp + widthDIV2;
			y4 = yp;
		};
		ss = formatPrecision(x2,1) + "," + formatPrecision(y2,1) + "  " + formatPrecision(x3,1) + "," + formatPrecision(y3,1) + "  " + formatPrecision(x4,1) + "," + formatPrecision(y4,1);
		ss = format("<polygon fill=\"black\" stroke=\"black\" stroke-width=\"1\" points=\"{}\" />", ss);
		polygonList.push_back(ss);
	};

	void svgPolymerSave(const TSimpleMolecule & sm, std::vector<std::string> & svgData, bool numerationOutput, int imgWidth, int imgHeight) {
		MolStruct::TSVGMolecule test;
		std::vector<std::vector<int>*> ringList;
		int nTotalCycles, nAromFive, nAromSixs, nCondensed;
		int i,n;
		double xU, yU, d;

		test.moleculeCopy(sm);
		test.removeExplicitHydrogens();
		test.defineAtomConn();
		test.options.fIOPT1 = 1;
		d = test.averageBondLength();
		if (d == 0) d = 30;
		for (i = test.nAtoms() - 1; i>=0; i--) {
			n = test.getAtom(i)->iz / 3;
			while (n > 0) {
				test.unitVector(i, xU, yU);
				xU = test.getAtom(i)->rx + xU*d;
				yU = test.getAtom(i)->ry + yU*d;
				test.addAtom(ID_ZVEZDA, 0, xU, yU);
				test.addBond(1, i, test.nAtoms() - 1);
				test.defineAtomConn();
				n--;
			}
		}

		test.allAboutCycles();
		test.cyclesCalculate(nTotalCycles, nAromFive, nAromSixs, nCondensed, &ringList);
		test.getSVG(imgWidth, imgHeight, ringList, svgData, true, numerationOutput);

		for (i = 0; i < ringList.size(); i++) delete(ringList[i]);
	}



	void TSVGMolecule::dotLineSVG(double x1, double y1, double x2, double y2, std::vector<std::string> & outBuffer) const {
		double	siX, siY, cs, si, tg, r1, r2, dx;
		bool test;
		std::string s, sX, sY, sR;
		std::vector<std::string> params;
		
		dx = (svgPixPerInch*fSVGLineWidth) + 1;
		r1 = abs(x2 - x1);
		r2 = abs(y2 - y1);
		if (r1 > 0.0000001) tg = r2 / r1; else tg = 1.0E10; //Calculation of tangens
		tg = tg*tg;
		if (tg > 0.0000001) si = sqrt(1 / (1 + 1 / tg)); else si = 0; //Calculation of sine and cosine
		cs = sqrt(1 / (1 + tg));
		if (x2 >= x1) siX = dx; else siX = -dx;
		if (y2 >= y1) siY = dx; else siY = -dx;
		r1 = x1; 
		r2 = y1;
		sX = formatPrecision(r1,1); //str(r1:0 : 1, sX);
		sY = formatPrecision(r2,1); //str(r2:0 : 1, sY);
		sR = formatPrecision(dx,2); //str(Dx:0 : 2, sR);
		test = true;
		while (test) { // Set pixel and four empy pixels until end of line will be achieved 
			//Canvas.Ellipse(round(R1-Dx),round(R2-Dx),round(R1+Dx),round(R2+Dx));
			params.clear();
			params.push_back(sX);
			params.push_back(sY);
			params.push_back(sR);
			s = format("<circle fill=\"black\" stroke=\"black\" stroke-width=\"1\" cx=\"{}\" cy=\"{}\" r=\"{}\" />", params);
			outBuffer.push_back(s);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 -y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
		};

	};

	void TSVGMolecule::dashedLineSVG(double x1, double y1, double x2, double y2 , std::vector<std::string> & outBuffer) const {
		//Draw dashed line from X1, Y1 to X2, Y2 with specified color
		double 	siX, siY, cs, si, tg, r1, r2, dx, t1, t2, t;
		bool test;
	    int polyX[4], polyY[4];
		std::string s;
		int i;
		std::string sX, sY, sR;
		std::vector<std::string> params;

		
		dx = (svgPixPerInch*fSVGLineWidth / 2) + 1;
		r1 = abs(x2 - x1);
		r2 = abs(y2 - y1);
		if (r1 > 0.0000001) tg = r2 / r1; else tg = 1.0E10; //Calculation of tangens
		tg = tg*tg;
		if (tg > 0.0000001) si = sqrt(1 / (1 + 1 / tg)); else si = 0; //Calculation of sine and cosine
		cs = sqrt(1 / (1 + tg));
		if (x2 >= x1) siX = dx; else siX = -dx;
		if (y2 >= y1) siY = dx; else siY = -dx;
		r1 = x1; r2 = y1;
		t1 = x2 - x1;
		t2 = y2 - y1;
		t = sqrt(t1*t1 + t2*t2);
		if (abs(t) < 1E-6) return;
		t1 = t1 / t; t2 = t2 / t;
		test = true;
		while (test) {
			polyX[0] = round(r1 - dx*t2);
			polyX[1] = round(r1 + dx*t2);
			polyY[0] = round(r2 + dx*t1);
			polyY[1] = round(r2 - dx*t1);

			polyX[2] = round(r1 + 3 * siX*cs + dx*t2);
			polyX[3] = round(r1 + 3 * siX*cs - dx*t2);
			polyY[2] = round(r2 + 3 * siY*si - dx*t1);
			polyY[3] = round(r2 + 3 * siY*si + dx*t1);

			test = (abs(r1 + 3 * siX*cs - x2) > dx) || (abs(r2 + 3 * siY*si - y2) > dx);
			if (test) {
				s = "";
				for (i = 0; i < 3; i++) {
					if (s.length() > 0) s = s + " ";
					s = s + std::to_string(polyX[i]) + "," + std::to_string(polyY[i]);
				};
				s = format("<polygon fill=\"black\" stroke=\"black\" stroke-width=\"1\" points=\"{}\" />", s);
				outBuffer.push_back(s);
				
				sX = formatPrecision(r1,1); 
				sY = formatPrecision(r2,1); 
				sR = formatPrecision(dx,1); 
				params.clear();
				params.push_back(sX);
				params.push_back(sY);
				params.push_back(sR);
				s = format("<circle fill=\"black\" stroke=\"black\" stroke-width=\"1\" cx=\"{}\" cy=\"{}\" r=\"{}\" />",params);
				outBuffer.push_back(s);
				//Canvas.Ellipse(round(R1-Dx),round(R2-Dx),round(R1+Dx),round(R2+Dx));
				sX = formatPrecision(r1 + 3 * siX*cs,1); 
				sY = formatPrecision(r2 + 3 * siY*si,1); 
				params.clear();
				params.push_back(sX);
				params.push_back(sY);
				params.push_back(sR);

				s = format("<circle fill=\"black\" stroke=\"black\" stroke-width=\"1\" cx=\"{}\" cy=\"{}\" r=\"{}\" />", params);
				outBuffer.push_back(s);
			};

			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
			r1 = r1 + siX*cs;
			r2 = r2 + siY*si;
			if (test) test = (abs(r1 - x2) > dx) || (abs(r2 - y2) > dx);
		};

	};

	std::string  TSVGMolecule::solidLineSVG(double x1, double y1, double x2, double y2) const {
		std::string result = "";
		std::string sX1, sY1, sX2, sY2;

		sX1 = formatPrecision(x1,1);
		sY1 = formatPrecision(y1,1);
		sX2 = formatPrecision(x2,1);
		sY2 = formatPrecision(y2,1);
		result = " M " + sX1 + " " + sY1 + " L " + sX2 + " " + sY2;

		return result;
	};

	int TSVGMolecule::getTextWidthLarge(const std::string & data) const {
		int arial16Widths[256] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 4, 3, 5, 7, 7, 12, 9, 2, 4, 4, 5, 8, 4, 4, 4, 4, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 4, 4, 8, 8, 8, 7, 13, 9, 9, 9, 9, 9, 8, 10, 9, 3, 6, 9, 7, 11, 9, 10, 9, 10, 9, 9, 7, 9, 9, 13, 7, 9, 7, 4, 4, 4, 5, 7, 4, 7, 7, 7, 7, 7, 3, 7, 7, 3, 3, 7, 3, 11, 7, 7, 7, 7, 4, 7, 4, 7, 5, 9, 7, 7, 7, 4, 3, 4, 8, 10, 11, 7, 3, 5, 4, 13, 7, 7, 7, 14, 14, 4, 13, 8, 11, 9, 7, 3, 3, 4, 4, 5, 7, 13, 0, 13, 12, 4, 11, 6, 7, 7, 4, 8, 7, 6, 7, 6, 3, 7, 9, 10, 9, 7, 8, 4, 10, 3, 5, 7, 3, 3, 5, 7, 7, 4, 7, 14, 7, 7, 3, 9, 7, 3, 9, 9, 9, 7, 9, 9, 11, 8, 9, 9, 8, 9, 11, 9, 10, 9, 9, 9, 7, 8, 11, 7, 10, 9, 11, 11, 10, 11, 8, 9, 13, 9, 7, 7, 7, 5, 8, 7, 9, 6, 7, 7, 6, 7, 9, 7, 7, 7, 7, 7, 5, 7, 9, 7, 7, 7, 11, 12, 8, 9, 7, 7, 10, 7 };
		int result = 0;
		int i;
		char ch;

		for (i = 0; i < data.length(); i++) {
			ch = data[i];
			result = result + arial16Widths[ch];
		}
		return result;
	};

	int TSVGMolecule::getTextWidthSmall(const std::string & data) const {
		int arial10Widths[256] = { 6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,2,3,3,4,4,7,5,2,3,3,3,5,2,3,2,2,4,4,4,4,4,4,4,4,4,4,2,2,5,5,5,4,8,5,5,6,6,5,5,6,6,2,4,5,4,7,6,6,5,6,6,5,5,6,5,7,5,5,5,2,2,2,3,4,3,4,4,4,4,4,2,5,4,2,2,3,2,6,4,4,4,5,3,4,2,4,3,5,3,3,4,3,3,3,5,6,7,4,2,3,3,8,4,4,4,7,8,3,8,5,7,6,4,2,2,3,3,3,4,8,0,8,7,3,7,4,4,4,2,5,3,4,4,4,3,4,5,6,6,4,5,3,6,2,3,4,2,2,3,5,4,3,4,9,4,4,2,5,4,3,5,5,5,4,5,5,7,5,6,6,5,5,7,6,6,6,5,6,5,5,6,5,6,5,7,8,6,7,5,6,8,6,4,5,4,3,5,4,5,4,4,4,4,5,6,4,4,4,4,4,4,3,7,3,5,4,6,7,5,6,4,4,6,4};
		int result = 0;
		int i;
		char ch;

		for (i = 0; i < data.length(); i++) {
			ch = data[i];
			result = result + arial10Widths[ch];
		}
		return result;
	};

	void TSVGMolecule::svgTextOut(double x, double y, double fontSize, int nH, int charge, int iz, int rl, const std::string sData, std::vector<std::string> & outBuffer) const {
		std::string	s, ss;
		double textWidth, textHeight;
		double xDLeft, xDRight, yD;
		double boxY, textY;
		std::string s1, s11, s2, s3, s4, anchor;
		std::string sMain;
		int smallHeight;

		bool leftDefinition;
		std::vector<std::string> boxStrings, letterStrings, params;
		int i,n;
		std::string textStringSpan;
		double boxWidth;
		int nr;
		int hSuperscriptShift, dataHeight, dataWidth;
		double boxRad, r;

		leftDefinition= nH>0;
		nH = abs(nH);
		if (nH>100) nH = 0;
		hSuperscriptShift = fontSize -round(fontSize*svgFontSmallRatio);  //font size == fontHeight
		smallHeight = round(fontSize*svgFontSmallRatio);

		dataWidth = getTextWidthLarge(sData);
		boxWidth = dataWidth + 2;
		dataHeight = fontSize;
		yD = y;
		xDLeft = x - dataWidth / 2 - 1;        //XD
		xDRight = xDLeft + boxWidth;       //XD

		boxY = yD - fontSize / 2;
		textY = yD + 0.38*fontSize;
		

		s1 = formatPrecision(xDLeft,1);
		if (svgDefaultAtomBox == 0) {
			s11 = formatPrecision(boxY,1); //str(boxY:1 : 1, s11);
			s3 = formatPrecision(boxWidth + 1,1); 
			s4 = formatPrecision(dataHeight,1); 
			params.clear();
			params.push_back(s1);
			params.push_back(s11);
			params.push_back(s3);
			params.push_back(s4);
			ss = format("<rect  x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" />", params);  //[x-textWidth div 2,2+y+textHeight div 2,4+textWidth,textHeight]);
		}
		else {
			r = (boxWidth + 1) / 2;
			if ((dataHeight / 2) > r) r = dataHeight / 2;
			if (sData.length() == 1) s4 = formatPrecision(xDLeft + r - 1,1); else s4 = formatPrecision(xDLeft + r,1);
			s11 = formatPrecision(boxY + dataHeight / 2,1); 
			s3 = formatPrecision(r + svgDeltaR,1);  
			params.clear();
			params.push_back(s4);
			params.push_back(s11);
			params.push_back(s3);

			ss = format("<circle  cx=\"{}\" cy=\"{}\" r=\"{}\" />",params);  //[x-textWidth div 2,2+y+textHeight div 2,4+textWidth,textHeight]);
		};
		boxStrings.push_back(ss);
		//Canvas.FillRect(R);
		s2 = formatPrecision(textY,1);
		params.clear();
		params.push_back(s1);
		params.push_back(s2);
		params.push_back(sData);

		ss = format("<text x=\"{}\" y=\"{}\">{}</text>", params);  //[X-textWidth div 2,Y+(3*textHeight)div 2,s]);
		letterStrings.push_back(ss);

		//radical addition-at up
		if (rl != 0) {    //Ne poluchilos radical
			boxRad = yD;
			s1 = formatPrecision(xDLeft,1); //str(xDLeft:1 : 1, s1);
			s11 = formatPrecision(boxRad,1);  //str(boxRad:1 : 1, s11);
			s3 = formatPrecision(boxWidth + 1,1); //str((boxWidth + 1) : 1 : 1, s3);
			s4 = "14";
		};

		if (iz != 0) { //from left
			s = intToStr(iz);
			n = getTextWidthSmall(s);
			dataHeight = round(fSVGFontHeight*svgFontSmallRatio);//getTextHeightSmall(s);
			boxWidth = n + 2;
			xDLeft = xDLeft - boxWidth - 1;     //XD

			s1=formatPrecision(xDLeft,1);
			if (svgDefaultAtomBox == 0) {

				s11=formatPrecision(boxY,1);
				s3 = formatPrecision(boxWidth + 1, 1);
				s4=std::to_string(dataHeight);
				params.clear();
				params.push_back(s1);
				params.push_back(s11);
				params.push_back(s3);
				params.push_back(s4);

				ss = format("<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" />", params);  //[x-textWidth div 2,2+y+textHeight div 2,4+textWidth,textHeight]);
			}
			else {
				r = (boxWidth + 1) / 2;
				if ((dataHeight / 2) > r) r = dataHeight / 2;
				s4=formatPrecision(xDLeft + r,1);
				s11=formatPrecision(boxY + r,1);
				s3=formatPrecision(r + svgDeltaR,2);
				params.clear();
				params.push_back(s4);
				params.push_back(s11);
				params.push_back(s3);

				ss = format("<circle  cx=\"{}\" cy=\"{}\" r=\"{}\" />",params);  //[x-textWidth div 2,2+y+textHeight div 2,4+textWidth,textHeight]);
			};



			boxStrings.push_back(ss);
			s2=formatPrecision(textY - hSuperscriptShift,1);
			params.clear();
			params.push_back(s1);
			params.push_back(s2);
			params.push_back(intToStr(smallHeight));
			params.push_back(s);

			ss = format(" <text x=\"{}\" y=\"{}\" font-size=\"{}\">{}</text>", params);  //[X-textWidth div 2,Y+(3*textHeight)div 2,s]);
			letterStrings.push_back(ss);
		};
		

		if (charge != 0) { //from left
			if (abs(charge) > 1) s = intToStr(abs(charge)); else s = "";
			if (charge < 0) s = s + '-'; else s = s + '+';

			n = getTextWidthSmall(s);
			dataHeight = round(fSVGFontHeight*svgFontSmallRatio);//getTextHeightSmall(s);
			boxWidth = n + 2;
			s1=formatPrecision(xDRight,1);
			if (svgDefaultAtomBox == 0) {
				s11=formatPrecision(boxY,1);
				s3=formatPrecision(boxWidth + 1,1);
				s4=intToStr(dataHeight);
				params.clear();
				params.push_back(s1);
				params.push_back(s11);
				params.push_back(s3);
				params.push_back(s4);

				ss = format("<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" />", params);
			}
			else {
				r = (boxWidth + 1) / 2;
				if ((dataHeight / 2) > r) r = dataHeight / 2;
				s4=formatPrecision(xDRight + r,1);
				s11=formatPrecision(boxY + r,1);
				s3=formatPrecision(r + svgDeltaR,1);
				params.clear();
				params.push_back(s4);
				params.push_back(s11);
				params.push_back(s3);

				ss = format("<circle  cx=\"{}\" cy=\"{}\" r=\"{}\" />", params);  //[x-textWidth div 2,2+y+textHeight div 2,4+textWidth,textHeight]);
			};

			boxStrings.push_back(ss);
			s2=formatPrecision(textY - hSuperscriptShift,1);
			params.clear();
			params.push_back(s1);
			params.push_back(s2);
			params.push_back(intToStr(smallHeight));
			params.push_back(s);

			ss = format(" <text x=\"{}\" y=\"{}\" font-size=\"{}\">{}</text>", params);
			letterStrings.push_back(ss);
			xDRight = xDRight + boxWidth + 1;     //XD
		};
		
		if (nH != 0) {
			textStringSpan = "";
			boxWidth = 0;
			if (leftDefinition) {
				if (nH > 1) {
					s = intToStr(nH);
					params.clear();
					params.push_back(intToStr(smallHeight));
					params.push_back(s);

					textStringSpan = format("<tspan font-size=\"{}\">{}</tspan>", params) + textStringSpan;
					n = getTextWidthSmall(s);
					boxWidth = boxWidth + n + 1;
				};
				s = "H";

				params.clear();
				params.push_back(intToStr(round(fontSize)));
				params.push_back(s);
				textStringSpan = format("<tspan font-size=\"{}\">{}</tspan>",params) + textStringSpan;
				n = getTextWidthLarge(s);
				boxWidth = boxWidth + n + 2;
				dataHeight = round(fontSize);


				s1=formatPrecision(xDLeft - boxWidth,1);
				if (svgDefaultAtomBox == 0) {
					s11 = formatPrecision(boxY, 1);
					s3=formatPrecision(boxWidth + 1,1);
					s4=intToStr(dataHeight);
					params.clear();
					params.push_back(s1);
					params.push_back(s11);
					params.push_back(s3);
					params.push_back(s4);

					ss = format("<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" />",params);
				}
				else {
					r = (boxWidth + 1) / 2;
					if ((dataHeight / 2) > r)  r = dataHeight / 2;

					s4=formatPrecision(xDLeft - boxWidth + r,1);
					s11=formatPrecision(boxY + dataHeight / 2,1);
					s3=formatPrecision(r + svgDeltaR,2);
					params.clear();
					params.push_back(s4);
					params.push_back(s11);
					params.push_back(s3);

					ss = format("<circle  cx=\"{}\" cy=\"{}\" r=\"{}\" />", params);

				};

				boxStrings.push_back(ss);
				s2 = formatPrecision(textY, 1);
				params.clear();
				params.push_back(s1);
				params.push_back(s2);
				params.push_back(intToStr(smallHeight));
				params.push_back(textStringSpan);

				ss = format(" <text x=\"{}\" y=\"{}\" font-size=\"{}\">{}</text>",params);
				letterStrings.push_back(ss);
				xDLeft = xDLeft - boxWidth - 1;     //XD
			}
			else {
				if (nH > 1) {
					s = intToStr(nH);
					params.clear();
					params.push_back(intToStr(smallHeight));
					params.push_back(s);

					textStringSpan = format("<tspan font-size=\"{}\">{}</tspan>",params) + textStringSpan;
					n = getTextWidthSmall(s);
					boxWidth = boxWidth + n + 1;
				};
				s = "H";
				params.clear();
				params.push_back(intToStr(fontSize));
				params.push_back(s);

				textStringSpan = format("<tspan font-size=\"{}\">{}</tspan>", params) + textStringSpan;
				n = getTextWidthLarge(s);
				boxWidth = boxWidth + n + 2;
				dataHeight = fontSize;


				s1=formatPrecision(xDRight,1);
				if (svgDefaultAtomBox == 0) {
					s11 = formatPrecision(boxY, 1);
					s3 = formatPrecision(boxWidth + 2, 1);
					s4 = formatPrecision(dataHeight,1);
					params.clear();
					params.push_back(s1);
					params.push_back(s11);
					params.push_back(s3);
					params.push_back(s4);

					ss = format("<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" />", params);
				}
				else {
					r = (boxWidth + 1) / 2;
					if ((dataHeight / 2) > r)  r = dataHeight / 2;
					s4 = formatPrecision(xDRight + r, 1);
					s11 = formatPrecision(boxY + r, 1);
					s3=formatPrecision(r + svgDeltaR,1);
					params.clear();
					params.push_back(s4);
					params.push_back(s11);
					params.push_back(s3);

					ss = format("<circle  cx=\"{}\" cy=\"{}\" r=\"{}\" />", params);

				};
				boxStrings.push_back(ss);
				s2=formatPrecision(textY,1);
				params.clear();
				params.push_back(s1);
				params.push_back(s2);
				params.push_back(intToStr(smallHeight));
				params.push_back(textStringSpan);

				ss = format(" <text x=\"{}\" y=\"{}\" font-size=\"{}\">{}</text>", params);
				letterStrings.push_back(ss);
				xDRight = xDRight + boxWidth + 1;     //XD
			};
		};
		

		for (i = 0; i < boxStrings.size(); i++) outBuffer.push_back(boxStrings[i]);
		for (i = 0; i<letterStrings.size(); i++) outBuffer.push_back(letterStrings[i]);
		
	};

	double TSVGMolecule::cosB(int b1, int b2) const {
		double r1, r2, x1, x2, y1, y2;
		double result = 0;
		int j1, j2;
		
		j1 = getBond(b1)->at[0];
		j2 = getBond(b1)->at[1];
		x1 = getAtom(j1)->rx - getAtom(j2)->rx;
		y1 = getAtom(j1)->ry - getAtom(j2)->ry;
		j1 = getBond(b2)->at[0];
		j2 = getBond(b2)->at[1];
		x2 = getAtom(j1)->rx - getAtom(j2)->rx;
		y2 = getAtom(j1)->ry - getAtom(j2)->ry;

		r1 = sqrt(x1*x1 + y1*y1);
		r2 = sqrt(x2*x2 + y2*y2);
		if ((r1 != 0) && (r2 != 0)) result = (x1*x2 + y1*y2) / (r1*r2);
		return result;
	};


	std::string TSVGMolecule::bDrawerSVG(std::vector<std::string> & polygonData, int bondN, double ml, const std::vector<singleCenter> & center) const {
		std::string result = "";
		int i, si;
		double x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, dx, dy, xN1, xN2, yN1, yN2;
		int lineType, j1, j2, j3, j4, j5, j6, j7, j8;  
		double cR1, r, rx, ry, lx, rx1, ry1, s1, s2, xc, yc, dist, dist1, r1, r2, cR, rR, xX1, yY1;
		std::string des;
		bool test;
		int polyDX[3];
		int polyDY[3];
		std::string s,ss;
		double rRMin, rRMax;

		j1 = getBond(bondN)->at[0]; j2 = getBond(bondN)->at[1];
		x1 = getAtom(j1)->rx; y1 = getAtom(j1)->ry;
		x2 = getAtom(j2)->rx; y2 = getAtom(j2)->ry;
		lineType = getBond(bondN)->tb;
		
		rR = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
		cR = 0;
		if (center[bondN].cSize > 0) cR = svgCenterRatio; else {
			if (getBond(bondN)->tb == 3) cR = svgCenterRatio2; else cR = svgCenterRatio1;
		};
		dist= abs(cR);
		if (cR > 0) { //double bond 
			dist = ml*dist / 100; dist = dist / 2;
		}
		else dist = dist*svgPixPerInch;
		rRMin = fMinLineDistance*svgPixPerInch;
		rRMax = fMaxLineDistance*svgPixPerInch;
		if (dist<rRMin / 2) dist = rRMin / 2;
		if (dist>rRMax / 2) dist = rRMax / 2;
		switch (lineType) {
			case 1:
			case 6: 
			case 8:
				
				if (center[bondN].cSize == 0) {
					j3 = getAtom(j1)->nb;
					if (j3 == 2) {
						j4 = 0;
						test = false;
						while (!test) {
							if (j4 != bondN) test = (getBond(j4)->at[0] == j1) || (getBond(j4)->at[1] == j1);
							if (!test) {
								j4++;
								if (j4 == nBonds()) test = true;
							}
						}; // until test or (j4 = fBond.nBonds);
						if ((getBond(j4)->tb >= 2) && (getBond(j4)->tb <= 5)) {
							r = cosB(bondN, j4);
							r = sqrt(abs(1 - r*r));
							//!!!!
							if (center[j4].cSize > 0)  cR1 = svgCenterRatio; else {
								if (getBond(j4)->tb == 3) cR1 = svgCenterRatio2; else cR1 = svgCenterRatio1;
							};
							dist1 = abs(cR1);
							if (cR1 > 0) {  // double bond
								dist1 = ml*dist1 / 100; dist1 = dist1 / 2;
							}
							else dist1 = dist1*svgPixPerInch;
							if (dist1 < rRMin / 2)  dist1 = rRMin / 2;
							if (dist1 > rRMax / 2)  dist1 = rRMax / 2;

							//!!!!
							if (r > 0.3) { // correction required 
								r1 = dist1 / r;
								dx = x1 - x2; 
								dy = y1 - y2;
								r2 = sqrt(dx*dx + dy*dy);
								dx = (r1 + r2)*dx / r2;
								dy = (r1 + r2)*dy / r2;
								x1 = x2 + dx; 
								y1 = y2 + dy;
							};
						};
					};
					j3 = getAtom(j2)->nb;
					if (j3 == 2) {
						j4 = 0;
						test = false;
						while (!test) {
						
							if (j4 != bondN) test = (getBond(j4)->at[0] == j2) || (getBond(j4)->at[1] == j2);
							if (!test) {
								j4++;
								if (j4 == nBonds()) test = true;
							}
						}; // until test or (j4 = fBond.nBonds);
						if ((getBond(j4)->tb >= 2) && (getBond(j4)->tb <= 5)) {

							cR1 = 0;
							if (center[j4].cSize > 0) cR1 = svgCenterRatio; else {
								if (getBond(j4)->tb == 3)  cR1 = svgCenterRatio2; else cR = svgCenterRatio1;
							};
							dist1 = abs(cR1);
							if (cR1 > 0) { //double bond 
								dist1 = ml*dist1 / 100; 
								dist1 = dist1 / 2;
							}
							else dist1 = dist1*svgPixPerInch;
							//!!!!
							r = cosB(bondN, j4);
							r = sqrt(abs(1 - r*r));
							if (r > 0.3) { //correction required 
								r1 = dist1 / r;
								dx = x2 - x1; 
								dy = y2 - y1;
								r2 = sqrt(dx*dx + dy*dy);
								dx = (r1 + r2)*dx / r2;
								dy = (r1 + r2)*dy / r2;
								x2 = x1 + dx; 
								y2 = y1 + dy;
							};
						};
					};
				};
				if (lineType == 1) result = result + solidLineSVG(x1, y1, x2, y2); else
				if (lineType == 6) dotLineSVG(x1, y1, x2, y2, polygonData); else
				if (lineType == 8) dashedLineSVG(x1, y1, x2, y2, polygonData);
			
			
			break;
		case 2: 
		case 3: 
		case 4: 
		case 5:
		case 7:
			
			if ((center[bondN].cSize > 0) && (lineType != 3)) {
				r = abs(cR);
				if (cR > 0)  r = 1 - r / 100; else r = 2 * r*svgPixPerInch;
				if (lineType == 5) dashedLineSVG(x1, y1, x2, y2, polygonData); else result = result + solidLineSVG(x1, y1, x2, y2);
				x3 = x1 - center[bondN].centerX;
				y3 = y1 - center[bondN].centerY;
				x4 = x2 - center[bondN].centerX;
				y4 = y2 - center[bondN].centerY;
				if (cR < 0) {
					x5 = x4 - x3;
					y5 = y4 - y3;
					x6 = x4*y3 - x3*y4;
					y6 = x5*x5 + y5*y5;
					y6 = sqrt((x6*x6) / y6);
					r = 1 - r / y6;
				};
				//dist checking..
				rR = sqrt(x3*x3 + y3*y3)*(1 - r);
				if (rR < rRMin) r = 1 - rRMin / sqrt(x3*x3 + y3*y3);
				if (rR > rRMax)  r = 1 - rRMax / sqrt(x3*x3 + y3*y3);
				if (r < 0.5) r = 0.5;
				if (r > 10) r = 0.9;
				x3 = center[bondN].centerX + r*x3;
				y3 = center[bondN].centerY + r*y3;
				x4 = center[bondN].centerX + r*x4;
				y4 = center[bondN].centerY + r*y4;
				if ((lineType == 4) || (lineType == 5)) dashedLineSVG(x3, y3, x4, y4, polygonData); else result = result + solidLineSVG(x3, y3, x4, y4);
			}
			else {
				x3 = x1; y3 = y1; x4 = x2; y4 = y2;
				x5 = x1; y5 = y1; x6 = x2; y6 = y2;
				dx = x2 - x1; dy = y2 - y1;
				if (getBond(bondN)->tb == 14) { dx = abs(dx); dy = abs(dy); };
				r = sqrt(dx*dx + dy*dy);
				if (r > 0) {
					dx = dist*dx / r; 
					dy = dist*dy / r;
				};
				x3 = x3 + dy; y3 = y3 - dx; x4 = x4 + dy; y4 = y4 - dx;
				x5 = x5 - dy; y5 = y5 + dx; x6 = x6 - dy; y6 = y6 + dx;
				j3 = getAtom(j1)->nb;
				if (j3 == 2) {
					j4 = 0;
					test = false;
					while (!test) {
						if (j4 != bondN) test = (getBond(j4)->at[0] == j1) || (getBond(j4)->at[1] == j1); 
						if (!test) {
							j4++;
							if (j4 == nBonds()) test = true;
						};
					}; 
					r = cosB(bondN, j4);
					r = sqrt(abs(1 - r*r));
					if (r > 0.3) {
						j5 = getBond(j4)->at[0];
						xN1 = getAtom(j5)->rx; yN1 = getAtom(j5)->ry;
						j5 = getBond(j4)->at[1];
						xN2 = getAtom(j5)->rx; yN2 = getAtom(j5)->ry;
						crossLines(x3, y3, x4, y4, xN1, yN1, xN2, yN2, rx, ry);
						x3 = rx; y3 = ry;
						crossLines(x5, y5, x6, y6, xN1, yN1, xN2, yN2, rx, ry);
						x5 = rx; y5 = ry;
						crossLines(x1, y1, x2, y2, xN1, yN1, xN2, yN2, rx, ry);
						x1 = rx; y1 = ry;
					};
				}
				else if (j3 == 3) {
					j4 = 0;
					j5 = -1;
					j6 = -1;
					test = false;
					while ((j6<0) && (j4<nBonds())) {
						//j4: = j4 + 1;
						if (j4 != bondN) test = (getBond(j4)->at[0] == j1) || (getBond(j4)->at[1] == j1);
						if (test) {
							if (j5 < 0) j5 = j4; else j6 = j4;
						};
					   test = false;
					   j4++;
					}	//until(j6 < >0) or (j4 = fBond.nBonds);
					r = cosB(bondN, j5);
					r = sqrt(abs(1 - r*r));
					if (r > 0.3) {
						j7 = getBond(j5)->at[0];
						xN1 = getAtom(j7)->rx; yN1 = getAtom(j7)->ry;
						j7 = getBond(j5)->at[1];
						xN2 = getAtom(j7)->rx; yN2 = getAtom(j7)->ry;
						crossLines(x3, y3, x4, y4, xN1, yN1, xN2, yN2, rx, ry);
						crossLines(x5, y5, x6, y6, xN1, yN1, xN2, yN2, rx1, ry1);
						if (xN1 < xN2) { lx = xN1; xN1 = xN2; xN2 = lx; };
						if (yN1 < yN2) { lx = yN1; yN1 = yN2; yN2 = lx; };
						if ((rx >= (xN2 - 0.01)) && (rx <= (xN1 + 0.01)) && (ry >= (yN2 - 0.01)) && (ry <= (yN1 + 0.01))) {
							x3 = rx; y3 = ry;
						}
						else if ((rx1 >= (xN2 - 0.01)) && (rx1 <= (xN1 + 0.01)) && (ry1 >= (yN2 - 0.01)) && (ry1 <= (yN1 + 0.01))) {
							x5 = rx1; y5 = ry1;
						};//	else return;;
					};
					r = cosB(bondN, j6);
					r = sqrt(abs(1 - r*r));
					if (r > 0.3) {
						j7 = getBond(j6)->at[0];
						xN1 = getAtom(j7)->rx; yN1 = getAtom(j7)->ry;
						j7 = getBond(j6)->at[1];
						xN2 = getAtom(j7)->rx; yN2 = getAtom(j7)->ry;
						crossLines(x3, y3, x4, y4, xN1, yN1, xN2, yN2, rx, ry);
						crossLines(x5, y5, x6, y6, xN1, yN1, xN2, yN2, rx1, ry1);
						if (xN1 < xN2) { lx = xN1; xN1 = xN2; xN2 = lx; };
						if (yN1 < yN2) { lx = yN1; yN1 = yN2; yN2 = lx; };
						if ((rx >= (xN2 - 0.01)) && (rx <= (xN1 + 0.01)) && (ry >= (yN2 - 0.01)) && (ry <= (yN1 + 0.01))) {
							x3 = rx; y3 = ry;
						}
						else if ((rx1 >= (xN2 - 0.01)) && (rx1 <= (xN1 + 0.01)) && (ry1 >= (yN2 - 0.01)) && (ry1 <= (yN1 + 0.01))) {
							x5 = rx1; y5 = ry1;
						};	//else return;
					};
				};
				j3 = getAtom(j2)->nb;
				if (j3 == 2) {
					j4 = 0;
					test = false;
					while (!test) {					
						if (j4 != bondN) test = (getBond(j4)->at[0] == j2) || (getBond(j4)->at[1] == j2);
						if (!test) {
							j4++;
							if (j4 == nBonds()) test = true;
						};

					}; //until test or (j4 = fBond.nBonds);
					r = cosB(bondN, j4);
					r = sqrt(abs(1 - r*r));
					if (r > 0.3) {
						j5 = getBond(j4)->at[0];
						xN1 = getAtom(j5)->rx; yN1 = getAtom(j5)->ry;
						j5 = getBond(j4)->at[1];
						xN2 = getAtom(j5)->rx; yN2 = getAtom(j5)->ry;
						crossLines(x3, y3, x4, y4, xN1, yN1, xN2, yN2, rx, ry);
						x4 = rx; y4 = ry;
						crossLines(x5, y5, x6, y6, xN1, yN1, xN2, yN2, rx, ry);
						x6 = rx; y6 = ry;
						crossLines(x1, y1, x2, y2, xN1, yN1, xN2, yN2, rx, ry);
						x2 = rx; y2 = ry;
					};
				}
				else if (j3 == 3) {
					j4 = 0;
					j5 = -1;
					j6 = -1;
					test = false;
					while ((j6 < 0) && (j4 < nBonds())) {
						if (j4 != bondN) test = (getBond(j4)->at[0] == j2) || (getBond(j4)->at[1] == j2);
						if (test) {
							if (j5 < 0) j5 = j4; else j6 = j4;
						};
						test = false;
						j4++;
					}; // until(j6 < >0) or (j4 = fBond.nBonds);
					r = cosB(bondN, j5);
					r = sqrt(abs(1 - r*r));
					if (r > 0.3) {
						j7 = getBond(j5)->at[0];
						xN1 = getAtom(j7)->rx; yN1 = getAtom(j7)->ry;
						j7 = getBond(j5)->at[1];
						xN2 = getAtom(j7)->rx; yN2 = getAtom(j7)->ry;
						crossLines(x3, y3, x4, y4, xN1, yN1, xN2, yN2, rx, ry);
						crossLines(x5, y5, x6, y6, xN1, yN1, xN2, yN2, rx1, ry1);
						if (xN1 < xN2) { lx = xN1; xN1 = xN2; xN2 = lx; };
						if (yN1 < yN2) { lx = yN1; yN1 = yN2; yN2 = lx; };
						if ((rx >= (xN2 - 0.01)) && (rx <= (xN1 + 0.01)) && (ry >= (yN2 - 0.01)) && (ry <= (yN1 + 0.01))) {
							x4 = rx; y4 = ry;
						}
						else if ((rx1 >= (xN2 - 0.01)) && (rx1 <= (xN1 + 0.01)) && (ry1 >= (yN2 - 0.01)) && (ry1 <= (yN1 + 0.01))) {
							x6 = rx1; y6 = ry1;
						}; //	else return;
					};
					r = cosB(bondN, j6);
					r = sqrt(abs(1 - r*r));
					if (r > 0.3) {
						j7 = getBond(j6)->at[0];
						xN1 = getAtom(j7)->rx; yN1 = getAtom(j7)->ry;
						j7 = getBond(j6)->at[1];
						xN2 = getAtom(j7)->rx; yN2 = getAtom(j7)->ry;
						crossLines(x3, y3, x4, y4, xN1, yN1, xN2, yN2, rx, ry);
						crossLines(x5, y5, x6, y6, xN1, yN1, xN2, yN2, rx1, ry1);
						if (xN1 < xN2) { lx = xN1; xN1 = xN2; xN2 = lx; };
						if (yN1 < yN2) { lx = yN1; yN1 = yN2; yN2 = lx; };
						if ((rx >= (xN2 - 0.01)) && (rx <= (xN1 + 0.01)) && (ry >= (yN2 - 0.01)) && (ry <= (yN1 + 0.01))) {
							x4 = rx; y4 = ry;
						}
						else if ((rx1 >= (xN2 - 0.01)) && (rx1 <= (xN1 + 0.01)) && (ry1 >= (yN2 - 0.01)) && (ry1 <= (yN1 + 0.01))) {
							x6 = rx1; y6 = ry1;
						};	//	else return;
					};
				};
				if (lineType == 5)  dashedLineSVG(x3, y3, x4, y4, polygonData); else result = result + solidLineSVG(x3, y3, x4, y4);
				if ((lineType == 4) || (lineType == 5)) dashedLineSVG(x5, y5, x6, y6, polygonData); else result = result + solidLineSVG(x5, y5, x6, y6);
				if (lineType == 3) result = result + solidLineSVG(x1, y1, x2, y2);
			};
			
			break;
			case 9 :
				
				dx = abs(x1 - x2);
				dy = abs(y1 - y2);
				if (dx > dy) {
					dx = 0; dy = 1 / 120; //120 = pi/360
				}
				else {
					dx = 1 / 120; dy = 0;
				};
				dx = dx*fSVGAngle*rR;
				dy = dy*fSVGAngle*rR;
				x3 = x2 + dx; y3 = y2 + dy;
				x2 = x2 - dx; y2 = y2 - dy;
				polyDX[0] = x1; polyDY[0] = y1;
				polyDX[1] = x2; polyDY[1] = y2;
				polyDX[3] = x3; polyDY[2] = y3;
				ss = intToStr(polyDX[0]) + "," + intToStr(polyDY[0]) + "  " + intToStr(polyDX[1]) + "," + intToStr(polyDY[1]) + "  " + intToStr(polyDX[2]) + "," + intToStr(polyDY[2]);
				ss = format("<polygon fill=\"black\" stroke=\"black\" stroke-width=\"1\" points=\"{}\" />", ss);

				polygonData.push_back(ss);
			
			break;
			case 10 :
				
				s1 = (x2 - x1);
				s2 = (y2 - y1);
				r = sqrt(s1*s1 + s2*s2);
				r1 = 0.05*svgPixPerInch;
				si = trunc(r / r1);
				if (si < 5) { si = 5; r1 = r / si; };
				for (i = 0; i <= si; i++) {
					rx = r1*i*s1 / r;
					ry = r1*i*s2 / r;
					lx = r1*i / si;
					lx = lx*fSVGAngle*rR / svgPixPerInch / 4;
					xN1 = x1 + rx - lx*s2 / r;
					xN2 = x1 + rx + lx*s2 / r;
					yN1 = y1 + ry + lx*s1 / r;
					yN2 = y1 + ry - lx*s1 / r;
					result = result + solidLineSVG(xN1, yN1, xN2, yN2);
				};
			
			break;
			case 11:
				s1  = x2 - x1; 
				s2 = y2 - y1;
				r = sqrt(s1*s1 + s2*s2);
				r1 = fSVGSOR*svgPixPerInch / 2;
				si = trunc(r / r1); 
				if ((si % 2) == 1)  si = si + 1; 
				if (si == 0) si = 2;
				r = sqrt(s1*s1 + s2*s2); 
				dx = r / si;
				r1 = 2 * dx; 
				si = si / 2;
				rx = r1*s1 / r / 2; 
				ry = r1*s2 / r / 2;
				xX1 = x1 - rx; 
				yY1 = y1 - ry;
				if (si > 0) for (i = 1; i < si; i++) {
					rx = r1*i*s1 / r;
					ry = r1*i*s2 / r;
					rx1 = r1*(i + 1)*s1 / r;
					ry1 = r1*(i + 1)*s2 / r;
					xN1 = xX1 + rx - dx*s2 / r;
					xN2 = xX1 + rx1 + dx*s2 / r;
					yN1 = yY1 + ry + dx*s1 / r;
					yN2 = yY1 + ry1 - dx*s1 / r;
					result = result + solidLineSVG(xN1, yN1, xN2, yN2);
					dx = -dx;
					if (i == 1) result = result + solidLineSVG(xN1, yN1, x1, y1);
					if (i == (si - 1)) result = result + solidLineSVG(xN2, yN2, x2, y2);
				};
			break;
		}; //switch


		return result;
	};



	void TSVGMolecule::aDrawerSVG(std::vector<std::string> & dataOut, int atomN, int fontSize, const std::vector<std::string> & atomProperties, bool  arrowDrawIsotope, const std::string aNum) const {
		const std::string valDes[11] = { "(0)", "(I)", "(II)", "(III)", "(IV)", "(V)", "(VI)", "(VII)", "(VIII)",	"(IX)", "(X)" };

		int	c2, k, j;
		std::string s;
		int i;
		bool test1, test;
		double cSF, xu, yu, xu1, yu1, r1, r2, r3;
		int nH;
		int c1;
		double right, left, yout;
		TSingleAtom  aT;
		std::string asym, asymLeft;
		double rr;
		int iz;

		int nnh;
		


		r2 = -1;
		if (getAtom(atomN)->na == 0) return;
		asym = "";
		if (getAtom(atomN)->na < NELEMMCDL) asym = aSymb[getAtom(atomN)->na];
		if (asym.length() == 0) return;
		aT = *getAtom(atomN);
		test1 = (aT.na == 6) && (aT.nb == 2);
		if (test1) { //checking if linear fragment is used for the atom 
			c1 = aT.ac[0];
		  //it is necessary to draw carbon in linear fragment)
			xu = getAtom(c1)->rx - aT.rx;
			yu = getAtom(c1)->ry - aT.ry;
			c1 = aT.ac[1];
			xu1 = getAtom(c1)->rx - aT.rx;
			yu1 = getAtom(c1)->ry - aT.ry;
			cSF = sqrt(xu*xu + yu*yu)*sqrt(xu1*xu1 + yu1*yu1);
			if (cSF > 0.00001) {
				cSF = 1 + (xu*xu1 + yu*yu1) / cSF;
				if (cSF < 0.1) test1 = true; else test1 = false;
			}
			else test1 = false;
		};
		if (! test1) test1 = aT.astereo != 0;
		test = (aT.na == 6) && (aT.nv == 4) && (aT.nc == 0) && (aT.iz == 0)  &&(aNum.length()==0) && ((aT.astereo == 0) || (options.fIOPT3 == 3)) && (aT.nb != 0) && (! options.fIOPT4) && (aT.special == 0);
		test = test && (! test1);
		
		//carbon without attribute will not be drawn

		if (!((!test) || ((aT.nb <= 1) && (options.fIOPT1 == 2)))) {
			if (atomProperties.size() > atomN)	if (atomProperties[atomN].length() > 0) {
				k = round(fontSize*svgFontSmallRatio);
				yout = round(aT.ry) - k;
				left = round(aT.rx) - getTextWidthLarge(atomProperties[atomN]) / 2; 
				//   !!!!!!!!!     sVGTextOut(round(left),round(yOut),atomProperties[atomN],);
			};
			return;
		};
		nnh = 0;
		
		if (((aT.nb <= 1) && (options.fIOPT1 == 2)) || ((aT.na != 6) && ((options.fIOPT1 == 2) || (options.fIOPT1 == 3)))) {

			//implicitly defined hydrogens atoms are drawn.iOPT[1] values:
			//= 1 - no hydrigen draw. = 2 - for atoms, which have one or less neighbour, =3 - for hetero - atoms

			unitVector(atomN, xu, yu);
			nH = aT.nv;
			nH = nH - aT.currvalence - abs(aT.nc) - aT.rl;
		
			if (nH > 0) {
				if (nH > 1) s = intToStr(nH);
				if ((abs(xu) < 0.1) || (aT.nb == 1)) { // definition must be more sertain 
					xu = 0;
					if (aT.nb > 0) for (i = 0; i < aT.nb; i++) xu = xu - (getAtom(aT.ac[i])->rx - aT.rx);
					if (abs(xu) < 0.1) {
						r1 = 100000;
						for (i = 0; i < nBonds(); i++) {
							r3 = xDist(atomN, i);
							if (abs(r3) < r1) { r1 = abs(r3); r2 = r3; };
						};
						if (r2 < 0) xu = 1; else xu = -1;
					}
				};

				if (xu <= 0) { //left definition 
					nnh = nH;
				}
				else {
					nnh = -nH;
				};
			};
		};

		//if (fAtom[atomN].NA=6) and (fAtom[atomN].NB=1) then Asym:='CH3';
		unitVector(atomN, xu, yu);
		if ((abs(xu) < 0.1) || (aT.nb == 1)) { //Definition must be more sertain 
			xu = 0;
			if (aT.nb > 0) for (i=0; i<aT.nb; i++) xu = xu - (getAtom(aT.ac[i])->rx - aT.rx);
		};
		if (nnh == 0) {
			if (xu <= 0) nnh = 1000; else nnh = -1000;
		};
		iz = 0;
		if (aT.iz != 0) iz  = round(aMass[aT.na] + aT.iz);
		if (arrowDrawIsotope) iz = 0;
		s = asym;
		if (aNum.length() > 0) s = s + "(" + aNum + ")";
		svgTextOut(aT.rx, aT.ry, fontSize, nnh, aT.nc, iz, aT.rl, s, dataOut);


	};

	double TSVGMolecule::xDist(int aN, int bN) const {
		//the function finds distance between fAtom aN and bond bN(including sign).
		//	the bond is treated as segment, not as straight line!
		double d, r, x1, y1, x2, y2, x0, y0, yMin, yMax;
		double result=1e9;
		
		if ((aN == getBond(bN)->at[0]) || (aN == getBond(bN)->at[1])) return result;
		x1 = getAtom(getBond(bN)->at[0])->rx; 
		y1 = getAtom(getBond(bN)->at[0])->ry;
		x2 = getAtom(getBond(bN)->at[1])->rx;
		y2 = getAtom(getBond(bN)->at[1])->ry;

  	    x0= getAtom(aN)->rx; 
		y0 = getAtom(aN)->ry;
		if (y1 < y2) { yMin = y1; yMax = y2; }else { yMin = y2; yMax = y1; };
		r = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
		yMin = yMin - 0.1*r; yMax = yMax + 0.1*r;
		d = y2 - y1;
		if (abs(d) < 1E-8) return result;
		if ((y0 > yMin) && (y0 < yMax)) {
			r = x1 + (y0 - y1)*(x2 - x1) / d;
			result = r - x0;
		}
		return result;
	};


	std::string TSVGMolecule::arrowDrawSVG(std::vector<std::string> & polygonList) {
		std::string result = "";
		int i, n, k, lenNak, widNak;
		double rX, rY, d, x, y, a45;
		
		a45 = sqrt(2) / 2;
		for (i = 0; i < nAtoms(); i++) {
			n = getAtom(i)->iz % 3;
			if (n != 0) {
				unitVector(i, rX, rY);
				if (getAtom(i)->nb > 0) {
					k = getAtom(i)->ac[0];
					x = getAtom(i)->rx - getAtom(k)->rx;
					y = getAtom(i)->ry - getAtom(k)->ry;
					d = sqrt(x*x + y*y);
				} else {
					d = averageBondLength();
					if (d == 0) d = 30;//recommendedBondLength;
				};
				lenNak = round(d / 3);
				widNak = round(d / 12);
				if (n == 1) {
					x = rX*d;
					y = rY*d;
					drawArrowSVG(result, polygonList, getAtom(i)->rx, getAtom(i)->ry, getAtom(i)->rx + x, getAtom(i)->ry + y, lenNak, widNak);
				}
				else {
					x = rX*a45 + rY*a45;
					y = -rX*a45 + rY*a45;
					x = x*d;
					y = y*d;
					drawArrowSVG(result, polygonList, getAtom(i)->rx, getAtom(i)->ry, getAtom(i)->rx + x, getAtom(i)->ry + y, lenNak, widNak);
					x = rX*a45 - rY*a45;
					y = rX*a45 + rY*a45;
					x = x*d;
					y = y*d;
					drawArrowSVG(result, polygonList, getAtom(i)->rx, getAtom(i)->ry, getAtom(i)->rx + x, getAtom(i)->ry + y, lenNak, widNak);
				};
			};
		};

		return result;
	}


	std::string TSVGMolecule::svgSaveInternal(const std::vector<std::string> & atomProperties, const std::vector<int> & redBonds, const std::vector<std::vector<int>*> & ringList, 
		bool  arrowDrawIsotope, bool numerationDraw) {

		int w;
		std::vector<singleCenter> moleculeCenters;
	    double rr, rX, rY;
		std::vector<std::string> polygonData;
		std::vector<int> * bondList;
		std::vector<int> atomList;
		std::string s, s1;
	    int i,j,n;
		std::string	result = "";
		singleCenter bondCenter;

		if (nAtoms() == 0) return result;
		s = "";

		bondCenter.centerX = 0;
		bondCenter.centerY = 0;
		bondCenter.cSize = 0;
		moleculeCenters.reserve(nBonds());
		for (i = 0; i < nBonds(); i++) moleculeCenters.push_back(bondCenter);
		for (i = 0; i < ringList.size(); i++) {
			rX = 0; rY = 0;
			bondList = ringList[i];
			bondListToAtomList(*bondList, atomList);
			for (j = 0; j < atomList.size(); j++) {
				n = atomList[j];
				rX = rX + getAtom(n)->rx;
				rY = rY + getAtom(n)->ry;
			};
			rX = rX / atomList.size();
			rY = rY / atomList.size();
			for (j = 0; j < bondList->size(); j++) {
				n = (*bondList)[j];
				bondCenter = moleculeCenters[n];
				if (bondList->size() > bondCenter.cSize) {
					bondCenter.centerX = rX;
					bondCenter.centerY = rY;
					bondCenter.cSize = bondList->size();
					moleculeCenters[n] = bondCenter;
				};
			};
		};


		rr = averageBondLength();
		if (nBonds() == 0) rr = averageAtomDistance() / 2;

		s1= formatPrecision(fSVGLineWidth*svgPixPerInch + 1,2);

		s = "";

		for (w = 0; w < nBonds(); w++) {
			s = s + bDrawerSVG(polygonData, w, rr, moleculeCenters);
		};

		if (arrowDrawIsotope) {
			s = s + arrowDrawSVG(polygonData);
		};


		s = "<path stroke=\"black\" stroke-width=\"" + s1 + "\" d=\"" + s + "\" /" + ">\n";
		for (i=0; i<polygonData.size(); i++) s = s + polygonData[i] + "\n";

		polygonData.clear();
		n = 0;
		for (w = 0; w < nAtoms(); w++) {
			if (getAtom(w)->na != ID_ZVEZDA) {
				n++;
				s1 = std::to_string(n);
			}	else s1 = "";
			if (!numerationDraw) s1 = "";
			aDrawerSVG(polygonData, w, fSVGFontHeight, atomProperties, arrowDrawIsotope, s1);
		};
		for (i=0; i<polygonData.size(); i++) s = s + polygonData[i]+"\n";
//		fCharSeparate: = false;
		result = s;
		return result;
	};

	std::string TSVGMolecule::getSVG(int bmWidth, int bmHeight, const std::vector<std::vector<int>*> & ringList, std::vector<std::string> & outBuffer, bool  arrowDrawIsotope, bool numerationDraw) {
		double dNorm, dd, dScale, xMin, xMax, yMin, yMax;
		double xCorr, yCorr;
		int i;
		int bmWInternal, bmHInternal;
		std::string s;
		std::string s1, s2;
		std::vector<std::string> parList, atomProperties;
		std::vector<int> redBonds;

		std::string result = "";

		if (nAtoms() == 0) return result;
		bmWInternal = bmWidth - 2 * svgMarginXPix;
		bmHInternal = bmHeight - 2 * svgMarginYPix;
		if ((bmWInternal < 30) || (bmHInternal < 20)) return result;


		//defineConn;
		//determineFormula;
		//allAboutCycles;

		dNorm = averageBondLength();
		if (nBonds() == 0) dNorm = 0;
		if (dNorm <= 0) dNorm = averageAtomDistance() / 2;
		if (dNorm > 0)  for (i = 0; i < nAtoms(); i++) {  //Pixels coordinates
			getAtom(i)->rx = getAtom(i)->rx * fRecommendedBondLengthPix / dNorm;
			getAtom(i)->ry = getAtom(i)->ry * fRecommendedBondLengthPix / dNorm;
		};
		xMin = getAtom(0)->rx;
		xMax = getAtom(0)->rx;
		yMin = getAtom(0)->ry;
		yMax = getAtom(0)->ry;
		for (i = 0; i < nAtoms(); i++) {
			dd = getAtom(i)->rx;
			if (dd < xMin) xMin = dd;
			if (dd > xMax) xMax = dd;
			dd = getAtom(i)->ry;
			if (dd < yMin) yMin = dd;
			if (dd > yMax) yMax = dd;
		};
		//scale and shift
		if (dNorm == 0) {
			//single atom
			getAtom(0)->rx = svgMarginXPix + bmWInternal / 2;
			getAtom(0)->ry = svgMarginYPix + bmHInternal / 2;
			dScale = 1;
		}
		else {
			//shift atoms
			xCorr = svgMarginXPix;
			yCorr = svgMarginYPix;
			if (((xMax - xMin) <= bmWInternal) && ((yMax - yMin) <= bmHInternal)) {
				//MUST NOT be norm-all bond lengthes were calculated accoedinly recommented bond lenght
				xCorr = xCorr + (bmWInternal - (xMax - xMin)) / 2;
				yCorr = yCorr + (bmHInternal - (yMax - yMin)) / 2;
				for (i = 0; i < nAtoms(); i++) {
					getAtom(i)->rx = (getAtom(i)->rx - xMin) + xCorr;
					getAtom(i)->ry = (getAtom(i)->ry - yMin) + yCorr;
				};
				dScale = 1;
			}
			else {
				//MUST NOT be norm-all bond lengthes were calculated accoedinly recommented bond lenght
				if (((xMax - xMin) / bmWInternal) > ((yMax - yMin) / bmHInternal)) {
					//X range is too high- squizeeng by X
					dScale = (xMax - xMin) / bmWInternal;        //small scale  - will be divider...
					xCorr = svgMarginXPix;
					yCorr = (yMax - yMin) / dScale;
					yCorr = (bmHeight - yCorr) / 2;
				}
				else {
					dScale = (yMax - yMin) / bmHInternal;
					yCorr = svgMarginYPix;
					xCorr = (xMax - xMin) / dScale;
					xCorr = (bmWidth - xCorr) / 2;
				};
				for (i=0; i<nAtoms(); i++){
					getAtom(i)->rx = (getAtom(i)->rx - xMin) + xCorr;
					getAtom(i)->ry = (getAtom(i)->ry - yMin) + yCorr;
				};
			};
		};
		parList.clear();
		parList.push_back(std::to_string(bmWidth));
		parList.push_back(std::to_string(bmHeight));
		parList.push_back(std::to_string(bmWidth));
		parList.push_back(std::to_string(bmHeight));

		s = format("<svg width=\"{}\" height=\"{}\" viewbox=\"0 0 {} {}\"  xmlns=\"http://www.w3.org/2000/svg\">", parList); //[bmWidth, bmHeight, bmWidth, bmHeight]);
		outBuffer.push_back(s);

		parList.clear();
		parList.push_back(svgDefaultAtomBackColor);
		parList.push_back(svgDefaultAtomBackColor);
		parList.push_back(svgDefaultAtomBackColor);
		parList.push_back(svgDefaultAtomBackColor);
		parList.push_back(svgDefaultAtomFontColor);

		s = format("<style type=\"text/css\"><![CDATA[ circle { stroke: {}; fill: {}; stroke-width: 1.0;} rect { stroke: {}; fill: {}; stroke-width: 1.0;} text { stroke-width: 0; fill: {}} line { stroke: black;} ]]></style>", parList); //[SVGDefaultAtomBackColor, SVGDefaultAtomBackColor, SVGDefaultAtomBackColor, SVGDefaultAtomBackColor, SVGDefaultAtomFontColor]);
		outBuffer.push_back(s);

		s = "";
		if (dScale != 1) {
			dScale = 1 / dScale;
			s= formatPrecision(dScale,5);
			s = format("transform=\"scale({})\"", s);
		};
		parList.clear();
		parList.push_back(std::to_string(fSVGFontHeight));
		parList.push_back(fSVGFontFamilyName);
		parList.push_back(s);
		s = format("<g font-size=\"{}\" font-family=\"{}\" stroke=\"black\" stroke-width=\"1\" fill=\"black\" {}>", parList); //[RecommendedFontHeight, fontFamilyName, s]);
		outBuffer.push_back(s);

		parList.clear();
		parList.push_back(std::to_string(bmWidth));
		parList.push_back(std::to_string(bmHeight));
		s = format("<rect x=\"0\" y=\"0\" width=\"{}\" height=\"{}\" />", parList);
		outBuffer.push_back(s);



		s = svgSaveInternal(atomProperties, redBonds, ringList, arrowDrawIsotope, numerationDraw);
		outBuffer.push_back(s);
		outBuffer.push_back("</g>");
		/*
		if (true) {//(numerationDraw) {
			parList[0] = std::to_string(fSVGFontHeight * 2 / 3);
			s = format("<g font-size=\"{}\" font-family=\"{}\" stroke=\"black\" stroke-width=\"1\" fill=\"black\" {}>", parList);
			//outBuffer.push_back(s);
			for (i = 0; i < nAtoms(); i++) {
				s= "("+std::to_string(i+1)+")";
				parList.clear();
				s1 = formatPrecision(getAtom(i)->rx, 1);
				s2 = formatPrecision(getAtom(i)->ry, 1);
				parList.push_back(s1);
				parList.push_back(s2);
				parList.push_back(s);

				s = format("<text x=\"{}\" y=\"{}\">{}</text>", parList);  //[X-textWidth div 2,Y+(3*textHeight)div 2,s]);
				outBuffer.push_back(s);

			}
			//outBuffer.push_back("</g>");
		}
		*/
		outBuffer.push_back("</svg>");

		return result;
	};


};//namespace