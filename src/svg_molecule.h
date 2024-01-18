#ifndef _HDR_SVG_MOLECULE_
#define _HDR_SVG_MOLECULE_

#include <string>
#include<vector>
#include "simple_molecule.h"
#include "mol_struct_common.h"

namespace MolStruct {

#define SVG_IMAGE_WIDTH 400
#define SVG_IMAGE_HEIGHT 300

struct singleCenter {
	double centerX, centerY;
    int cSize;
};


class TSVGMolecule: public TSimpleMolecule {
protected:
	void dotLineSVG(double x1, double y1, double x2, double y2, std::vector<std::string> & outBuffer) const;
	void dashedLineSVG(double x1, double y1, double x2, double y2, std::vector<std::string> & outBuffer) const;
	std::string  solidLineSVG(double x1, double y1, double x2, double y2) const;
	void svgTextOut(double x, double y, double fontSize, int nH, int charge, int iz, int rl, const std::string sData, std::vector<std::string> & outBuffer) const;
	std::string bDrawerSVG(std::vector<std::string> & polygonData, int bondN, double ml, const std::vector<singleCenter> & center) const;
	void aDrawerSVG(std::vector<std::string> & dataOut, int atomN, int fontSize, const std::vector<std::string> & atomProperties, bool  arrowDrawIsotope, const std::string aNum) const;
	std::string svgSaveInternal(const std::vector<std::string> & atomProperties, const std::vector<int> & redBonds, const std::vector<std::vector<int>*> & ringList, bool  arrowDrawIsotope, bool numerationDraw);
	int getTextWidthLarge(const std::string & data) const;
	int getTextWidthSmall(const std::string & data) const;
	double cosB(int b1, int b2) const;
	double xDist(int aN, int bN) const;
	std::string  arrowDrawSVG(std::vector<std::string> & polygonList);
	double	_getMW() const;

public:
	const static int		svgMarginXPix = 15;
	const static int		svgMarginYPix = 10;
	const static int		svgCenterRatio = 20;  //cyclic double bond
	const static int		svgCenterRatio2 = 18; //triple bond
	const static int		svgCenterRatio1 = 12; //acyclic double bond
	const static int		svgPixPerInch = 96;


	double			svgFontSmallRatio;
	std::string		svgDefaultAtomBackColor;
	std::string		svgDefaultAtomFontColor;
	int				svgDefaultAtomBox;  //=0-box, =1 - circle
	double			svgDeltaR;
	double			fSVGLineWidth;
    int				fSVGFontHeight;
	std::string		fSVGFontFamilyName;
	int				fRecommendedBondLengthPix;
	double			fMinLineDistance, fMaxLineDistance;
	double		    fSVGAngle; //angle up / down bond
	double			fSVGSOR; //amplitude of Either bond inch
	cfIOPT			options;


	std::string    getSVG(int bmWidth, int bmHeight, const std::vector<std::vector<int>*> & ringList, std::vector<std::string> & outBuffer, bool  arrowDrawIsotope, bool numerationDraw);
	std::string cb_getSVG(int bmWidth, int bmHeight, const std::vector<std::vector<int>*> & ringList, std::vector<std::string> & outBuffer, bool  arrowDrawIsotope, bool numerationDraw);
	void 		cb_svgSave(std::vector<std::string> & svgData, bool numerationOutput, int imgWidth, int imgHeight);
	TSVGMolecule() : TSimpleMolecule() {
		svgDefaultAtomBackColor = "white";// RGB(255, 255, 255)";
		svgDefaultAtomFontColor = "black";// RGB(0, 0, 0)";
		svgDefaultAtomBox = 0;  //=0-box, =1 - circle
		svgDeltaR = 1;
		svgFontSmallRatio = 2.0 / 3.0;
		fSVGFontHeight=16;
		fSVGFontFamilyName = "Arial";
		fSVGLineWidth = 0.010;
		fRecommendedBondLengthPix = 40;
		fMinLineDistance = 1.0/96;
		fMaxLineDistance = 6.0/96;
		fSVGAngle = 10; 
		fSVGSOR = 0.05; 


	};
};
void svgPolymerSave(const TSimpleMolecule & sm, std::vector<std::string> & svgData, bool numerationOutput, int imgWidth= SVG_IMAGE_WIDTH, int imgHeight= SVG_IMAGE_HEIGHT);

} // namespace MolStruct
#endif
