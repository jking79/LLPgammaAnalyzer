// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// basic C++ types
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cmath>

//Root includes
#include "TTree.h"

#ifndef KUCMSBranchHeader
#define KUCMSBranchHeader

typedef unsigned int uInt;

class KUCMSBranch {

	public:

	// BranchType List
	// --------------------------
	// vector<vector<uInt> 	VVUI
	// vector<uInt>			VUI 
	// vector<int>			VI
	// vector<float>		VF
	// vector<string>		VSTR
	// uInt					UI
	// int					SI
	// float				FL
	// string				STR
	
	enum BType{ VVUI, VUI, VI, VF, VSTR, UI, SI, FL, STR };

	KUCMSBranch();
	KUCMSBranch( BType type, std::string name, std::string doc );
	//~KUCMSBranch();
	
	void initBranch( TTree* fOutTree );
	void clearBranch();

	void fillBranch( std::vector<uInt> val );
    void fillBranch( uInt val );
    void fillBranch( int val );
    void fillBranch( float val );
    void fillBranch( std::string val );

	void getValue( std::vector<uInt>& val, uInt index );
    void getValue( uInt& val, uInt index );
    void getValue( int& val, uInt index );
    void getValue( float& val, uInt index );
    void getValue( std::string& val, uInt index );

    std::vector<uInt> getVUIValue( uInt index );
    uInt getUIValue( uInt index );
    int getSIValue( uInt index );
    float getFLValue( uInt index );
    std::string getSTRValue( uInt index );

	float getMaxValue();
	uInt getLeadIndex();
	uInt getSubLeadIndex();

	private:

	BType BranchType;
	std::string BranchName;
	std::string BranchDoc;
	std::vector<std::vector<uInt>> VVUIBranch;
	std::vector<uInt> VUIBranch;
	std::vector<int> VIBranch;
    std::vector<float> VFBranch;
    std::vector<std::string> VSTRBranch;
	uInt UIBranch;
	int SIBranch;
	float FLBranch;
	std::string STRBranch;

	float max( std::vector<float> x );
	uInt leadIdx( std::vector<float> x );
	uInt subldIdx( std::vector<float> x, uInt ldx );

};//<<>>class KUCMSBranch

KUCMSBranch::KUCMSBranch(){}

KUCMSBranch::KUCMSBranch( KUCMSBranch::BType type, std::string name, std::string doc = "" ):

	BranchType(type),
	BranchName(name),
	BranchDoc(doc)

{}//<<>>KUCMSBranch::KUCMSBranch( BType type, std::string name, std::string doc )


void KUCMSBranch::initBranch( TTree* fOutTree ){

	switch( BranchType ){

		case VVUI	:	fOutTree->Branch( BranchName.c_str(), &this->VVUIBranch )->SetTitle( BranchDoc.c_str() ); break;
		case VUI	:	fOutTree->Branch( BranchName.c_str(), &this->VUIBranch )->SetTitle( BranchDoc.c_str() ); break;
		case VI		:	fOutTree->Branch( BranchName.c_str(), &this->VIBranch )->SetTitle( BranchDoc.c_str() ); break;
		case VF		:   fOutTree->Branch( BranchName.c_str(), &this->VFBranch )->SetTitle( BranchDoc.c_str() ); break;
		case VSTR	:   fOutTree->Branch( BranchName.c_str(), &this->VSTRBranch )->SetTitle( BranchDoc.c_str() ); break;
		case UI		:   fOutTree->Branch( BranchName.c_str(), &this->UIBranch )->SetTitle( BranchDoc.c_str() ); break;
		case SI		:   fOutTree->Branch( BranchName.c_str(), &this->SIBranch )->SetTitle( BranchDoc.c_str() ); break;
		case FL		:   fOutTree->Branch( BranchName.c_str(), &this->FLBranch )->SetTitle( BranchDoc.c_str() ); break;
		case STR	:   fOutTree->Branch( BranchName.c_str(), &this->STRBranch )->SetTitle( BranchDoc.c_str() ); break;
		default	: std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in setBranch!!!! " << std::endl;

	}//<<>>switch( BranchType )

}//<<>>void setBranch( TTree* fOutTree )

void KUCMSBranch::clearBranch(){

    switch( BranchType ){

        case VVUI   : VVUIBranch.clear(); break;
        case VUI    : VUIBranch.clear(); break;
        case VI     : VIBranch.clear(); break;
        case VF     : VFBranch.clear(); break;
        case VSTR   : VSTRBranch.clear(); break;
        case UI     : UIBranch = 0; break;
        case SI     : SIBranch = -99999; break;
        case FL     : FLBranch = -99999.f; break;
        case STR    : STRBranch = ""; break;
		default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in clearBranch!!!! " << std::endl;    

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::clearBranch()

void KUCMSBranch::fillBranch( std::vector<uInt> val ){

    switch( BranchType ){

        case VVUI   : VVUIBranch.push_back(val); break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( std::vector<uInt> val )

void KUCMSBranch::fillBranch( uInt val ){

    switch( BranchType ){

        case VUI    : VUIBranch.push_back(val); break;
        case UI     : UIBranch = val; break;
   		default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;
	 
    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( uInt val )

void KUCMSBranch::fillBranch( int val ){

    switch( BranchType ){

        case VI		: VIBranch.push_back(val); break;
        case SI    	: SIBranch = val; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( int val )

void KUCMSBranch::fillBranch( float val ){

    switch( BranchType ){

        case VF     : VFBranch.push_back(val); break;
        case FL     : FLBranch = val; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( float val )

void KUCMSBranch::fillBranch( std::string val ){

    switch( BranchType ){

        case VSTR   : VSTRBranch.push_back(val); break;
        case STR    : STRBranch = val; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with fillBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::fillBranch( std::string val )

void KUCMSBranch::getValue( std::vector<uInt>& val, uInt index = 0 ){

    switch( BranchType ){

        case VVUI   : val = VVUIBranch[index]; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )	

}//<<>>void KUCMSBranch::getBranch( std::vector<uInt>& val, uInt index = 0 )

void KUCMSBranch::getValue( uInt& val, uInt index = 0 ){

    switch( BranchType ){

        case VUI    : val = VUIBranch[index]; break;
        case UI     : val = UIBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( uInt& val, uInt index = 0 )

void KUCMSBranch::getValue( int& val, uInt index = 0 ){

    switch( BranchType ){

        case VI     : val = VIBranch[index]; break;
        case SI     : val = SIBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( int& val, uInt index = 0 )

void KUCMSBranch::getValue( float& val, uInt index = 0 ){

    switch( BranchType ){

        case VF     : val = VFBranch[index]; break;
        case FL     : val = FLBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( float& val, uInt index = 0 )

void KUCMSBranch::getValue( std::string& val, uInt index = 0 ){

    switch( BranchType ){

        case VSTR   : val = VSTRBranch[index]; break;
        case STR    : val = STRBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( std::string& val, uInt index = 0 )

std::vector<uInt> KUCMSBranch::getVUIValue( uInt index = 0 ){ std::vector<uInt> val; getValue( val, index ); return val; }
uInt KUCMSBranch::getUIValue( uInt index = 0 ){ uInt val; getValue( val, index ); return val; }
int KUCMSBranch::getSIValue( uInt index = 0 ){ int val; getValue( val, index ); return val; }
float KUCMSBranch::getFLValue( uInt index = 0 ){ float val; getValue( val, index ); return val; }
std::string KUCMSBranch::getSTRValue( uInt index = 0 ){ std::string val; getValue( val, index ); return val; }

float KUCMSBranch::max( std::vector<float> x ){float m(x[0]); for(auto ix : x ){ if( ix > m ) m = ix; } return m;}

uInt KUCMSBranch::leadIdx( std::vector<float> x ){

	float m(x[0]); uInt idx(0), it(0); 
	if( x.size() == 0 ) return -9; 
	if( x.size() == 1 ) return 0; 
	for(auto ix : x ){ if( ix > m ){ m = ix; idx = it; } it++; } 
	return idx;

}//<<>>uInt KUCMSBranch::leadIdx( std::vector<float> x )

uInt KUCMSBranch::subldIdx( std::vector<float> x, uInt ldx ){

	float m(x[0]); uInt idx(0), it(0); 
	if( x.size() <= 1 ) return -9; 
	if( ldx == 0 ){ m = x[1]; idx = 1;} 
	for(auto ix : x ){ if( ix > m && ix < x[ldx] ){ m = ix; idx = it; } it++; } 
	return idx;

}//<<>>uInt KUCMSBranch::subldIdx( std::vector<float> x, uInt ldx )


float KUCMSBranch::getMaxValue(){

	switch( BranchType ){

        case VF     : return max( VFBranch );
        case FL     : return FLBranch;
        default : 
			std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in getMaxValue!!!! " << std::endl;
			return -9999.f;
		
    }//<<>>switch( BranchType )

}//<<>>float KUCMSBranch::getMaxValue()

uInt KUCMSBranch::getLeadIndex(){

    switch( BranchType ){

        case VF     : return leadIdx( VFBranch );
        default : 
            std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in getLeadIndex!!!! " << std::endl;
            return 9999;

    }//<<>>switch( BranchType )

}//<<>>uInt KUCMSBranch::getLeadIndex()

uInt KUCMSBranch::getSubLeadIndex(){

    switch( BranchType ){

        case VF     : return subldIdx( VFBranch, leadIdx( VFBranch ) );
        default :
            std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType error in getSubLeadIndex!!!! " << std::endl;
            return 9999;

    }//<<>>switch( BranchType )

}//<<>>uInt KUCMSBranch::getSubLeadIndex()

// -------------------------------------------------------------------------------------------------
// -----------  Manager for KUCMS Branches ---------------------------------------------------------
// -------------------------------------------------------------------------------------------------

class KUCMSBranchManager {

	public:

	void makeBranch( std::string key, std::string name, KUCMSBranch::BType type, std::string doc );
    void makeBranch( std::string name, KUCMSBranch::BType type, std::string doc );
    void initBranches( TTree* fOutTree );
	void clearBranches();

    void fillBranch( std::string key, std::vector<uInt> val );
	void fillBranch( std::string key, uInt val );
    void fillBranch( std::string key, int val );
    void fillBranch( std::string key, float val );
    void fillBranch( std::string key, std::string val );

	void getBranchValue( std::string key, std::vector<uInt>& val, uInt index );
    void getBranchValue( std::string key, uInt& val, uInt index );
    void getBranchValue( std::string key, int& val, uInt index );
    void getBranchValue( std::string key, float& val, uInt index );
    void getBranchValue( std::string key, std::string& val, uInt index );

    std::vector<uInt> getVUIBranchValue( std::string key, uInt index );
    uInt getUIBranchValue( std::string key, uInt index );
	int getSIBranchValue( std::string key, uInt index );
    float getFLBranchValue( std::string key, uInt index );
    std::string getSTRBranchValue( std::string key, uInt index );

	float getBranchMaxValue( std::string key );
	uInt getBranchLeadIndex( std::string key );
	uInt getBranchSubLeadIndex( std::string key );

	private:

    std::map< std::string, KUCMSBranch > theBranches;
	bool valid( std::string key );

};//<<>>class KUCMSBranchManager 

void KUCMSBranchManager::makeBranch( std::string key, std::string name, KUCMSBranch::BType type, std::string doc = "" ){

	KUCMSBranch newBranch( type, name, doc );
	theBranches[key] = newBranch;
	
}//<<>>void KUCMSBranchManager::makeBranch( std::string key, std::string name, KUCMSBranch::BT type, std::string doc )

void KUCMSBranchManager::makeBranch( std::string name, KUCMSBranch::BType type, std::string doc = "" ){

    KUCMSBranch newBranch( type, name, doc );
    theBranches[name] = newBranch;

}//<<>>void KUCMSBranchManager::makeBranch( std::string name, KUCMSBranch::BT type, std::string doc )

bool KUCMSBranchManager::valid( std::string key ){ if( theBranches.find(key) == theBranches.end() ){ 
		std::cout << " -- Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; return false; } else return true;}

void KUCMSBranchManager::clearBranches(){ for( auto & branch : theBranches ){ (branch.second).clearBranch();}}
void KUCMSBranchManager::initBranches( TTree* fOutTree ){ for( auto & branch : theBranches ){ branch.second.initBranch( fOutTree );}}

void KUCMSBranchManager::fillBranch( std::string key, std::vector<uInt> val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, uInt val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, int val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, float val ){ if(valid(key)) theBranches[key].fillBranch( val );}
void KUCMSBranchManager::fillBranch( std::string key, std::string val ){ if(valid(key)) theBranches[key].fillBranch( val );}

void KUCMSBranchManager::getBranchValue( std::string key, std::vector<uInt>& val, uInt index = 0 ){
							if(valid(key)) theBranches[key].getValue(val,index); else { std::vector<uInt> c{0}; val = c;}}
void KUCMSBranchManager::getBranchValue( std::string key, uInt& val, uInt index = 0 ){ 
							if(valid(key)) theBranches[key].getValue( val, index ); else val = 0;}
void KUCMSBranchManager::getBranchValue( std::string key, int& val, uInt index = 0 ){ 
							if(valid(key)) theBranches[key].getValue( val, index ); else val = 0;}
void KUCMSBranchManager::getBranchValue( std::string key, float& val, uInt index = 0 ){ 
							if(valid(key)) theBranches[key].getValue( val, index ); else val = 0;}
void KUCMSBranchManager::getBranchValue( std::string key, std::string& val, uInt index = 0 ){ 
							if(valid(key)) theBranches[key].getValue( val, index ); }
std::vector<uInt> KUCMSBranchManager::getVUIBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getVUIValue(index); else { std::vector<uInt> c{0}; return c; }}
uInt KUCMSBranchManager::getUIBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getUIValue(index); else return 0;}
int KUCMSBranchManager::getSIBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getSIValue(index); else return 0;}
float KUCMSBranchManager::getFLBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getFLValue(index); else return 0;}
std::string KUCMSBranchManager::getSTRBranchValue( std::string key, uInt index = 0 ){ 
							if(valid(key)) return theBranches[key].getSTRValue(index); else return "";}

float KUCMSBranchManager::getBranchMaxValue( std::string key ){ if(valid(key)) return theBranches[key].getMaxValue(); else return 0;}
uInt KUCMSBranchManager::getBranchLeadIndex( std::string key ){ if(valid(key)) return theBranches[key].getLeadIndex();else return 0;}
uInt KUCMSBranchManager::getBranchSubLeadIndex( std::string key ){ if(valid(key)) return theBranches[key].getSubLeadIndex(); else return 0;}

#endif
