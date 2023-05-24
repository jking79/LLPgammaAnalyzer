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

	// BranchType vector<vector<uInt>,
	//            vector<uInt>, vector<int>, vector<float>
	//            uInt, int, float
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
	void getBranch( std::vector<uInt>& val, uInt index );
    void getBranch( uInt& val, uInt index );
    void getBranch( int& val, uInt index );
    void getBranch( float& val, uInt index );
    void getBranch( std::string& val, uInt index );

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

};//<<>>class KUCMSBranch

KUCMSBranch::KUCMSBranch(){}

KUCMSBranch::KUCMSBranch( KUCMSBranch::BType type, std::string name, std::string doc = "" ):

	BranchType(type),
	BranchName(name),
	BranchDoc(doc)

{}//<<>>KUCMSBranch::KUCMSBranch( BType type, std::string name, std::string doc )


void KUCMSBranch::initBranch( TTree* fOutTree ){

	switch( BranchType ){

		case VVUI	:	fOutTree->Branch( BranchName.c_str(), &this->VVUIBranch ); break;
		case VUI	:	fOutTree->Branch( BranchName.c_str(), &this->VUIBranch ); break;
		case VI		:	fOutTree->Branch( BranchName.c_str(), &this->VIBranch ); break;
		case VF		:   fOutTree->Branch( BranchName.c_str(), &this->VFBranch ); break;
		case VSTR	:   fOutTree->Branch( BranchName.c_str(), &this->VSTRBranch ); break;
		case UI		:   fOutTree->Branch( BranchName.c_str(), &this->UIBranch ); break;
		case SI		:   fOutTree->Branch( BranchName.c_str(), &this->SIBranch ); break;
		case FL		:   fOutTree->Branch( BranchName.c_str(), &this->FLBranch ); break;
		case STR	:   fOutTree->Branch( BranchName.c_str(), &this->STRBranch ); break;
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

void KUCMSBranch::getBranch( std::vector<uInt>& val, uInt index = 0 ){

    switch( BranchType ){

        case VVUI   : val = VVUIBranch[index]; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )	

}//<<>>void KUCMSBranch::getBranch( std::vector<uInt>& val, uInt index = 0 )

void KUCMSBranch::getBranch( uInt& val, uInt index = 0 ){

    switch( BranchType ){

        case VUI    : val = VUIBranch[index]; break;
        case UI     : val = UIBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( uInt& val, uInt index = 0 )

void KUCMSBranch::getBranch( int& val, uInt index = 0 ){

    switch( BranchType ){

        case VI     : val = VIBranch[index]; break;
        case SI     : val = SIBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( int& val, uInt index = 0 )

void KUCMSBranch::getBranch( float& val, uInt index = 0 ){

    switch( BranchType ){

        case VF     : val = VFBranch[index]; break;
        case FL     : val = FLBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( float& val, uInt index = 0 )

void KUCMSBranch::getBranch( std::string& val, uInt index = 0 ){

    switch( BranchType ){

        case VSTR   : val = VSTRBranch[index]; break;
        case STR    : val = STRBranch; break;
        default : std::cout << " -- KUCMSBranch " << BranchName << " Error : BranchType mismatch with getBranch type!!!! " << std::endl;

    }//<<>>switch( BranchType )

}//<<>>void KUCMSBranch::getBranch( std::string& val, uInt index = 0 )

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
	void getBranch( std::string key, std::vector<uInt>& val, uInt index );
    void getBranch( std::string key, uInt& val, uInt index );
    void getBranch( std::string key, int& val, uInt index );
    void getBranch( std::string key, float& val, uInt index );
    void getBranch( std::string key, std::string& val, uInt index );

	private:

    std::map< std::string, KUCMSBranch > theBranches;

};//<<>>class KUCMSBranchManager 

void KUCMSBranchManager::makeBranch( std::string key, std::string name, KUCMSBranch::BType type, std::string doc = "" ){

	KUCMSBranch newBranch( type, name, doc );
	theBranches[key] = newBranch;
	
}//<<>>void KUCMSBranchManager::makeBranch( std::string key, std::string name, KUCMSBranch::BT type, std::string doc )

void KUCMSBranchManager::makeBranch( std::string name, KUCMSBranch::BType type, std::string doc = "" ){

    KUCMSBranch newBranch( type, name, doc );
    theBranches[name] = newBranch;

}//<<>>void KUCMSBranchManager::makeBranch( std::string name, KUCMSBranch::BT type, std::string doc )

void KUCMSBranchManager::clearBranches(){ for( auto & branch : theBranches ){ (branch.second).clearBranch(); }}
void KUCMSBranchManager::initBranches( TTree* fOutTree ){ for( auto & branch : theBranches ){ branch.second.initBranch( fOutTree ); }}
void KUCMSBranchManager::fillBranch( std::string key, std::vector<uInt> val ){ theBranches[key].fillBranch( val ); }
void KUCMSBranchManager::fillBranch( std::string key, uInt val ){ theBranches[key].fillBranch( val ); }
void KUCMSBranchManager::fillBranch( std::string key, int val ){ theBranches[key].fillBranch( val ); }
void KUCMSBranchManager::fillBranch( std::string key, float val ){ theBranches[key].fillBranch( val ); }
void KUCMSBranchManager::fillBranch( std::string key, std::string val ){ theBranches[key].fillBranch( val ); }
void KUCMSBranchManager::getBranch( std::string key, std::vector<uInt>& val, uInt index = 0 ){ theBranches[key].getBranch( val, index ); }
void KUCMSBranchManager::getBranch( std::string key, uInt& val, uInt index = 0 ){ theBranches[key].getBranch( val, index ); }
void KUCMSBranchManager::getBranch( std::string key, int& val, uInt index = 0 ){ theBranches[key].getBranch( val, index ); }
void KUCMSBranchManager::getBranch( std::string key, float& val, uInt index = 0 ){ theBranches[key].getBranch( val, index ); }
void KUCMSBranchManager::getBranch( std::string key, std::string& val, uInt index = 0 ){ theBranches[key].getBranch( val, index ); }

#endif
