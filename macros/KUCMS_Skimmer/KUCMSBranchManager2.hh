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

//KUCMS includes
#include "KUCMSItemManager.hh"

//Root includes
#include "TTree.h"

#ifndef KUCMSBranchHeader
#define KUCMSBranchHeader

typedef unsigned int uInt;

template < class T >
class ScalarBranch : public Item<T> {

    public:
    
    void attachBranch( TTree* fOutTree ){ 
           fOutTree->Branch( Item<T>::iName.c_str(), &this->Item<T>::iValue )->SetTitle( Item<T>::iDoc.c_str() ); };

};//<<>> class Branch : Item<T> 


template < class T >
class VectorBranch : public VectorItem<T> {

    public:

    void attachBranch( TTree* fOutTree ){
           fOutTree->Branch( VectorItem<T>::iName.c_str(), &this->VectorItem<T>::iVector )->SetTitle( VectorItem<T>::iDoc.c_str() ); };

};//<<>> class Branch : Item<T> 

enum BType{ VVUINT, VUINT, VINT, VFLOAT, VSTR, VBOOL, UINT, INT, FLOAT, STR, BOOL };

class KUCMSBranchBase {

    public:

    KUCMSBranchBase(){};
    virtual ~KUCMSBranchBase(){};

    virtual void attachBranch( TTree* fOutTree ){};
    virtual void clear(){};

    virtual void fill( std::vector<uInt> val ){};
    virtual void fill( uInt val ){};
    virtual void fill( int val ){};
    virtual void fill( float val ){};
    virtual void fill( std::string val ){};
    virtual void fill( bool val ){};

};//<<>> class KUCMSBranchBase

template <class T>
class KUCMSBranch : public KUCMSBranchBase {

    public:

    KUCMSBranch( std::string name, BType type, std::string doc  );
    ~KUCMSBranch(){};

    void attachBranch( TTree* fOutTree );
    void clear();
    void fill( T val );
    
    private:

    BType vartype;
    ScalarBranch<T> bscalar;
    VectorBranch<T> bvector;
 
};//class Branch 

template <class T>
KUCMSBranch<T>::KUCMSBranch( std::string name, BType type, std::string doc  ){

    vartype = type;
    ( vartype >= UINT ) ? bscalar.make( name, doc ) : bvector.make( name, doc );

}//<<>>KUCMSBranch<T>::KUCMSBranch( std::string name, BType type, std::string doc  )

template <class T>
void KUCMSBranch<T>::attachBranch( TTree* fOutTree ){ 

    ( vartype >= UINT ) ? bscalar.attachBranch( fOutTree ) : bvector.attachBranch( fOutTree );

}//<<>>void KUCMSBranch<T>::attachBranch()

template <class T>
void KUCMSBranch<T>::clear(){

    ( vartype >= UINT ) ? bscalar.clear() : bvector.clear();

}//<<>>void KUCMSBranch<T>::clear()

template <class T>
void KUCMSBranch<T>::fill( T val ){

    ( vartype >= UINT ) ? bscalar.fill( val ) : bvector.fill( val );

}//<<>>void KUCMSBranch<T>::fill( T val )

class KUCMSBranchManager {

    public:

    KUCMSBranchManager(){};
    ~KUCMSBranchManager();

    void makeBranch( std::string key, std::string name, BType type, std::string doc = "" );
    void makeBranch( std::string name, BType type, std::string doc = "" );
    void attachBranches( TTree* fOutTree );
    void clearBranches();

    void fillBranch( std::string key, std::vector<uInt> val );
    void fillBranch( std::string key, uInt val );
    void fillBranch( std::string key, int val );
    void fillBranch( std::string key, float val );
    void fillBranch( std::string key, std::string val );
    void fillBranch( std::string key, bool val );

    private:

    std::map< std::string, KUCMSBranchBase* > theBranches;
    bool valid( std::string key );

};//<<>>class KUCMSBranchManager 

KUCMSBranchManager::~KUCMSBranchManager(){ 

    for( auto & branch : theBranches ){ delete branch.second; }

}//<<>>KUCMSBranchManager::~KUCMSBranchManager()

void KUCMSBranchManager::makeBranch( std::string key, std::string name, BType type, std::string doc ){

    switch( type ){

        case VVUINT : theBranches[key] = new KUCMSBranch<std::vector<uInt>>( name, type, doc ); break;
        case VUINT  : theBranches[key] = new KUCMSBranch<uInt>( name, type, doc ); break;
        case VINT   : theBranches[key] = new KUCMSBranch<int>( name, type, doc ); break;
        case VFLOAT : theBranches[key] = new KUCMSBranch<float>( name, type, doc ); break;
        case VSTR   : theBranches[key] = new KUCMSBranch<std::string>( name, type, doc ); break;
        case VBOOL  : theBranches[key] = new KUCMSBranch<bool>( name, type, doc ); break;
        case UINT   : theBranches[key] = new KUCMSBranch<uInt>( name, type, doc ); break;
        case INT    : theBranches[key] = new KUCMSBranch<int>( name, type, doc ); break;
        case FLOAT  : theBranches[key] = new KUCMSBranch<float>( name, type, doc ); break;
        case STR    : theBranches[key] = new KUCMSBranch<std::string>( name, type, doc ); break;
        case BOOL   : theBranches[key] = new KUCMSBranch<bool>( name, type, doc ); break;
        default : std::cout << " -- KUCMSBranch " << name << " Error : BranchType error in makeBranch!!!! " << std::endl;

    }//<<>>switch( type )
    
}//<<>>void KUCMSBranchManager::makeBranch( std::string key, std::string name, BT type, std::string doc )

void KUCMSBranchManager::makeBranch( std::string name, BType type, std::string doc ){

    makeBranch( name, name, type, doc );

}//<<>>void KUCMSBranchManager::makeBranch( std::string name, BType type, std::string doc )

void KUCMSBranchManager::attachBranches( TTree* fOutTree ){ 

    for( auto & branch : theBranches ){ branch.second->attachBranch( fOutTree ); }

}//<<>>void KUCMSBranchManager::attachBranches( TTree* fOutTree )

bool KUCMSBranchManager::valid( std::string key ){ 

        if( theBranches.find(key) == theBranches.end() ){ 
            std::cout << " -- BM Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; 
            return false; 
        } else return true;

}//<<>>bool KUCMSBranchManager::valid( std::string key )

void KUCMSBranchManager::clearBranches(){ for( auto & branch : theBranches ){ (branch.second)->clear();}}
void KUCMSBranchManager::fillBranch( std::string key, std::vector<uInt> val ){ if(valid(key)) theBranches[key]->fill( val );}
void KUCMSBranchManager::fillBranch( std::string key, uInt val ){ if(valid(key)) theBranches[key]->fill( val );}
void KUCMSBranchManager::fillBranch( std::string key, int val ){ if(valid(key)) theBranches[key]->fill( val );}
void KUCMSBranchManager::fillBranch( std::string key, float val ){ if(valid(key)) theBranches[key]->fill( val );}
void KUCMSBranchManager::fillBranch( std::string key, std::string val ){ if(valid(key)) theBranches[key]->fill( val );}
void KUCMSBranchManager::fillBranch( std::string key, bool val ){ if(valid(key)) theBranches[key]->fill( val );}

#endif
