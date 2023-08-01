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
#include <iostream>
#include <limits> 

//Root includes
#include "TTree.h"

#ifndef KUCMSItemHeader
#define KUCMSItemHeader


//.............................................................................................

typedef unsigned int uInt;

template <class T> 
class Item {

    public:

    Item();
    Item( std::string name, std::string doc = "" );
    Item( std::string name, T value, std::string doc = "" );
    //~Item();

    void make( std::string name, std::string doc = "" );
    void fill( T val );
    void clear();
    T getvalue();
    std::string getdoc();

    //protected:

    std::string iName;
    std::string iDoc;
    T iValue;

};//<<>>class Item 

template <class T>
Item<T>::Item() : iName(""), iDoc(""), iValue(T()) {}

template <class T>
Item<T>::Item( std::string name, std::string doc ):

    iName(name),
    iDoc(doc)

{}//<<>>Item::Item( std::string name, std::string doc )

template <class T>
Item<T>::Item( std::string name, T value, std::string doc ):
    
    iName(name),
    iDoc(doc),
    iValue(value)

{}//<<>>Item::Item( std::string name, std::string doc, T value  )

template <class T>
std::string Item<T>::getdoc(){

    return iDoc;

}//<<>>std::string Item::getDoc

template <class T>
T Item<T>::getvalue(){

    return iValue;

}//<<>>Item::getVal<T>()

template <class T>
void Item<T>::clear(){

    iValue = -1.0*std::numeric_limits<T>::max();

}//<<>>void Item::clear()

template <>
void Item<std::string>::clear(){

    iValue = "";

}//<<>>void Item::clear()

template <>
void Item<bool>::clear(){

    iValue = false;

}//<<>>void Item::clear()

template <>
void Item<uInt>::clear(){

    iValue = 0;

}//<<>>void Item::clear()

template <>
void Item<std::vector<unsigned int>>::clear(){

    iValue.clear();

}//<<>>void Item<std::vector<unsigned int>>::clear()

template <class T>
void Item<T>::make( std::string name, std::string doc ){

    iName = name;
    iDoc = doc;

}//<<>>void Item::make(  std::string name, std::string doc  )

template <class T>
void Item<T>::fill( T val ){

    iValue = val;

}//<<>>void Item::fill( T val )

//..................................................................................
template < class T>
class VectorItem {

    public:

    VectorItem();
    VectorItem( std::string name, std::string doc );
    VectorItem( std::string name, T value, std::string doc );
    //~VectorItem();

    void make( std::string name, std::string doc = "" );
    void fill( T val );
    void clear();
    std::vector<T> getvalue();
    std::string getdoc();

    protected:

    std::string iName;
    std::string iDoc;
    std::vector<T> iVector;

};//<<>>class VectorItem

template < class T >
VectorItem<T>::VectorItem() : iName(""), iDoc("") {}

template < class T >
VectorItem<T>::VectorItem( std::string name, std::string doc ):

    iName(name),
    iDoc(doc)

{}//<<>>VectorItem::VectorItem( std::string name, std::string doc )

template < class T >
VectorItem<T>::VectorItem( std::string name, T value, std::string doc ):
 
    iName(name),
    iDoc(doc)

{ this.fill(value); }//<<>>VectorItem::VectorItem( std::string name, std::string doc, T value  )

template < class T>
std::string VectorItem<T>::getdoc(){

    return iDoc;

}//<<>>std::string VectorItem::getDoc

template < class T >
std::vector<T> VectorItem<T>::getvalue(){

    return iVector;

}//<<>>VectorItem::getVal<T>()

template < class T >
void VectorItem<T>::clear(){

    iVector.clear();

}//<<>>void VectorItem::clear()

template <class T>
void VectorItem<T>::make( std::string name, std::string doc ){

    iName = name;
    iDoc = doc;

}//<<>>void VectorItem::make(  std::string name, std::string doc )

template < class T >
void VectorItem<T>::fill( T val ){

    iVector.push_back( val );

}//<<>>void VectorItem::clear()

//......................................................................................
//   ItemManager
//......................................................................................

template < class T , template< class > class C = Item >
class ItemManager {

    public:
    
    void set( std::string key, std::string name, std::string doc );
    void set( std::string name );
    void set( std::string name, T val );

    void fill( std::string key, T val );
    T get( std::string key );

    void clear();
    void clear( std::string key );
    void reset();

    T operator()( std::string key );

    private:

    std::map< std::string, C<T> > items;
    bool valid( std::string key );

};//<<>>class ItemManager 

template < class T , template< class > class C >
void ItemManager<T,C>::set( std::string key, std::string name, std::string doc ){

    C<T> newitem( name, doc );
    items[key] = newitem;

}//>><<void ItemManager::set( std::string key, std::string name, std::string doc )

template < class T , template< class > class C >
void ItemManager<T,C>::set( std::string name ){ set( name, name, "" ); }

template < class T , template< class > class C >
void ItemManager<T,C>::set( std::string name, T val ){

    C<T> newitem( name, val, "" );
    items[name] = newitem;

}//>><<void ItemManager::set( std::string key, std::string name, std::string doc, T val )

template < class T , template< class > class C >
void ItemManager<T,C>::fill( std::string key, T val ){

    if( valid( key ) ) items[key].fill( val );
    
}//<<>>void ItemManager<T>::fill( std::string key, T val )

template < class T , template< class > class C >
T ItemManager<T,C>::get( std::string key ){

    return valid( key ) ? items[key].getvalue() : T();

}//<<>>T ItemManager::get( std::string key )

template < class T , template< class > class C >
void ItemManager<T,C>::clear(){

    for( auto & item : items ){ item.second.clear(); }

}//<<>>T ItemManager::clear()

template < class T , template< class > class C >
void ItemManager<T,C>::reset(){

    for( auto & item : items ){ item.erase(); }

}//<<>>T ItemManager::reset()

template < class T , template< class > class C >
void ItemManager<T,C>::clear( std::string key ){
    
    if( valid( key ) ) items[key].clear();

}//<<>>T ItemManager::clear( std::string key )

template < class T , template< class > class C >
T ItemManager<T,C>::operator()( std::string key ){

    return valid( key ) ? items[key].getvalue() : T();

}//<<>>T ItemManager<T,C>::operator()( std::string key )

template < class T , template< class > class C >
bool ItemManager<T,C>::valid( std::string key ){ 

    if( items.find(key) == items.end() ){
        std::cout << " -- IM Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl; 
        return false; 
    } else return true;

}//<<>>bool ItemManager<T>::::valid( std::string key )

#endif
