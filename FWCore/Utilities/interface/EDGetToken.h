#ifndef FWCore_Utilities_EDGetToken_h
#define FWCore_Utilities_EDGetToken_h
// -*- C++ -*-
//
// Package:     FWCore/Utilities
// Class  :     EDGetToken
// 
/**\class EDGetToken EDGetToken.h "FWCore/Utilities/interface/EDGetToken.h"

 Description: A Token used to get data from the EDM

 Usage:
    A EDGetToken is created by calls to 'consumes' or 'mayConsume' from an EDM module.
 The EDGetToken can then be used to quickly retrieve data from the edm::Event, edm::LuminosityBlock or edm::Run.
 
The templated form, EDGetTokenT<T>, is the same as EDGetToken except when used to get data the framework
 will skip checking that the type being requested matches the type specified during the 'consumes' or 'mayConsume' call.

*/
//
// Original Author:  Chris Jones
//         Created:  Wed, 03 Apr 2013 17:54:11 GMT
// $Id: EDGetToken.h,v 1.1 2013/04/14 19:01:14 chrjones Exp $
//

// system include files

// user include files

// forward declarations
namespace edm {
  class EDConsumerBase;
  template <typename T> class EDGetTokenT;
  
  class EDGetToken
  {
    friend class EDConsumerBase;
    
  public:
    static const unsigned int s_uninitializedValue = 0xFFFFFFFF;
    
    EDGetToken() : m_value{s_uninitializedValue} {}
    template<typename T>
    EDGetToken(EDGetTokenT<T> iOther): m_value{iOther.value()} {}
    
    // ---------- const member functions ---------------------
    unsigned int value() const {return m_value;}
    

  private:
    explicit EDGetToken(unsigned int iValue) : m_value{iValue} {}
    
    // ---------- member data --------------------------------
    unsigned int m_value;
  };

  template<typename T>
  class EDGetTokenT
  {
    friend class EDConsumerBase;
    
  public:
    static const unsigned int s_uninitializedValue = 0xFFFFFFFF;
    EDGetTokenT() : m_value{s_uninitializedValue} {}
    
    
    // ---------- const member functions ---------------------
    unsigned int value() const {return m_value;}
    
    
  private:
    explicit EDGetTokenT(unsigned int iValue) : m_value{iValue} {}
    
    // ---------- member data --------------------------------
    unsigned int m_value;
  };

  
}


#endif