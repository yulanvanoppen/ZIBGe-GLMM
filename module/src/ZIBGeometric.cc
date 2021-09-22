#include <module/Module.h>      // include JAGS module base class
#include <distributions/DZIBG.h>// include DZIBG distribution class
#include <distributions/DBEEG.h>// include DBEEG distribution class

namespace jags::ZIBGeometric {  // start defining the module namespace

                                // module class
    class ZIBGModule : public Module {  
        public:
            ZIBGModule();       // constructor
            ~ZIBGModule();      // destructor
    };
    
                                // constructor implementation
    ZIBGModule::ZIBGModule() : Module("ZIBGeometric") {
        insert(new DZIBG);      // inherited function to load objects into JAGS
        insert(new DBEEG);
    }

    ZIBGModule::~ZIBGModule() { // destructor implementation
        std::vector<Distribution*> const &dvec = distributions();
        for (unsigned int i = 0; i < dvec.size(); ++i)
            delete dvec[i];     // delete all instantiated distribution objects
    }

}                               // end namespace definition

jags::ZIBGeometric::ZIBGModule _ZIBGeometric_module;
