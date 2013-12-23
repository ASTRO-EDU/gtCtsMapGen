#ifndef __AstroI_h__
#define __AstroI_h__

#include <Astro.h>

namespace Astro
{

class AstroInterfaceI : virtual public AstroInterface
{
public:

    virtual ::Astro::Matrix calculateMapKey(const ::Astro::SeqEvtKey&,
                                            const Ice::Current&);

    virtual ::Astro::Matrix calculateMapVector(const ::Astro::Ra&,
                                               const ::Astro::Dec&,
                                               const Ice::Current&);

    virtual void shutdown(const Ice::Current&);
};

}

#endif
