#include "MyTauClass/interface/MyTauClass.h"
#include "MyTauClass/interface/TauMatch.h"

#include <iostream>

int main(int argv, char* argc[])
{
    MyTauClass tc;
    TauMatch   tm;
    tc.MegaLoop(argc[0]);
    // need to hadd those files
    tc.show();
    tm.Loop();
    return 0;
}