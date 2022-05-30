#include "MyTauClass/TauClass/interface/MyTauClass.h"
#include "MyTauClass/TauClass/interface/TauMatch.h"

#include <iostream>

int main(int argc, char** argv)
{
    if(argc > 1){
    MyTauClass tc;
    TauMatch   tm(argv[1]);
        tc.show(argv[1]);
        tm.Loop();
    }
    return 0;
}