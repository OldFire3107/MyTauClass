#include "MyTauClass/TauClass/interface/MyTauClass.h"
#include "MyTauClass/TauClass/interface/TauMatch.h"

#include <iostream>

int main(int argc, char** argv)
{
    if(argc > 2){
        MyTauClass tc(argv[1], argv[2]);
        tc.Loop();
    }
    return 0;
}
