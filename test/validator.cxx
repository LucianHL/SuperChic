#include "HepMC3/ReaderFactory.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
int main (int argc, char** argv) {
    if (argc!=2) return 1;
    auto inputA = HepMC3::deduce_reader(argv[1]);
    if (!inputA||inputA->failed()) return 1;

    while( !inputA->failed() )
    {
        HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
        bool res = inputA->read_event(evt);
        if (!res) return 2;
        if( inputA->failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }

        if (evt.particles().size()==0)return 2;
        evt.clear();
    }
    inputA->close();
    return 0;
}