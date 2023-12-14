#include "HepMC3/ReaderFactory.h"
#include "HepMC3/ReaderAsciiHepMC2.h"
int main (int argc, char** argv) {
    if (argc!=2) return 1;
    auto inputA = HepMC3::deduce_reader(argv[1]);
    if (!inputA||inputA->failed()) return 1;
    size_t events=0;
    while( !inputA->failed() ) {
        HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
        bool res = inputA->read_event(evt);
        if( inputA->failed() )  {
            printf("End of file reached. Exit.\n");
            break;
        }
        events++;
        if (evt.particles().size() < 2) {
            printf("Too few particles\n");
            return 2;
        }
        evt.clear();
    }
    inputA->close();
    if (events == 0) {
        printf("Too few events\n");
        return 3;
    }
    return 0;
}
