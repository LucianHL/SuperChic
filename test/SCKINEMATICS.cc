#include "Rivet/Analysis.hh"
const size_t MAXPARTICLES  = 15;
namespace Rivet {
class SCKINEMATICS : public Analysis {
public:

    RIVET_DEFAULT_ANALYSIS_CTOR(SCKINEMATICS);
    void init() {
        book(_h["size"], "size", 40, 0.0, 20.0);
        for (size_t i = 0; i < MAXPARTICLES+1; i++) {
          book(_h["pt" + std::to_string(100+i+1)], "pt" + std::to_string(100+i+1), 50, 0.0, 10.0);
          book(_h["eta" + std::to_string(100+i+1)], "eta" + std::to_string(100+i+1), 16, -8.0, 8.0);
        }
    }
    void analyze(const Event& event) {
        auto e = event.hepmcEventPtr();
        const auto& particles = e->particles();
        _h["size"]->fill(1.0*particles.size());
        for (size_t i = 0; i < std::min(MAXPARTICLES, particles.size()); i++) {
          const auto& mom  = particles.at(i)->momentum();
          _h["pt" + std::to_string(100+i+1)]->fill(mom.pt());
          _h["eta" + std::to_string(100+i+1)]->fill(mom.eta());
        }
    }
    void finalize() {
//        const double sf = crossSection()/picobarn/sumOfWeights();
        const double sf = 1.0/sumOfWeights();
        scale(_h,sf);
    }
    map<string, Histo1DPtr> _h;
};
RIVET_DECLARE_PLUGIN(SCKINEMATICS);
}
