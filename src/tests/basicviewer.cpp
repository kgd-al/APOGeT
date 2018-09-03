#include "../core/crossover.hpp"
#include "../visu/standaloneviewer.hpp"

struct Genome {
  struct Alignment {};
  using CData = genotype::BailOutCrossover<Genome, Alignment>::Data;
  friend void to_json (nlohmann::json&, const Genome&) {}
  friend void from_json (const nlohmann::json&, Genome&) {}
};

namespace config {
#define CFILE Basic
struct Basic : public ConfigFile<CFILE> {
  DECLARE_SUBCONFIG(Crossover, crossoverConfig)
  DECLARE_SUBCONFIG(PTree, phenotypicTreeConfig)
};
DEFINE_SUBCONFIG(Crossover, crossoverConfig)
DEFINE_SUBCONFIG(PTree, phenotypicTreeConfig)
#undef CFILE
} // end of namespace config

int main(int argc, char *argv[]) {
  return run<Genome, config::Basic>(argc, argv);
}
