#include "../core/crossover.hpp"
#include "../visu/standaloneviewer.hpp"

struct Genome {
  struct Alignment {};
  using CData = genotype::BailOutCrossover<Genome, Alignment>::Data;
  friend void to_json (nlohmann::json&, const Genome&) {}
  friend void from_json (const nlohmann::json&, Genome&) {}
};

#define CFILE BasicConfig
struct BasicConfig : public ConfigFile<CFILE> {
  DECLARE_SUBCONFIG(CrossoverConfig, crossoverConfig)
  DECLARE_SUBCONFIG(PTreeConfig, phenotypicTreeConfig)

};
DEFINE_SUBCONFIG(CrossoverConfig, crossoverConfig)
DEFINE_SUBCONFIG(PTreeConfig, phenotypicTreeConfig)
#undef CFILE

int main(int argc, char *argv[]) {
  return run<Genome, BasicConfig>(argc, argv);
}
