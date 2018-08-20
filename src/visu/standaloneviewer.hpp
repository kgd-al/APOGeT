#include <QApplication>
#include <QMainWindow>

#include "kgd/external/cxxopts.hpp"
#include "kgd/external/json.hpp"
#include "kgd/utils/utils.h"

#include "kgd/settings/configfile.h"

#include "phylogenyviewer.h"

template <typename GENOME, typename CONFIG>
int run(int argc, char *argv[]) {
  using PTree = phylogeny::PhylogenicTree<GENOME>;
  using PViewer = gui::PhylogenyViewer<GENOME>;

  cxxopts::Options options("PTreeViewer", "Loads and displays a phenotypic tree for \""
                           + utils::className<GENOME>() + "\" genomes");
  options.add_options()
    ("h,help", "Display help")
    ("c,config", "File containing configuration data", cxxopts::value<std::string>())
    ("v,verbosity", "Verbosity level. " + verbosityValues(), cxxopts::value<Verbosity>())
    ("minSurvival", "Minimal survival duration", cxxopts::value<uint>())
    ("minEnveloppe", "Minimal fullness for the enveloppe", cxxopts::value<float>())
    ("showNames", "Whether or not to show node names", cxxopts::value<bool>())
    ("circular", "Whether or not to render a circular p-tree", cxxopts::value<bool>())
    ("p,print", "Render p-tree into 'filename'", cxxopts::value<std::string>())
    ("t,tree", "File containing the phenotypic tree [MANDATORY]", cxxopts::value<std::string>())
    ;

  auto result = options.parse(argc, argv);
  options.parse_positional("tree");

  if (result.count("help")) {
      std::cout << options.help() << std::endl;
      return 0;
  }

  if (result.count("tree") != 1) {
    std::cerr << "Missing mandatory argument 'tree'" << std::endl;
    return 1;
  }

  std::string configFile;
  if (result.count("config"))    configFile = result["config"].as<std::string>();

  Verbosity verbosity = Verbosity::SHOW;
  if (result.count("verbosity")) verbosity = result["verbosity"].as<Verbosity>();

  CONFIG::setupConfig(configFile, verbosity);

  auto config = PViewer::defaultConfig();

  if (result.count("minSurvival"))  config.minSurvival = result["minSurvival"].as<uint>();
  if (result.count("minEnveloppe"))  config.minEnveloppe = result["minEnveloppe"].as<float>();
  if (result.count("showNames"))  config.showNames = result["showNames"].as<bool>();
  if (result.count("circular"))  config.circular = result["circular"].as<bool>();

  QString outfile;
  if (result.count("print"))  outfile = QString::fromStdString(result["print"].as<std::string>());

  QApplication a(argc, argv);
  setlocale(LC_NUMERIC,"C");

  PTree pt = nlohmann::json::parse(utils::readAll(result["tree"].as<std::string>()));
  PViewer pv (nullptr, pt, config);
//  std::cout << pt << std::endl;

  pv.update();

  if (outfile.isEmpty()) {
    pv.show();
    return a.exec();

  } else {
    pv.printTo(outfile);
    return 0;
  }

}
