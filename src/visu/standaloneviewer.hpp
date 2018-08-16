#include <QApplication>
#include <QMainWindow>

#include "../tools/external/cxxopts.hpp"
#include "../tools/utils/utils.h"

template <typename GENOME>
int run(int argc, char *argv[]) {
  cxxopts::Options options("P-Tree viewer", "Loads and displays a phenotypic tree for \""
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
    ("t,tree", "File containing the phenotypic tree", cxxopts::value<std::string>())
    ;

  auto result = options.parse(argc, argv);
  options.parse_positional("tree");

  if (result.count("help")) {
      std::cout << options.help() << std::endl;
      return 0;
  }

  if (result.count("tree") != 1) {
    std::cerr << "Wrong number of arguments" << std::endl;
    return 1;
  }

  std::string configFile;
  if (result.count("config"))    configFile = result["config"].as<std::string>();

  Verbosity verbosity = Verbosity::SHOW;
  if (result.count("verbosity")) verbosity = result["verbosity"].as<Verbosity>();

  SimuConfig::setupConfig(configFile, verbosity);

  auto config = gui::PhylogenyViewer::defaultConfig();

  if (result.count("minSurvival"))  config.minSurvival = result["minSurvival"].as<uint>();
  if (result.count("minEnveloppe"))  config.minEnveloppe = result["minEnveloppe"].as<float>();
  if (result.count("showNames"))  config.showNames = result["showNames"].as<bool>();
  if (result.count("circular"))  config.circular = result["circular"].as<bool>();

  QString outfile;
  if (result.count("print"))  outfile = QString::fromStdString(result["print"].as<std::string>());

  QApplication a(argc, argv);
  setlocale(LC_NUMERIC,"C");

  simulation::Simulation::PTree pt = json::parse(utils::readAll(result["tree"].as<std::string>()));
  gui::PhylogenyViewer pv (nullptr, pt, config);
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
