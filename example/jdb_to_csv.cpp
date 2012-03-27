
#include "jdb.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
//#include <boost/program_options.hpp>

using std::string;

#define LOG(message) { std::cout << message << std::endl; }

inline void write_equations (
    string outstem,
    string name,
    const Johann::Eqn * eqns,
    const size_t count)
{
  string outfile = outstem + "." + name + "s.csv";
  LOG("writing " << count << " " << name << " equations to " << outfile);
  std::ofstream file(outfile.c_str());

  file << "lhs,rhs," << name;
  for (const Johann::Eqn * i = eqns, * end = i + count; i != end; ++i) {
    file << '\n' << i->lhs << ',' << i->rhs << ',' << i->result;
  }
}

void jdb_to_csv (string infile, string outstem)
{
  Johann::Database db(infile);

  {
    string outfile = outstem + ".params.csv";
    LOG("writing parameters to " << outfile);
    std::ofstream file(outfile.c_str());

    file << "parameter,value";
    file << "\nob_count," << db.ob_count();
    file << "\napp_prob," << db.app_prob();
    file << "\ncomp_prob," << db.comp_prob();
    file << "\njoin_prob," << db.join_prob();
  }

  write_equations(outstem, "app", db.apps(), db.app_count());
  write_equations(outstem, "comp", db.comps(), db.comp_count());
  write_equations(outstem, "join", db.joins(), db.join_count());

  {
    const Johann::Database::AtomProbs & atom_probs = db.atom_probs();

    string outfile = outstem + ".weights.csv";
    LOG("writing " << atom_probs.size() << " atom weights to " << outfile);
    std::ofstream file(outfile.c_str());

    file << "ob,probability";
    typedef Johann::Database::AtomProbs::const_iterator Auto;
    for (Auto i = atom_probs.begin(); i != atom_probs.end(); ++i) {
      file << '\n' << i->first << ',' << i->second;
    }
  }

  {
    const Johann::Database::NameToOb & name_to_ob = db.name_to_ob();

    string outfile = outstem + ".names.csv";
    LOG("writing " << name_to_ob.size() << " atom names to " << outfile);
    std::ofstream file(outfile.c_str());

    file << "ob,name";
    typedef Johann::Database::NameToOb::const_iterator Auto;
    for (Auto i = name_to_ob.begin(); i != name_to_ob.end(); ++i) {
      file << '\n' << i->second << ',' << i->first;
    }
  }

  LOG("TODO write remaining structure");
}

int main (int argc, char ** argv)
{
  string infile, outstem;

  if (argc == 2) {

    infile = argv[1];
    size_t dotpos = infile.rfind(".");
    outstem = dotpos == infile.npos ? infile : infile.substr(dotpos + 1);

  } else if (argc == 3) {

    infile = argv[1];
    outstem = argv[2];

  } else {

    LOG("Usage: jdb_to_csv INFILE.jdb [OUTSTEM]");
    exit(1);
  }

  jdb_to_csv(infile, outstem);

  return 0;
}

