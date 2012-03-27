
#include "jdb.h"
#include <iostream>
#include <utility>
#include <cstdio>
#include <errno.h>

#define LOG(message) { std::cout << message << std::endl; }
#define ERROR(message) { LOG("ERROR " << message); abort(); }
#define ASSERT(condition, message) { if (!(condition)) { ERROR(message); } }

namespace Johann
{

//----------------------------------------------------------------------------
// WARNING this must agree with Version struct in johann/cpp/version.h

struct Version
{
  uint8_t a,b,c,d;
  uint32_t num () const { return a << 24 | b << 16 | c << 8 | d; }
};
const Version OLDEST_COMPATIBLE_VERSION = {0,9,1,0};
const Version NEWEST_COMPATIBLE_VERSION = {0,9,1,255};

//----------------------------------------------------------------------------
// WARNING these structures & constants must match those in johann/cpp/files.h

struct Header
{
  Version version;

  // integer properties
  uint32_t props1, props2; // various properties

  // brain age
  uint64_t age;

  // structure sizes
  uint32_t b_size, o_size;         // number of atoms/obs
  uint32_t a_size, c_size, j_size; // number of app/comp/join eqations
  uint32_t w_size, r_size;         // number of words/rules in language
  uint32_t z_size;                 // reserved for later use

  // data offsets, in units of BLOCK_SIZE_BYTES
  uint32_t b_data; // for the atomic basis
  uint32_t o_data; // for ob properties
  uint32_t a_data; // for app eqns
  uint32_t c_data; // for comp eqns
  uint32_t j_data; // for join eqns
  uint32_t l_data; // for the [= table
  uint32_t L_data; // for language
  uint32_t y_data, z_data; // reserved for later use
};

struct ObNamePair
{
  uint32_t ob;
  char name[16];
};

typedef std::pair<uint32_t, float> IntWMass;

const size_t BLOCK_SIZE_BYTES = 256;

//----------------------------------------------------------------------------

inline FILE * safe_fopen (std::string filename)
{
  FILE * file = fopen(filename.c_str(), "rb");
  ASSERT(file, "failed to open file " << filename << " for loading");
  return file;
}

inline FILE * safe_fclose (FILE * file)
{
  int info = fclose(file);
  ASSERT(info == 0, "fclose failed with errno = " << errno);
}

inline void safe_fseek (FILE * file, size_t block_pos)
{
  int info = fseek(file, pos * BLOCK_SIZE_BYTES, SEEK_SET);
  ASSERT(info == 0, "fseek failed with errno = " << errno);
}

template<class T>
inline void safe_fread (T * dest, FILE * file, size_t count = 1)
{
  size_t info = fread(dest, sizeof(T), count, file);
  ASSERT(info == count, "fread read " << info << " of " << count << " items");
}

template<class T>
inline T * safe_malloc (size_t count)
{
  T * result = malloc(count * sizeof(T));
  ASSERT(T, "failed to allocate " << count << " items");
}

// WARNING must match load(string) in johann/src/meas_lite.C
Database::Database (std::string filename)
{
  LOG("loading database from " << filename);

  FILE * file = safe_fopen(filename);

  Header header;
  safe_fread(&Header, file);

  TODO validate header

  m_ob_size = header.o_size;
  m_app_size = header.a_size;
  m_comp_size = header.c_size;
  m_join_size = header.j_size;

  LOG(" loading " << m_app_size << " application equations");
  m_app_data = safe_malloc<Ob>(m_app_size);
  safe_fseek(file, header.a_data);
  safe_fread(m_app_data, file, m_app_size);

  LOG(" loading " << m_comp_size << " composition equations");
  m_comp_data = safe_malloc<Ob>(m_comp_size);
  safe_fseek(file, header.a_data);
  safe_fread(m_comp_data, file, m_comp_size);

  LOG(" loading " << m_join_size << " join equations");
  m_join_data = safe_malloc<Ob>(m_join_size);
  safe_fseek(file, header.a_data);
  safe_fread(m_join_data, file, m_join_size);

  // WARNING must match Language::load_from_file in johann/src/languages.C
  LOG(" loading probabilistic grammar with " << header.w_size << " atoms");
  safe_fseek(file, header.L_data);
  safe_fread(& m_app_prob, file);
  safe_fread(& m_comp_prob, file);
  safe_fread(& m_join_prob, file);
  {
    float atom_prob = 1.0f - m_app_prob - m_comp_prob - m_join_prob;
    IntWMass * pairs = safe_malloc<IntWMass>(header.w_size);
    safe_fread(pairs, file, header.w_size);
    for (size_t i = 0; i < header.w_size; ++i) {
      IntWMass pair = pairs[i];
      m_atom_probs[pair.first] = atom_prob * pair.second;
    }
    free(pairs);
  }

  // WARNING must match Brain::load(filename) in johann/src/brain.C
  LOG(" loading names of " << header.b_size << " obs");
  safe_fseek(file, header.b_data);
  // WARNING must match CombinatoryStructure::load_from_file(filename,b_size)
  {
    ObNamePair * pairs = safe_malloc<ObNamePair>(header.b_size);
    safe_fread(pairs, file, header.b_size);
    for (size_t i = 0; i < header.b_size; ++i) {
      m_name_to_ob[pairs[i].name] = pairs[i].ob;
    }
    free(pairs);
  }
  ASSERT(m_name_to_ob.size() == header.b_size,
      "duplicated names in loading atoms");

  safe_fclose(file);
}

Database::~Database ()
{
  free(m_app_data);
  free(m_comp_data);
  free(m_join_data);
}

Database (const Database &)
{
  ERROR("Database copying is not supported");
}

void Databse::operator= (const Database &)
{
  ERROR("Database copying is not supported");
}

} // namespace Johann

