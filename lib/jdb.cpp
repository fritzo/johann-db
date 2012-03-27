
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
  friend std::ostream & operator<< (std::ostream & o, const Version & v)
  {
    return o << v.a << '.' << v.b << '.' << v.c << '.' << v.d;
  }
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

typedef std::pair<uint32_t, double> IntWMass;

const size_t BLOCK_SIZE_BYTES = 256;

//----------------------------------------------------------------------------

inline FILE * safe_fopen (std::string filename)
{
  FILE * file = fopen(filename.c_str(), "rb");
  ASSERT(file, "failed to open file " << filename << " for loading");
  return file;
}

inline void safe_fclose (FILE * file)
{
  int info = fclose(file);
  ASSERT(info == 0, "fclose failed with errno = " << errno);
}

inline void safe_fseek (FILE * file, size_t block_pos)
{
  int info = fseek(file, block_pos * BLOCK_SIZE_BYTES, SEEK_SET);
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
  T * result = static_cast<T *>(malloc(count * sizeof(T)));
  ASSERT(result, "failed to allocate " << count << " items");
  return result;
}

inline void validate_equations (
    const Eqn * eqns,
    const size_t eqn_count,
    const size_t ob_count)
{
  for (const Eqn * i = eqns, * end = i + eqn_count; i != end; ++i) {
    ASSERT(0 < i->lhs and i->lhs <= ob_count, " bad lhs: " << i->lhs);
    ASSERT(0 < i->rhs and i->rhs <= ob_count, " bad rhs: " << i->rhs);
    ASSERT(0 < i->result and i->result <= ob_count,
        " bad result: " << i->result);
  }
}

// WARNING must match load(string) in johann/src/meas_lite.C
Database::Database (std::string filename)
{
  LOG("loading database from " << filename);

  FILE * file = safe_fopen(filename);

  Header header;
  safe_fread(& header, file);

  ASSERT(OLDEST_COMPATIBLE_VERSION.num() <= header.version.num(),
      "jdb file is unreadably old: version = " << header.version);
  ASSERT(header.version.num() <= NEWEST_COMPATIBLE_VERSION.num(),
      "jdb file is unreadably new: version = " << header.version);

  m_ob_count = header.o_size;
  m_app_count = header.a_size;
  m_comp_count = header.c_size;
  m_join_count = header.j_size;
  const size_t weight_count = header.w_size;
  const size_t name_count = header.b_size;
  ASSERT(m_app_count <= m_ob_count * m_ob_count,
      "bad app_count: " << m_app_count);
  ASSERT(m_comp_count <= m_ob_count * m_ob_count,
      "bad comp_count: " << m_comp_count);
  ASSERT(m_join_count <= m_ob_count * (m_ob_count + 1) / 2,
      "bad join_count: " << m_join_count);
  ASSERT(weight_count <= m_ob_count, "bad weight count: " << weight_count);
  ASSERT(name_count <= m_ob_count, "bad name count: " << name_count);

  LOG(" reading " << m_app_count << " application equations");
  m_app_data = safe_malloc<Eqn>(m_app_count);
  safe_fseek(file, header.a_data);
  safe_fread(m_app_data, file, m_app_count);
  validate_equations(m_app_data, m_app_count, m_ob_count);

  LOG(" reading " << m_comp_count << " composition equations");
  m_comp_data = safe_malloc<Eqn>(m_comp_count);
  safe_fseek(file, header.c_data);
  safe_fread(m_comp_data, file, m_comp_count);
  validate_equations(m_comp_data, m_comp_count, m_ob_count);

  LOG(" reading " << m_join_count << " join equations");
  m_join_data = safe_malloc<Eqn>(m_join_count);
  safe_fseek(file, header.j_data);
  safe_fread(m_join_data, file, m_join_count);
  validate_equations(m_join_data, m_join_count, m_ob_count);

  // WARNING must match Language::load_from_file in johann/src/languages.C
  LOG(" reading probabilistic grammar with " << weight_count << " atoms");
  safe_fseek(file, header.L_data);
  safe_fread(& m_app_prob, file);
  safe_fread(& m_comp_prob, file);
  safe_fread(& m_join_prob, file);
  m_atom_prob = 1.0 - m_app_prob - m_comp_prob - m_join_prob;
  ASSERT(0 < m_app_prob and m_app_prob < 1, "bad app prob: " << m_app_prob);
  ASSERT(0 < m_comp_prob and m_comp_prob < 1, "bad comp prob: " << m_comp_prob);
  ASSERT(0 < m_join_prob and m_join_prob < 1, "bad join prob: " << m_join_prob);
  ASSERT(0 < m_atom_prob and m_atom_prob < 1, "bad atom prob: " << m_atom_prob);
  {
    IntWMass * pairs = safe_malloc<IntWMass>(weight_count);
    safe_fread(pairs, file, weight_count);
    for (size_t i = 0; i < weight_count; ++i) {
      IntWMass pair = pairs[i];
      Ob ob = pair.first;
      double prob = pair.second;
      m_atom_probs[ob] = prob;
      ASSERT(0 < prob and prob <= 1, "bad prob: " << prob);
      ASSERT(0 < ob and ob <= m_ob_count, "bad ob: " << ob);
    }
    free(pairs);
  }

  // WARNING must match Brain::load(filename) in johann/src/brain.C
  LOG(" reading names of " << name_count << " obs");
  safe_fseek(file, header.b_data);
  // WARNING must match CombinatoryStructure::load_from_file(filename,b_size)
  {
    ObNamePair * pairs = safe_malloc<ObNamePair>(name_count);
    safe_fread(pairs, file, name_count);
    for (size_t i = 0; i < name_count; ++i) {
      Ob ob = pairs[i].ob;
      m_name_to_ob[pairs[i].name] = ob;
      ASSERT(0 < ob and ob <= m_ob_count, "ob out of range: " << ob);
    }
    free(pairs);
  }
  ASSERT(m_name_to_ob.size() == name_count, "duplicated names");

  safe_fclose(file);

  LOG(" done loading");
}

Database::~Database ()
{
  free(m_app_data);
  free(m_comp_data);
  free(m_join_data);
}

Database::Database (const Database &)
{
  ERROR("Database copying is not supported");
}

void Database::operator= (const Database &)
{
  ERROR("Database copying is not supported");
}

} // namespace Johann

