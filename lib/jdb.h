/** jdb - a lightweight reader for Johann databases.
 * see http://github.com/fritzo/Johann
 * see http://github.com/fritzo/johann_lite
 */

#ifndef JOHANN_JDB_H
#define JOHANN_JDB_H

#include <stdint.h>
#include <cstdlib>
#include <string>
#include <map>

namespace Johann
{

typedef uint16_t Ob;
struct Eqn { Ob lhs, rhs, result; };

class Database
{
  size_t m_ob_size;
  size_t m_app_size;
  size_t m_comp_size;
  size_t m_join_size;

  Eqn * m_app_data;
  Eqn * m_comp_data;
  Eqn * m_join_data;

  float m_app_prob;
  float m_comp_prob;
  float m_join_prob;
public:
  typedef std::map<Ob, float> AtomProbs;
private:
  AtomProbs m_atom_probs;

public:
  typedef std::map<std::string, Ob> NameToOb;
private:
  NameToOb m_name_to_ob;

  // TODO load [= and [!= order tables

  // copying is not supported
  Database (const Database &);
  void operator= (const Database &);

public:

  Database (std::string filename);
  virtual ~Database ();

  size_t ob_size () const { return m_ob_size; }
  size_t app_size () const { return m_app_size; }
  size_t comp_size () const { return m_comp_size; }
  size_t join_size () const { return m_join_size; }

  const Eqn * apps () const { return m_app_data; }
  const Eqn * comps () const { return m_comp_data; }
  const Eqn * joins () const { return m_join_data; }

  float app_prob () const { return m_app_prob; }
  float comp_prob () const { return m_comp_prob; }
  float join_prob () const { return m_join_prob; }
  const AtomProbs & atom_probs () const { return m_atom_probs; }
  float atom_prob (Ob ob)
  {
    AtomProbs::iterator i = m_atom_probs.find(ob);
    return i == m_atom_probs.end() ? 0 : i->second;
  }

  const NameToOb & name_to_ob () const { return m_name_to_ob; }
  Ob name_to_ob (const char * name)
  {
    NameToOb::iterator i = m_name_to_ob.find(name);
    return i == m_name_to_ob.end() ? 0 : i->second;
  }
};

} // namespace Johann

#endif // JOHANN_JDB_H
