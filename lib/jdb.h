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
  size_t m_ob_count;
  size_t m_app_count;
  size_t m_comp_count;
  size_t m_join_count;

  Eqn * m_app_data;
  Eqn * m_comp_data;
  Eqn * m_join_data;

  double m_app_prob;
  double m_comp_prob;
  double m_join_prob;
  double m_atom_prob;
public:
  typedef std::map<Ob, double> AtomProbs;
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

  size_t ob_count () const { return m_ob_count; }
  size_t app_count () const { return m_app_count; }
  size_t comp_count () const { return m_comp_count; }
  size_t join_count () const { return m_join_count; }

  const Eqn * apps () const { return m_app_data; }
  const Eqn * comps () const { return m_comp_data; }
  const Eqn * joins () const { return m_join_data; }

  double app_prob () const { return m_app_prob; }
  double comp_prob () const { return m_comp_prob; }
  double join_prob () const { return m_join_prob; }
  double atom_prob () const { return m_atom_prob; }
  const AtomProbs & atom_probs () const { return m_atom_probs; }
  double atom_prob (Ob ob)
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
