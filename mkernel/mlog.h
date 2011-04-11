#ifndef __MLOG_H__
#define __MLOG_H__


#include <ctime>
#include <sstream>
#include <vector>


namespace mlog {


/// Logging levels and their descriptions
enum               LEVELi { ERROR, WARN, INFO, VERB, DEBUG };
extern const char* LEVELs[];
extern const char  LEVELc[];


/// Logging channel base (to derive from)
struct LogChannel
{
  virtual ~LogChannel() {}
  virtual void message(
    const LEVELi& lvl,
    const size_t& idt,
    const std::string& tim,
    const std::string& msg ) = 0;
  virtual void indent  (const std::string& name="") {}
  virtual void deindent(const float& duration=0.f)  {}
};


/// Example channel to log a (fairly) plain text message to a given stream
class LogChannelStream : public LogChannel
{
 public:
  LogChannelStream(std::ostream& stream) : s(stream) {}
  virtual ~LogChannelStream() {}
  virtual void message(
    const LEVELi& lvl,
    const size_t& idt,
    const std::string& tim,
    const std::string& msg );
 protected:
  std::ostream &s;
};


/// Example channel to log an xml-encoded text message to a given stream
class LogChannelXML : public LogChannelStream
{
 public:
  LogChannelXML(std::ostream& stream);
  ~LogChannelXML();
  void message(
    const LEVELi& lvl,
    const size_t& idt,
    const std::string& tim,
    const std::string& m );
  void indent  (const std::string &name="");
  void deindent(const float& duration=0.f);
};


/// Logging stream flush manipulator
class facility;
struct log_manip_nl
{
  void operator()(facility& l);
}
extern nl;


/// Logging scope
struct scope {
  scope(const std::string& name="");
  ~scope();
  clock_t begin;
};


/// Logging facility
class facility {

  // friend function(s), to access private members
  friend void log_manip_nl::operator()(facility& l);
  friend scope::scope(const std::string& name);
  friend scope::~scope();

 private:

  /// Facility instance access (Meyer singleton)
  static facility& instance() { static facility instance; return instance; }

  // channels access (private, because nobody else should touch this)
  typedef std::pair< LEVELi, LogChannel* > t_lvl_channel;
  typedef std::vector< t_lvl_channel >::iterator t_lvl_channel_it;

  /// Facility cons/destructor
  facility();
  ~facility();

  /// Apply indentation/deindentation (and notify channels)
  void indent(const std::string& name);
  void deindent(const float& duration=0.f);

 public:

  /// Access the facility at given message level
  static facility& log(LEVELi lvl=INFO);

  /// Control channel registration
  static void channel(LEVELi lvl, LogChannel* chn);

  /// Defer operator<< to the internal stream
  template< typename T >
  facility& operator<<(const T& t) { m_stream << (t); return instance(); }

 private:

  // members
  LEVELi                       m_level;     /// Message level
  std::ostringstream           m_stream;    /// Message buffer stream
  std::vector< std::string >   m_indent;    /// Message indentation scopes
  std::vector< t_lvl_channel > m_channels;  /// Channels

};


/// Logging interfaces declarations, according to message level
facility& error();
facility& warn();
facility& info();
facility& verb();
facility& debug();


/// Logging manipulator: declare stream flushing manipulator to facility
facility& operator<<(facility& l, log_manip_nl& m);


}  // namespace mlog


#endif // __MLOG_H__
