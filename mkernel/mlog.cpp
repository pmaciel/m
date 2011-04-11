
#include "mlog.h"


// text format: number of spaces per indentation and time stamp
// ('clock' [0] is over 50% faster than 'date' [1])
#define MLOG_INDENT_SIZE    2
#define MLOG_TIMESTAMP_DATE 0

// error/warning messages exception throw integer value ([0] doesn't throw)
#define MLOG_MESSAGE_ERROR_THROWS 42
#define MLOG_MESSAGE_WARN_THROWS   0


namespace mlog {


  // implementation of logging levels descriptions
  const char* LEVELs[] = { (char*) "error", (char*) "warning", (char*) "info", (char*) "verbose", (char*) "debug" };
  const char  LEVELc[] =   "E"              "W"                "I"             "V"                "D";


  // implementation of facility member functions

  facility::facility() : m_level(INFO) {}
  facility::~facility() {
    while (m_indent.size())
      deindent();
    for (t_lvl_channel_it c=m_channels.begin(); c!=m_channels.end(); ++c)
      delete c->second;
  }

  void facility::indent(const std::string& name)
  {
    if (m_channels.size())
      m_indent.push_back(name.length()? name:"scope");
    for (t_lvl_channel_it c=m_channels.begin(); c!=m_channels.end(); ++c)
      c->second->indent(m_indent.back());
  }

  void facility::deindent(const float& duration)
  {
    if (m_indent.size()) {
      for (t_lvl_channel_it c=m_channels.begin(); c!=m_channels.end(); ++c)
        c->second->deindent(duration);
      m_indent.pop_back();
    }
  }

  facility& facility::log(LEVELi lvl)
  {
    instance().m_level = lvl;
    return instance();
  }

  void facility::channel(LEVELi lvl, LogChannel* chn)
  {
    instance().m_channels.push_back(t_lvl_channel(lvl,chn));
  }


  // implementation of logging interfaces according to message level
  facility& error() { return facility::log(ERROR); }
  facility& warn()  { return facility::log(WARN);  }
  facility& info()  { return facility::log(INFO);  }
  facility& verb()  { return facility::log(VERB);  }
  facility& debug() { return facility::log(DEBUG); }


  // logging manipulator: stream flush object
  struct log_manip_nl nl;


  // logging manipulator: implementation of linking stream to manipulator
  facility& operator<<(facility& l, log_manip_nl& m) { m(l); return l; }


  // logging manipulator: implementation of stream flushing (new line)
  void log_manip_nl::operator()(facility& l)
  {
    std::ostringstream& s = l.m_stream;

    // provide indentation and message
    const size_t indent(l.m_indent.size() * MLOG_INDENT_SIZE);
    const std::string msg(s.str());
    s.str("");

    // provide a time stamp
#if MLOG_TIMESTAMP_DATE
    const time_t systime_time(time(NULL));
    char systime_cptr[22];
    strftime(systime_cptr,22,"[%Y/%m/%d %H:%M:%S]",localtime(&systime_time));
    const std::string timestamp(systime_cptr);
#else
    const std::istream::fmtflags old_flags = s.setf(std::istream::right | std::istream::fixed);
    const std::streamsize        old_prec  = s.precision(2),
                                 old_width = s.width(8);

    s << float(clock())/CLOCKS_PER_SEC;
    const std::string timestamp("[" + s.str() + "s]");
    s.str("");

    s.flags    (old_flags);
    s.precision(old_prec);
    s.width    (old_width);
#endif

    // broadcast on registered channels, according to their levels
    for (facility::t_lvl_channel_it c=l.m_channels.begin(); c!=l.m_channels.end(); ++c)
      if (l.m_level <= c->first)
        c->second->message(l.m_level,indent,timestamp,msg);

    // throw (integer) exception if indicated to do so
    if (MLOG_MESSAGE_ERROR_THROWS && l.m_level==ERROR)
      throw MLOG_MESSAGE_ERROR_THROWS;
    else if (MLOG_MESSAGE_WARN_THROWS && l.m_level==WARN)
      throw MLOG_MESSAGE_WARN_THROWS;
  }


  // implementation of scope notification cons/destructor
  scope::scope(const std::string& name) : begin(clock())
  {
    facility::instance().indent(name);
  }
  scope::~scope()
  {
    facility::instance().deindent(float(clock()-begin)/CLOCKS_PER_SEC);
  }


  // channels: implementation of included stream channel member function
  void LogChannelStream::message(const LEVELi &lvl, const size_t &idt, const std::string &tim, const std::string &msg)
  {
    s << tim << ' ' << LEVELc[lvl] << ": " << std::string(idt,' ') << msg << std::endl;
  }


  // channels: implementation of included xml stream channel member functions
  LogChannelXML::LogChannelXML(std::ostream& stream) : LogChannelStream(stream)
  {
    s << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl
      << "<log>" << std::endl;
  }
  LogChannelXML::~LogChannelXML()
  {
    /*FIXME this segfaults for some reason... maybe stream is already closed?
    s << "</log>" << std::endl; */
  }
  void LogChannelXML::message(const LEVELi& lvl, const size_t& idt, const std::string& tim, const std::string& m)
  {
    s << '<' << LEVELs[lvl] << " time=\"" << tim << "\">" << m << "</" << LEVELs[lvl] << '>' << std::endl;
  }
  void LogChannelXML::indent(const std::string &name)
  {
    s <<  "<scope" << (name.length()? " name=\""+name+'"':"") << '>' << std::endl;
  }
  void LogChannelXML::deindent(const float& duration)
  {
    s << "<scope duration=\"" << duration << "s\"/>" << '\n'
      << "</scope>" << std::endl;
  }


}  // namespace mlog


#undef MLOG_INDENT_SIZE
#undef MLOG_TIMESTAMP_DATE
#undef MLOG_MESSAGE_ERROR_THROWS
#undef MLOG_MESSAGE_WARN_THROWS


