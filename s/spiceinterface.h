#ifndef SPICEINTERFACE_H
#define SPICEINTERFACE_H


namespace {
extern "C" {
#include "ext/ngspice/src/include/ngspice/sharedspice.h"
}
}


namespace spice {


/// @brief SPICE interface_t, for managing circuit simulations contexts
struct interface_t{

  /**
   * @brief Simulation logging
   * @param msg  message
   * @param id   identification number of calling ngspice shared lib
   * @param data return pointer received from caller
   */
  static int SendChar(char* msg, int id, void *data);

  /**
   * @brief Simulation simulation status
   * @param msg  message
   * @param id   identification number of calling ngspice shared lib
   * @param data return pointer received from caller
   */
  static int SendStatus(char* msg, int id, void *data);

  /**
   * @brief Simulation controlled exit
   * @param st     exit status
   * @param unload if true: immediate unloading dll, if false: just set flag, unload is done when function has returned
   * @param exit   if true: exit upon 'quit', if false: exit due to ngspice.dll error
   * @param id     identification number of calling ngspice shared lib
   * @param data   return pointer received from caller
   */
  static int ControledExit(int st, bool unload, bool exit, int id, void *data);

  /**
   * @brief Data communication: manage simulation data
   * @param v    pointer to array of structs containing actual values from all vectors
   * @param n    number of structs (one per vector)
   * @param id   identification number of calling ngspice shared lib
   * @param data return pointer received from caller
   */
  static int SendData(pvecvaluesall v, int n, int id, void *data)
  {
    // FIXME not implemented
    return 0;
  }

  /**
   * @brief Data communication: manage simulation initialization data
   * @param v    pointer to array of structs containing data from all vectors right after initialization
   * @param id   identification number of calling ngspice shared lib
   * @param data return pointer received from caller
   */
  static int SendInitData(pvecinfoall v, int id, void *data)
  {
    // FIXME not implemented
    return 0;
  }

  /**
   * @brief Simulation thread running check
   * @param running true if background thread is running
   * @param id      identification number of calling ngspice shared lib
   * @param data    return pointer received from caller
   */
  static int BGThreadRunning(bool running, int id, void *data);

};


}  // namespace spice


#endif // SPICEINTERFACE_H
