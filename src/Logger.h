#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>

#if DEBUG_ENABLE
#define DEBUG(msg) msg
#else
#define DEBUG(msg)
#endif

#endif /* LOGGER_H */
