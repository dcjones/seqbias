
#include "logger.hpp"
#include <R.h>


/* take care of R namespace pollution */
#if defined(error)
#undef error
#endif

#if defined(ERROR)
#undef ERROR 
#endif

#if defined(WARN)
#undef WARN
#endif



extern "C" {
    void logger_abort(const char* fmt, ...)
    {
        va_list args;
        va_start(args, fmt);

        logger::abort(fmt, args);

        va_end(args);
    }
}


const size_t logger::msg_buffer_len = 4096;

logger::logger()
{
    msg_buffer = new char [msg_buffer_len];
}

logger::~logger()
{
    delete [] msg_buffer;
}


logger& logger::instance()
{
    static logger S;
    return S;
}

void logger::set_level(level l)
{
    logger::instance().L = l;
}


void logger::put(level l, const char* fmt, va_list args)
{
    if (L > l) return;

    int len = vsnprintf(msg_buffer, msg_buffer_len, fmt, args);

    /* make sure there is exactly one trailing newlines */
    while (len > 0 && msg_buffer[len - 1] == '\n') {
        msg_buffer[len - 1] = '\0';
        --len;
    }
    msg_buffer[len] = '\n';
    msg_buffer[len + 1] = '\0';

    /* special case for warnings */
    if (l == WARN) {
        Rf_warning("%s", msg_buffer);
    }
    else {
        Rprintf("%s", msg_buffer);
    }
}


void logger::debug (const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    logger::instance().put(DEBUG, fmt, args);

    va_end(args);
}

void logger::info  (const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    logger::instance().put(INFO, fmt, args);

    va_end(args);
}

void logger::warn  (const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    logger::instance().put(WARN, fmt, args);

    va_end(args);
}

void logger::error (const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    logger::instance().put(ERROR, fmt, args);

    va_end(args);
}

void logger::print(const char* msg)
{
    Rprintf("%s", msg);
}


void logger::abort(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    vsnprintf(logger::instance().msg_buffer, msg_buffer_len, fmt, args);
    va_end(args);

    Rf_error("%s", logger::instance().msg_buffer);
}


void logger::start()
{
    // this space intentionally left blank
}

void logger::suspend()
{
    // this space intentionally left blank
}

void logger::resume()
{
    // this space intentionally left blank
}



