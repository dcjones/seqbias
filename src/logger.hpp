
#ifndef SEQBIAS_LOGGER_HPP
#define SEQBIAS_LOGGER_HPP

#include <cstdarg>
#include <cstdlib>

/** This is meant to be striped-down, more R-friendly mimic of the logger class
 * in isolator.
 */
class logger
{
    public:
        static logger& instance();

        static void debug (const char* fmt, ...);
        static void info  (const char* fmt, ...);
        static void warn  (const char* fmt, ...);
        static void error (const char* fmt, ...);

        enum level 
        {
            DEBUG,
            INFO,
            WARN,
            ERROR
        };

        static void set_level(level);

        /** print a message and exit. */
        static void abort (const char* fmt, ...);

        static void print(const char* msg);

        static void start();
        static void suspend();
        static void resume();

    private:
        logger();
        ~logger();

        /** Print a log message at the given level. */
        void put(level, const char* fmt, va_list args); 


        level L;
        char* msg_buffer;
        static const size_t msg_buffer_len;
};


#endif


