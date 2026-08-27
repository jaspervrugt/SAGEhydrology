#pragma once

#include <io.h>
#include <cstdio>
#include <string>
#include <mutex>
#include <initializer_list>

#if defined(_WIN32)
// #define NOMINMAX
#include <windows.h>
#else
#include <sys/file.h>
#include <unistd.h>
#endif

// ----------------------------------------------------------------------------
// Header-only fail-parameter logger for MEX models
// Usage:
//   #include "log_fail_params.hpp"
//   log_fail_params("gr4j", x1, x2, x3, x4, tt, ddf);
// ----------------------------------------------------------------------------

namespace fail_log {

inline std::mutex& process_mutex()
{
    static std::mutex m;
    return m;
}

// Cross-process advisory lock helpers
inline bool lock_file(FILE* fp)
{
    if (!fp) {
        return false;
    }

#if defined(_WIN32)
    // Lock the whole file (best-effort)
    HANDLE h = (HANDLE)_get_osfhandle(_fileno(fp));
    if (h == INVALID_HANDLE_VALUE) {
        return false;
    }

    OVERLAPPED ov = {};
    // Lock 0..MAXDWORD (whole file)
    return LockFileEx(h, LOCKFILE_EXCLUSIVE_LOCK, 0, MAXDWORD, MAXDWORD, &ov) != 0;
#else
    int fd = fileno(fp);
    if (fd < 0) {
        return false;
    }
    return flock(fd, LOCK_EX) == 0;
#endif
}

inline void unlock_file(FILE* fp)
{
    if (!fp) {
        return;
    }

#if defined(_WIN32)
    HANDLE h = (HANDLE)_get_osfhandle(_fileno(fp));
    if (h == INVALID_HANDLE_VALUE) {
        return;
    }

    OVERLAPPED ov = {};
    UnlockFileEx(h, 0, MAXDWORD, MAXDWORD, &ov);
#else
    int fd = fileno(fp);
    if (fd >= 0) {
        flock(fd, LOCK_UN);
    }
#endif
}

inline std::string default_path()
{
    // If you want per-model location, change this, or pass a path into log_fail_params().
    return "fail_params.txt";
}

} // namespace fail_log

// Variadic template: supports different models with different #parameters
template <typename... Ts> inline void log_fail_params(const char* model_name, Ts... vals)
{
    // 1) Within-process safety
    std::lock_guard<std::mutex> guard(fail_log::process_mutex());

    const std::string path = fail_log::default_path();
    FILE* fp = std::fopen(path.c_str(), "a");
    if (!fp) {
        return;
    }

    // 2) Cross-process lock (best-effort)
    (void)fail_log::lock_file(fp);

    // Write: model_name followed by values (tab-separated), one row
    std::fprintf(fp, "%s", (model_name ? model_name : "UNKNOWN"));

    // Print values with good precision for debugging
    // (casts to double to keep formatting consistent for numeric inputs)
    const double arr[] = {static_cast<double>(vals)...};
    for (double v : arr) {
        std::fprintf(fp, "\t%.17g", v);
    }
    std::fprintf(fp, "\n");
    std::fflush(fp);

    fail_log::unlock_file(fp);
    std::fclose(fp);
}
