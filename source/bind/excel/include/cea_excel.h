#ifndef CEA_EXCEL_H
#define CEA_EXCEL_H

#ifdef _WIN32
    #ifdef CEA_EXCEL_BUILD
        #define CEA_EXCEL_API __declspec(dllexport)
    #else
        #define CEA_EXCEL_API __declspec(dllimport)
    #endif
#else
    #define CEA_EXCEL_API __attribute__((visibility("default")))
#endif

#define CEA_EXCEL_SUCCESS 0
#define CEA_EXCEL_ERROR_NULL_BUFFER 1
#define CEA_EXCEL_ERROR_INVALID_LENGTH 2
#define CEA_EXCEL_ERROR_TRUNCATED 3

#ifdef __cplusplus
extern "C" {
#endif

CEA_EXCEL_API int cea_excel_version(char *version, int version_len);
CEA_EXCEL_API int cea_excel_test_add(int a, int b);

#ifdef __cplusplus
}
#endif

#endif
