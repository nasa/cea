#include <stdio.h>
#include <string.h>

#include "cea_excel.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

typedef int (*cea_excel_version_fn)(char *version, int version_len);
typedef int (*cea_excel_test_add_fn)(int a, int b);

static int expect_true(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "%s\n", message);
        return 1;
    }
    return 0;
}

#ifdef _WIN32
static void *load_library_handle(void)
{
    HMODULE handle;

    handle = LoadLibraryA(CEA_EXCEL_TEST_LIBRARY);
    if (handle == NULL) {
        fprintf(stderr, "LoadLibrary failed for %s\n", CEA_EXCEL_TEST_LIBRARY);
        return NULL;
    }

    return (void *)handle;
}
#else
static void *load_library_handle(void)
{
    void *handle;

    handle = dlopen(CEA_EXCEL_TEST_LIBRARY, RTLD_NOW);
    if (handle == NULL) {
        fprintf(stderr, "dlopen failed for %s: %s\n",
            CEA_EXCEL_TEST_LIBRARY, dlerror());
        return NULL;
    }

    return handle;
}
#endif

int main(void)
{
    char version[256];
    char tiny[2];
    void *handle;
    cea_excel_test_add_fn cea_excel_test_add_ptr;
    cea_excel_version_fn cea_excel_version_ptr;
    int status;

    handle = load_library_handle();
    if (expect_true(handle != NULL, "Expected to load test library")) {
        return 1;
    }

#ifdef _WIN32
    cea_excel_test_add_ptr = (cea_excel_test_add_fn)GetProcAddress((HMODULE)handle, "cea_excel_test_add");
    cea_excel_version_ptr = (cea_excel_version_fn)GetProcAddress((HMODULE)handle, "cea_excel_version");
#else
    cea_excel_test_add_ptr = (cea_excel_test_add_fn)dlsym(handle, "cea_excel_test_add");
    cea_excel_version_ptr = (cea_excel_version_fn)dlsym(handle, "cea_excel_version");
#endif

    if (expect_true(cea_excel_test_add_ptr != NULL, "Expected to load cea_excel_test_add")) {
        return 1;
    }
    if (expect_true(cea_excel_test_add_ptr(2, 3) == 5, "Expected cea_excel_test_add(2,3) == 5")) {
        return 1;
    }

    if (expect_true(cea_excel_version_ptr != NULL, "Expected to load cea_excel_version")) {
        return 1;
    }

    memset(version, 0, sizeof(version));
    status = cea_excel_version_ptr(version, (int)sizeof(version));
    if (expect_true(status == CEA_EXCEL_SUCCESS, "Expected success for full buffer")) {
        return 1;
    }
    if (expect_true(version[0] != '\0', "Expected non-empty version string")) {
        return 1;
    }

    status = cea_excel_version_ptr(NULL, (int)sizeof(version));
    if (expect_true(status == CEA_EXCEL_ERROR_NULL_BUFFER, "Expected null-buffer error")) {
        return 1;
    }

    status = cea_excel_version_ptr(version, 0);
    if (expect_true(status == CEA_EXCEL_ERROR_INVALID_LENGTH, "Expected invalid-length error")) {
        return 1;
    }

    tiny[0] = 'X';
    tiny[1] = 'Y';
    status = cea_excel_version_ptr(tiny, (int)sizeof(tiny));
    if (expect_true(status == CEA_EXCEL_ERROR_TRUNCATED, "Expected truncation error")) {
        return 1;
    }
    if (expect_true(tiny[1] == '\0', "Expected tiny buffer to be null terminated")) {
        return 1;
    }

    return 0;
}
