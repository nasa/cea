#include "cea_fortran_api.h"

#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
#define CEA_THREAD_LOCAL __declspec(thread)
#else
#define CEA_THREAD_LOCAL _Thread_local
#endif

#define CEA_ABORT_MESSAGE_CAPACITY 4096

typedef struct
{
  int active;
  jmp_buf env;
  char message[CEA_ABORT_MESSAGE_CAPACITY];
} cea_abort_state;

static CEA_THREAD_LOCAL cea_abort_state cea_abort_tls = {0};

void cea_bindc_abort_clear_message(void)
{
  cea_abort_tls.message[0] = '\0';
}

void cea_bindc_abort_set_active(int active)
{
  cea_abort_tls.active = active;
}

jmp_buf *cea_bindc_abort_jmp_buf(void)
{
  return &cea_abort_tls.env;
}

static void cea_abort_store_message(const char *message, cea_int message_len)
{
  size_t ncopy = 0;

  if (message != NULL && message_len > 0)
  {
    ncopy = (size_t)message_len;
    if (ncopy >= sizeof(cea_abort_tls.message))
    {
      ncopy = sizeof(cea_abort_tls.message) - 1u;
    }
    memcpy(cea_abort_tls.message, message, ncopy);
  }
  cea_abort_tls.message[ncopy] = '\0';
}

int cea_bindc_abort_recovery_is_active(void)
{
  return cea_abort_tls.active;
}

void cea_bindc_abort_recover(const char *message, cea_int message_len)
{
  cea_abort_store_message(message, message_len);
  if (cea_abort_tls.active)
  {
    longjmp(cea_abort_tls.env, 1);
  }
  abort();
}

cea_err cea_last_error_message_len(cea_int *message_len)
{
  if (message_len == NULL)
  {
    return CEA_INVALID_SIZE;
  }
  *message_len = (cea_int)strlen(cea_abort_tls.message);
  return CEA_SUCCESS;
}

cea_err cea_last_error_message_buf(char *message, const cea_int buf_len)
{
  size_t msg_len = strlen(cea_abort_tls.message);
  size_t ncopy;

  if (message == NULL || buf_len <= 0)
  {
    return CEA_INVALID_SIZE;
  }

  ncopy = msg_len;
  if (ncopy >= (size_t)buf_len)
  {
    ncopy = (size_t)buf_len - 1u;
  }

  if (ncopy > 0)
  {
    memcpy(message, cea_abort_tls.message, ncopy);
  }
  message[ncopy] = '\0';

  if (msg_len + 1u > (size_t)buf_len)
  {
    return CEA_INVALID_SIZE;
  }
  return CEA_SUCCESS;
}
