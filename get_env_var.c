#include <stdlib.h>

void get_env_var_(const char* env_var_name, char* value, int* length) {
    char* env_value = getenv(env_var_name);
    if (env_value) {
        int i;
        for (i = 0; env_value[i] != '\0' && i < *length; i++) {
            value[i] = env_value[i];
        }
        *length = i;  // Return the actual length of the environment variable value
    } else {
        value[0] = '\0';  // If not found, return an empty string
        *length = 0;
    }
}
