/* FastGA Embedded Library
 *
 * This wraps FastGA's main function to be callable as a library function
 * from Rust, eliminating the need for external binaries.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

// Declare FastGA's main function (renamed to avoid conflicts)
int fastga_main(int argc, char *argv[]);

// Buffer for capturing output
static char* output_buffer = NULL;
static size_t output_size = 0;
static size_t output_capacity = 0;
static FILE* output_stream = NULL;

// Custom write function to capture output
static ssize_t capture_write(void* cookie, const char* buf, size_t size) {
    (void)cookie;  // Unused

    if (output_size + size > output_capacity) {
        size_t new_capacity = (output_capacity + size) * 2;
        char* new_buffer = realloc(output_buffer, new_capacity);
        if (!new_buffer) return -1;
        output_buffer = new_buffer;
        output_capacity = new_capacity;
    }

    memcpy(output_buffer + output_size, buf, size);
    output_size += size;
    return size;
}

// Cookie functions for custom FILE stream
static cookie_io_functions_t capture_funcs = {
    .read = NULL,
    .write = capture_write,
    .seek = NULL,
    .close = NULL
};

// Initialize output capture
void fastga_init_capture() {
    output_buffer = malloc(1024 * 1024);  // Start with 1MB
    output_capacity = 1024 * 1024;
    output_size = 0;

    // Create custom FILE* that writes to our buffer
    output_stream = fopencookie(NULL, "w", capture_funcs);
}

// Get captured output
const char* fastga_get_output() {
    if (output_buffer && output_size < output_capacity) {
        output_buffer[output_size] = '\0';
    }
    return output_buffer;
}

// Clean up capture
void fastga_cleanup_capture() {
    if (output_stream) {
        fclose(output_stream);
        output_stream = NULL;
    }
    if (output_buffer) {
        free(output_buffer);
        output_buffer = NULL;
    }
    output_size = 0;
    output_capacity = 0;
}

// Main entry point callable from Rust
int fastga_run_embedded(int argc, char** argv, char** output) {
    // Save original stdout
    FILE* orig_stdout = stdout;

    // Initialize capture
    fastga_init_capture();

    // Redirect stdout to our capture stream
    if (output && output_stream) {
        stdout = output_stream;
    }

    // Call FastGA's main function
    int result = fastga_main(argc, argv);

    // Flush output
    if (output_stream) {
        fflush(output_stream);
    }

    // Restore original stdout
    stdout = orig_stdout;

    // Return captured output
    if (output) {
        const char* captured = fastga_get_output();
        if (captured) {
            *output = strdup(captured);
        } else {
            *output = NULL;
        }
    }

    // Clean up
    fastga_cleanup_capture();

    return result;
}

// Simple version without output capture
int fastga_run_simple(int argc, char** argv) {
    return fastga_main(argc, argv);
}

// Now include FastGA.c but rename main to fastga_main
#define main fastga_main
#include "FastGA.c"
#undef main