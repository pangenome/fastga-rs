// Rust-friendly wrappers for ONElib functions
// These provide clean accessors without macro complexity

#include "ONElib.h"
#include "alncode.h"

// Field accessors
int64_t one_int(OneFile *of, int index) {
    return oneInt(of, index);
}

double one_real(OneFile *of, int index) {
    return oneReal(of, index);
}

char one_char(OneFile *of, int index) {
    return oneChar(of, index);
}

// Get current line type
char one_line_type(OneFile *of) {
    return of->lineType;
}

// Get line count for this line type
int64_t one_line_count(OneFile *of) {
    return of->line;
}
