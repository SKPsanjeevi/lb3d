//! Extendable and human-readable verbosity values.
enum VerbosityLevel {
  NONE   = 0, //!< No writing to stdout whatsoever.
  RESULT = 1, //!< Print final result of a function.
  REPORT = 2, //!< Print intermediate steps and final result.
  DEBUG  = 3  //!< Print debug information, intermediate steps and final result.
};
