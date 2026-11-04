# Instructions for coding agents working on Phantom

Phantom's conventions are in the developer guide. Read these before changing code, and follow them:

- [docs/developer-guide/philosophy.rst](docs/developer-guide/philosophy.rst): what Phantom is and its design goals
- [docs/developer-guide/styleguide.rst](docs/developer-guide/styleguide.rst): Fortran style, including subroutines and functions at module level (never contained inside another routine) and no new `#ifdef`s
- [docs/developer-guide/testing.rst](docs/developer-guide/testing.rst): running the test suite; code must compile without warnings (`NOWARN=yes`)
- [docs/developer-guide/bots.rst](docs/developer-guide/bots.rst): run `cd scripts; ./bots.sh --apply` before opening a pull request

Match the surrounding code, reuse existing routines rather than writing new ones, and keep pull requests small and focused.
