#----------------------------------------------------------------
#   Parent Makefile for Phantom SPH code
#   This file is just a wrapper for the various sub-makes
#
#   See build/Makefile for the main Makefile
#
#   (c) 2007-2024 The Authors (see AUTHORS)
#
#----------------------------------------------------------------

.PHONY: phantom
phantom:
	@cd build; ${MAKE} ${MAKECMDGOALS}

%::
	@cd build; ${MAKE} "${MAKECMDGOALS}"

clean:
	@cd build; ${MAKE} clean
