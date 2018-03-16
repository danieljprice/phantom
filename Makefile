#----------------------------------------------------------------
#   Parent Makefile for Phantom SPH code
#   This file is just a wrapper for the various sub-makes
#
#   See build/Makefile for the main Makefile
#
#   (c) 2007-2018 The Authors (see AUTHORS)
#
#  $Id: Makefile,v 98b9fad01f38 2013/03/25 23:02:49 daniel $
#----------------------------------------------------------------

.PHONY: phantom
phantom:
	@cd build; ${MAKE} ${MAKECMDGOALS}

%::
	@cd build; ${MAKE} ${MAKECMDGOALS}

clean:
	@cd build; ${MAKE} clean
