Fortran style guide for Phantom
===============================

Unsure how to format your routine? How many spaces to indent? Want to
make your code fit in with the rest of Phantom? We present the Phantom
f90 style guide…

Conventions which are enforced by the nightly bots
--------------------------------------------------

You can run the bots yourself as follows:

::

   cd phantom/scripts
   ./bots.sh --apply

Indentation
~~~~~~~~~~~

Indentation is enforced automatically by the indent-bot, which uses the
`findent <https://sourceforge.net/projects/findent/>`__ tool. The
specific command used to indent the code is:

::

   findent -r1 -m1 -c3 -Rr -C- -k- -j1 < file.f90

Use modern Fortran
------------------

Do not use SHOUT CASE. THERE IS NO NEED FOR SHOUTING.

endif, enddo
------------

Use endif and enddo, not end if or end do

::

   if (blah) then
      ...
   endif

not

::

   if (blah) then
      ...
   end if

Conventions which are not enforced (but may be in future if we could do it safely)
----------------------------------------------------------------------------------

if statements
~~~~~~~~~~~~~

Use a single space between the if and the bracket. Also use a space
either side of logical operators like .and. or .or.

::

   if (rin >= rup .or. rin < rlow) then

not

::

   if(rin >= rup.or.rin < rlow) then

spacing
~~~~~~~

Use a single space either side of equals sign:

::

   x = 3

not

::

   x=3
   x= 3
   x =3

Use a single space between type declaration and variables:

::

   real :: x

not

::

   real::x
   real ::x
   real:: x

Line continuation
~~~~~~~~~~~~~~~~~

Continue lines with a single ampersand at the end of the line to be
continued. Do not put another ampersand at the start of the next line:

::

   real :: x,y,z,averylongvariable,anotherlongvariable, &
           b,c,d

not

::

   real :: x,y,z,averylongvariable,anotherlongvariable, &
         & b,c,d
