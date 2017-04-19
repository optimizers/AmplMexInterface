Matlab MEX Interface to AMPL
============================

This is a Matlab interface to the [AMPL/MP Library](https://github.com/ampl/mp). It is setup as a Matlab [package](http://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html), which means that the enclosing directory *must* be named `+ampl`.

Clone it this way:

    git clone git@github.com:optimizers/AmplMexInterface.git +ampl

Supports both dense and sparse models.

The easiest way to install AMPL/MP is via [Homebrew](https://brew.sh) (on OSX) or [Linuxbrew](http://linuxbrew.sh) (on Linux):

    brew tap homebrew/science
    brew install ampl-mp

Note: This interface will not compile against the former Netlib ASL.

Enjoy!

Michael Friedlander
(<mpfriedlander@gmail.com>)
<br>
Dominique Orban
(<dominique.orban@gmail.com>)
