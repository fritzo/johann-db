Johann-DB
---------

A lightweight tool for reading Johann dababases.

The [Johann](http://github.com/fritzo/Johann) system
is a toolset for reasoning about combinatory algebras.
It includes both logical deductive algorithms (theorem-proving)
and statistical analysis algorithms for large combinatory databases.
Johann can be used to build large database of combinators
--essentially multiplication tables for functional programs.

The [Johann-DB](http://github.com/fritzo/johann-db) library
allows other programs to use the combinator databases created by Johann.

Usage
=====

1. include lib/jdb.h in your C++ code

2. create a Database object

3. compile and link lib/jdb.cpp with your C++ target

4. download a database (.jdb file) at http://fritzo.org/johann

Examples
========

see johann-db/examples

