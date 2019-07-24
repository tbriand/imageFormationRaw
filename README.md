# Image Formation from a Large Sequence of RAW Images

## Summary

A C, C++ and bash implementation of the image formation method from RAW images
proposed by Thibaud Briand in 2018 (PhD thesis).

Link to the PDF: https://pastel.archives-ouvertes.fr/tel-01980492 

## Authors ##

* Thibaud Briand <thibaud.briand@enpc.fr>

Laboratoire d'Informatique Gaspard Monge (LIGM)/
Ecole des Ponts ParisTech

## License ##

This program is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

You should have received a copy of the GNU General Pulic License
along with this program. If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2014-2019, Thibaud Briand <thibaud.briand@enpc.fr>

All rights reserved.

## Build ##

Required environment: Any unix-like system with a standard compilation
environment (make and C, C++ compiler) and [Cmake](https://cmake.org/).

Required libraries:
[libpng](http://libpng.org/pub/png/libpng.html),
[lipjpeg](http://ijg.org/),
[libtiff](http://simplesystems.org/libtiff/),
[libfftw3](http://www.fftw.org/)

Optional libraries:
[libgomp](https://gcc.gnu.org/projects/gomp/)

All the build instructions are contained in the build.sh script.
It produces the "build/" directory with all the required binaries.

## Usage ##

The program can be used through the image_formation.sh script.
The input filenames must be formatted as FILENAME%i.EXT where %i represents the number of the image. 
The reference image must be the first one.

  usage:./image_formation.sh in_path ind_ini ind_end output_image raw zoom crop_size free
      
  example:./image_formation.sh PA%i.ORF 1 100 out.tiff 1 1 512 0

## Files in the repository ##

* build.sh           : script for the compilation
* CMakeLists.txt     : cmake file
* external           : directory containing external programs (dcraw, ponomarenko, modified inverse compositional)
* image_formation.sh : main script for using the method
* LICENSE            : license file
* Readme.md          : readme file
* scripts            : directory containing the scripts for the main steps of the method
* src                : directory containing the source code

## Acknowledgements ##

* Pascal Monasse <monasse@imagine.enpc.fr>
* Enric Meinhardt-Llopis <enric.meinhardt@cmla.ens-cachan.fr>
* Jérémy Anger <anger@cmla.ens-cachan.fr>
* Jean-Michel Morel <morel@cmla.ens-cachan.fr>
* Gabriele Facciolo <gfacciol@gmail.com>
