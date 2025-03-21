# Third-Party Licenses

This project includes the following third-party libraries and system dependencies. All licenses are attributed to their respective owners.

---

## htslib
Source: https://github.com/samtools/htslib  
License:

The MIT/Expat License

Copyright (C) 2012-2024 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.


[Files within the cram/ subdirectory in this distribution are distributed
according to the terms of the following Modified 3-Clause BSD license.]

The Modified-BSD License

Copyright (C) 2012-2024 Genome Research Ltd.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
   nor the names of its contributors may be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


[The use of a range of years within a copyright notice in this distribution
should be interpreted as being equivalent to a list of years including the
first and last year specified and all consecutive years between them.

For example, a copyright notice that reads "Copyright (C) 2005, 2007-2009,
2011-2012" should be interpreted as being identical to a notice that reads
"Copyright (C) 2005, 2007, 2008, 2009, 2011, 2012" and a copyright notice
that reads "Copyright (C) 2005-2012" should be interpreted as being identical
to a notice that reads "Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010,
2011, 2012".]

---

## libBigWig
Source: https://github.com/dpryan79/libBigWig  
License:

The MIT License (MIT)

Copyright (c) 2015 Devon Ryan

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## zlib
Source: https://zlib.net/  
License: 

zlib.h -- interface of the 'zlib' general purpose compression library
  version 1.3.1, January 22nd, 2024

  Copyright (C) 1995-2024 Jean-loup Gailly and Mark Adler

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Jean-loup Gailly        Mark Adler
  jloup@gzip.org          madler@alumni.caltech.edu
  


## libbz2 / bzip2
Source: https://sourceware.org/bzip2/  
License: BSD-like License  

> Redistribution and use in source and binary forms, with or without modification, are permitted...

---

## liblzma / xz
Source: https://tukaani.org/xz/  
License: 

XZ Utils Licensing
   3 ==================
   4 
   5     Different licenses apply to different files in this package. Here
   6     is a summary of which licenses apply to which parts of this package:
   7 
   8       - liblzma is under the BSD Zero Clause License (0BSD).
   9 
  10       - The command line tools xz, xzdec, lzmadec, and lzmainfo are
  11         under 0BSD except that, on systems that don't have a usable
  12         getopt_long, GNU getopt_long is compiled and linked in from the
  13         'lib' directory. The getopt_long code is under GNU LGPLv2.1+.
  14 
  15       - The scripts to grep, diff, and view compressed files have been
  16         adapted from GNU gzip. These scripts (xzgrep, xzdiff, xzless,
  17         and xzmore) are under GNU GPLv2+. The man pages of the scripts
  18         are under 0BSD; they aren't based on the man pages of GNU gzip.
  19 
  20       - Most of the XZ Utils specific documentation that is in
  21         plain text files (like README, INSTALL, PACKAGERS, NEWS,
  22         and ChangeLog) are under 0BSD unless stated otherwise in
  23         the file itself. The files xz-file-format.txt and
  24         lzma-file-format.xt are in the public domain but may
  25         be distributed under the terms of 0BSD too.
  26 
  27       - Translated messages and man pages are under 0BSD except that
  28         some old translations are in the public domain.
  29 
  30       - Test files and test code in the 'tests' directory, and
  31         debugging utilities in the 'debug' directory are under
  32         the BSD Zero Clause License (0BSD).
  33 
  34       - The GNU Autotools based build system contains files that are
  35         under GNU GPLv2+, GNU GPLv3+, and a few permissive licenses.
  36         These files don't affect the licensing of the binaries being
  37         built.
  38 
  39       - The 'extra' directory contains files that are under various
  40         free software licenses. These aren't built or installed as
  41         part of XZ Utils.
  42 
  43     The following command may be helpful in finding per-file license
  44     information. It works on xz.git and on a clean file tree extracted
  45     from a release tarball.
  46 
  47         sh build-aux/license-check.sh -v
  48 
  49     For the files under the BSD Zero Clause License (0BSD), if
  50     a copyright notice is needed, the following is sufficient:
  51 
  52         Copyright (C) The XZ Utils authors and contributors
  53 
  54     If you copy significant amounts of 0BSD-licensed code from XZ Utils
  55     into your project, acknowledging this somewhere in your software is
  56     polite (especially if it is proprietary, non-free software), but
  57     it is not legally required by the license terms. Here is an example
  58     of a good notice to put into "about box" or into documentation:
  59 
  60         This software includes code from XZ Utils <https://tukaani.org/xz/>.
  61 
  62     The following license texts are included in the following files:
  63       - COPYING.0BSD: BSD Zero Clause License
  64       - COPYING.LGPLv2.1: GNU Lesser General Public License version 2.1
  65       - COPYING.GPLv2: GNU General Public License version 2
  66       - COPYING.GPLv3: GNU General Public License version 3
  67 
  68     If you have questions, don't hesitate to ask for more information.
  69     The contact information is in the README file.
  
---

## libcurl
Source: https://curl.se/libcurl/  
License: MIT / curl license  
Copyright: Daniel Stenberg, et al.

> Permission to use, copy, modify, and distribute this software for any purpose with or without fee...

---

## git

Source: https://git-scm.com/  
icense: GNU General Public License v2  

> This program is free software; you can redistribute it and/or modify...

---

## g++ (GNU Compiler Collection)

Source: https://gcc.gnu.org/  
License: GNU General Public License v3  

> This program is free software; you can redistribute it and/or modify...

---

## make (GNU Make)

Source: https://www.gnu.org/software/make/  
License: GNU General Public License v3  

---

## tclap

Source: https://github.com/mirror/tclap
License:

Copyright (c) 2003 Michael E. Smoot 
Copyright (c) 2004 Daniel Aarno
Copyright (c) 2017 Google Inc.

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be 
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR 
IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
THE SOFTWARE.

---

## nlohmann/json

Source: https://github.com/nlohmann/json/tree/626e7d61e44dee32887126c8f437dd077dec09cf  
License: MIT License  

MIT License 

Copyright (c) 2013-2025 Niels Lohmann

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---



