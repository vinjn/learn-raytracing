[via](http://www.kevinbeason.com/smallpt/smallpt.txt)

`smallpt.cpp`
The original version.

`explicit.cpp`
Huge speedup, especially for small lights. Adds explicit light sampling with 23 additional lines of code and a small function signature change. Produces this image in 10 seconds on a Intel Core i7 920 quad-core CPU using 16 samples per pixel.

`forward.cpp`
Revision of radiance() function that removes all recursion and uses only a simple loop and no path branching. That is, the ray tree is always one ray wide.
