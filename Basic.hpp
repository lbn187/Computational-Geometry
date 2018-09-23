#ifndef BASIC_HPP
#define BASIC_HPP
#include "start.hpp"
using D=double;//¿ÉÒÔÐÞ¸Ä 
const D eps=1e-8;
const double pi=acos(-1);
int sgn(D x){return x<-eps?-1:x>eps;}
int cmp(D x,D y){return sgn(x-y);}
#endif
