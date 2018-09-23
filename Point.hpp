#ifndef POINT_HPP
#define POINT_HPP 
#include "Basic.hpp"
struct P{
	D x,y;
	P()=default;
	P(D x,D y):x(x),y(y){}
	D norm(){return sqrt(x*x+y*y);}
	D norm2(){return x*x+y*y;}
	P unit(){
		D l=norm();
		return P(x/l,y/l);
	}
	P rot90(){return P(-y,x);}
	P _rot90(){return P(y,-x);}
	P rotate(D t){//»¡¶È 
		D c=cos(t),s=sin(t);
		return P(x*c-y*s,x*s+y*c);
	}
	D angle(){return atan2(y,x);}
};
ostream&operator<<(ostream&os,const P&a){
	os<<"("<<a.x<<","<<a.y<<")"<<endl;
	return os;
}
istream&operator>>(istream&is,P&a){
	is>>a.x>>a.y;
	return is;
}
bool operator==(const P&a,const P&b){return !cmp(a.x,b.x)&&!cmp(a.y,b.y);}
bool operator!=(const P&a,const P&b){return !(a==b);}
bool operator<(const P&a,const P&b){return cmp(a.x,b.x)==0?cmp(a.y,b.y)<0:cmp(a.x,b.x)<0;}
bool operator>(const P&a,const P&b){return cmp(a.x,b.x)==0?cmp(a.y,b.y)>0:cmp(a.x,b.x)>0;}
bool operator<=(const P&a,const P&b){return !(a>b);}
bool operator>=(const P&a,const P&b){return !(a<b);}
P operator-(const P&a){return P(-a.x,-a.y);}
P operator+(const P&a,const P&b){return P(a.x+b.x,a.y+b.y);}
P operator-(const P&a,const P&b){return P(a.x-b.x,a.y-b.y);}
template<typename T>P operator*(const P&a,T b){return P(a.x*b,a.y*b);}
template<typename T>P operator/(const P&a,T b){return P(a.x/b,a.y/b);}
D operator*(const P&a,const P&b){return a.x*b.x+a.y*b.y;}
D operator^(const P&a,const P&b){return a.x*b.y-a.y*b.x;}
D dot(const P&a,const P&b){return a.x*b.x+a.y*b.y;}
D det(const P&a,const P&b){return a.x*b.y-a.y*b.x;}
D dis(const P&a,const P&b){return (a-b).norm();}
D Area(const P&a,const P&b,const P&c){return fabs(det(b-a,c-a)*.5);}
D Area2(const P&a,const P&b,const P&c){return fabs(det(b-a,c-a));}
bool LeftTest(const P&a,const P&b,const P&s){return sgn(det(b-a,s-a))>0;}
bool RightTest(const P&a,const P&b,const P&s){return sgn(det(b-a,s-a))<0;}
bool FrontTest(const P&a,const P&b,const P&s){return sgn((b-a)*(s-a))>0;}
bool BehindTest(const P&a,const P&b,const P&s){return sgn((b-a)*(s-a))<0;}
#endif
