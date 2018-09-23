#ifndef LINE_HPP
#define LINE_HPP
#include "Point.hpp"
struct L{
	P s,t;
	L(P s=P(),P t=P()):s(s),t(t){}
	D length()const{return dis(s,t);}
};
ostream&operator<<(ostream&os,const L&a){
	os<<a.s<<" --> "<<a.t<<endl;
	return os;
}
istream&operator>>(istream&is,L&a){
	is>>a.s>>a.t;
	return is;
}
bool LeftTest(const P&a,const L&b){return LeftTest(b.s,b.t,a);}
bool RightTest(const P&a,const L&b){return RightTest(b.s,b.t,a);}
bool FrontTest(const P&a,const L&b){return FrontTest(b.s,b.t,a);}
bool BehindTest(const P&a,const L&b){return BehindTest(b.s,b.t,a);}
bool PointOnLine(const P&a,const L&b){
	return !sgn((a-b.s)*(b.t-b.s));
}
bool PointOnSegment(const P&a,const L&b){
	return !sgn((a-b.s)^(b.t-b.s))&&sgn((b.s-a)*(b.t-a))<=0;
}
bool TwoSide(const P&a,const P&b,const L&c){
	return sgn((a-c.s)^(c.t-c.s))*sgn((b-c.s)^(c.t-c.s))<0;
}
bool SegmentIntersectJudge(const L&a,const L&b){
	if(PointOnSegment(b.s,a)||PointOnSegment(b.t,a))return 1;
	if(PointOnSegment(a.s,b)||PointOnSegment(a.t,b))return 1;
	return TwoSide(a.s,a.t,b)&&TwoSide(b.s,b.t,a);
}
P LineIntersect(const L&a,const L&b){
	D s1=(a.t-a.s)^(b.s-a.s),s2=(a.t-a.s)^(b.t-a.s);
	return (b.s*s2-b.t*s1)/(s2-s1);
}
D PointToLine(const P&a,const L&b){
	return fabs((b.t-b.s)^(a-b.s))/dis(b.s,b.t);
}
P ProjectToLine(const P&a,const L&b){//a在b上投影 
	return b.s+(b.t-b.s)*((a-b.s)*(b.t-b.s)/(b.t-b.s).norm2());
}
P SymmetryPointToLine(const P&a,const L&b){//a关于b对称点
	return ProjectToLine(a,b)*2-a;
}
D PointToSegment(const P&a,const L&b){
	if(sgn(dot(b.s-a,b.t-b.s)*dot(b.t-a,b.t-b.s))<=0)return PointToLine(a,b);
	return min(dis(a,b.s),dis(a,b.t));
}
D AngleCos(const L&a,const L&b){
	return dot(a.t-a.s,b.t-b.s)/a.length()/b.length();
}
#endif
