#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP
#include "Line.hpp"
bool InTriangle(const P&a,const P&b,const P&c,const P&s){
	return (LeftTest(a,b,s)&&LeftTest(b,c,s)&&LeftTest(c,a,s))||(LeftTest(b,a,s)&&LeftTest(c,b,s)&&LeftTest(a,c,s));
}
P CenterOfGravity(const P&a,const P&b,const P&c){//重心 
	return P((a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0);
}
P InCenter(const P&a,const P&b,const P&c){//内心 
	D p=dis(a,b)+dis(b,c)+dis(c,a);
	return (a*dis(b,c)+b*dis(c,a)+c*dis(a,b))/p;
}
P CircumCenter(const P&a,const P&b,const P&c){//外心 
	P p=b-a,q=c-a,s(p.norm2()/2,q.norm2()/2);
	D d=p^q;
	return a+P(s^P(p.y,q.y),P(p.x,q.x)^s)/d;
}
P Heart(const P&a,const P&b,const P&c){//垂心 
	return a+b+c-CircumCenter(a,b,c)*2.0;
}
tuple<P,P,P> EsCenter(const P&a,const P&b,const P&c){//旁心 
	D la=dis(b,c),lb=dis(a,c),lc=dis(a,b);
	return MT((-a*la+b*lb+c*lc)/(-la+lb+lc),(a*la-b*lb+c*lc)/(la-lb+lc),(a*la+b*lb-c*lc)/(la+lb-lc));
}
P FermatPoint(const P&a,const P&b,const P&c){//费马点 
	if(a==b)return a;
	if(b==c)return b;
	if(c==a)return c;
	D ab=dis(a,b),bc=dis(b,c),ca=dis(c,a),sq3=pi/3.0;
	D cosa=dot(b-a,c-a)/ab/ca,cosb=dot(a-b,c-b)/ab/bc,cosc=dot(b-c,a-c)/ca/bc;
	P md;
	if(sgn(cosa+0.5)<0)md=a;
	else if(sgn(cosb+0.5)<0)md=b;
	else if(sgn(cosc+0.5)<0)md=c;
	else if(sgn(dot(b-a,c-a))<0)md=LineIntersect(L(a,b+(c-b).rotate(sq3)),L(b,c+(a-c).rotate(sq3)));
	else md=LineIntersect(L(a,c+(b-c).rotate(sq3)),L(c,b+(a-b).rotate(sq3)));
	return md;
}
#endif
