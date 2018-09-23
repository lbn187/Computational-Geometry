#ifndef CIRCLE_HPP
#define CIRCLE_HPP
#include "Line.hpp"
#include "Triangle.hpp"
struct C{
	P c;D r;
	C(P c=P(),D r=0):c(c),r(r){}//圆心和半径做圆 
	C(P a,P b):c((a+b)/2),r(dis(a,b)/2){}//两点为直径做圆 
	C(P x,P y,P z){//三点做圆 
		c=CircumCenter(x,y,z);
		r=dis(c,x);
	}
	D area()const{return pi*r*r;}
};
ostream&operator<<(ostream&os,const C&a){
	os<<a.c<<",  r="<<a.r<<endl;
	return os;
}
istream&operator>>(istream&is,C&a){
	is>>a.c>>a.r;
	return is;
}
bool operator==(const C&a,const C&b){return a.c==b.c&&cmp(a.r,b.r)==0;}
bool operator!=(const C&a,const C&b){return !(a==b);}
bool InCircle(const P&a,const C&b){return cmp(dis(a,b.c),b.r)<=0;}
int CircleRelationship(const C&a,const C&b){//两圆关系 
	if(a==b)return 0;//相同 
	D d=dis(a.c,b.c),r1=a.r+b.r,r2=fabs(a.r-b.r);
	if(sgn(d-r1)==1)return 5;//相离 
	if(sgn(d-r1)==0)return 4;//外切 
	if(sgn(d-r2)==1)return 3;//相交 
	if(sgn(d-r2)==0)return 2;//内切 
	return 1;//相含 
}
vector<P> LineCircleIntersect(const L&a,const C&b){//直线与圆交
	D dis=PointToLine(b.c,a),x=sqrt(sqr(b.r)-sqr(dis));
	if(cmp(dis,b.r)>0)return{};
	P p1=ProjectToLine(b.c,a)+(a.s-a.t).unit()*x;
	P p2=ProjectToLine(b.c,a)-(a.s-a.t).unit()*x;
	if(cmp(dis,b.r)==0)return{p1};
	return{p1,p2};
}
vector<P> SegmentCircleIntersect(const L&a,const C&b){//线段与圆交
	vector<P>p=LineCircleIntersect(a,b),an;
	for(auto o:p)if(PointOnSegment(o,a))an.PB(o);
	return an;
}
D CircleIntersectArea(C a,C b){//两圆相交的面积 
	D d=dis(a.c,b.c);
	if(sgn(d-(a.r+b.r))>=0)return 0;
	if(sgn(d-abs(a.r-b.r))<=0)return sqr(min(a.r,b.r))*pi;
	D x=(d*d+a.r*a.r-b.r*b.r)/(d*2),t1=acos(min(1.,max(-1.,x/a.r))),t2=acos(min(1.,max(-1.,(d-x)/b.r)));
	return a.r*a.r+t1+b.r*b.r*t2-d*a.r*sin(t1);
}
vector<P> CircleIntersect(const C&a,const C&b){//圆交
	int tmp=CircleRelationship(a,b);
	if(tmp<2||tmp>4)return{};
	P r=(b.c-a.c).unit();
	D d=dis(a.c,b.c),x=((sqr(a.r)-sqr(b.r))/d+d)/2,h=sqrt(sqr(a.r)-sqr(x));
	if(tmp!=3)return{a.c+r*x};
	return{a.c+r*x+r.rot90()*h,a.c+r*x-r.rot90()*h};
}
vector<P> Tangent(const P&a,const C&b){//点关于圆切点 
	return CircleIntersect(C(a,b.c),b);
}
vector<L> ExTangent(const C&a,const C&b){//两圆外公切线
	if(cmp(dis(a.c,b.c),fabs(a.r-b.r))<=0)return{};
	if(sgn(a.r-b.r)==0){
		P c=b.c-a.c;
		c=(c.unit()*a.r).rot90();
		return{L(a.c+c,b.c+c),L(a.c-c,b.c-c)};
	}else{
		P p=(b.c*a.r-a.c*b.r)/(a.r-b.r);
		vector<P> pp=Tangent(p,a),qq=Tangent(p,b);
		if(SZ(pp)==2&&SZ(qq)==2){
			if(cmp(a.r,b.r)<0)swap(pp[0],pp[1]),swap(qq[0],qq[1]);
			return{L(pp[0],qq[0]),L(pp[1],qq[1])};
		}else return{};
	}
}
vector<L> InTangent(const C&a,const C&b){//两圆内公切线
	P p=(a.c*b.r+b.c*a.r)/(a.r+b.r);
	vector<P> pp=Tangent(p,a),qq=Tangent(p,b);
	if(SZ(pp)==2&&SZ(qq)==2)return{L(pp[0],qq[0]),L(pp[1],qq[1])};
	return{};
}
D SectorArea(const C&c,const P&a,const P&b){//不大于180°的扇形面积 
	D u=(c.c-a).norm2(),v=(c.c-b).norm2(),w=(a-b).norm2();
	return sqr(c.r)*acos((u+v-w)*.5/sqrt(u)/sqrt(v))*.5;
}
C Appollo(const P&a,const P&b,const D&k){//dis(a,o)*k=dis(b,o)
	C o;P p=(a-b).unit();D d=dis(a,b);
	return C(b-p*d/(1.0/k-1.),b+p*d*k/(k+1.));
}
C MinimumCircle(vector<P>p){//最小圆覆盖 
	C a;int n=SZ(p);
	random_shuffle(ALL(p));
	FR(i,0,n-1)if(!InCircle(p[i],a)){
		a=C(p[i],0);
		FR(j,0,i-1)if(!InCircle(p[j],a)){
			a=C(p[j],p[i]);
			FR(k,0,j-1)if(!InCircle(p[k],a))a=C(p[i],p[j],p[k]);
		}
	}
	return a;
}
#endif
