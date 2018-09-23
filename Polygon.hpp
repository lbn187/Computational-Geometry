#ifndef POLYGON_HPP
#define POLYGON_HPP 
#include "Line.hpp"
bool InPolygon(P p,vector<P>E){
	int n=SZ(E),cnt=0,i;
	FR(i,0,n-1){
		P a=E[i],b=E[(i+1)%n];
		if(PointOnSegment(p,L(a,b)))return 1;
		int x=sgn(det(p-a,b-a)),y=sgn(a.y-p.y),z=sgn(b.y-p.y);
		if(x>0&&y<=0&&z>0)cnt++;
		if(x<0&&z<=0&&y>0)cnt--;
	}
	return cnt!=0;
}
D Area(vector<P> a){
	D an=0.0;int n=SZ(a);
	FR(i,0,n-1)an+=det(a[i],a[(i+1)%n])*0.5;
	return an;
}
P CenterOfGravity(vector<P> a){
	int n=SZ(a);P G=P(0,0);
	D s=Area(a);
	FR(i,0,n-1)G=G+(a[i]+a[(i+1)%n])*det(a[i],a[(i+1)%n]);
	return G/6.0;
}
#endif
