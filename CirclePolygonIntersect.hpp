#include "Circle.hpp"
#include "Line.hpp"
D CircleTriangleIntersect(const C&c,const P&a,const P&b){//圆与一点在圆心的三角形的交的面积 
	if(sgn(det(a-c.c,b-c.c))==0)return 0;
	bool f1=InCircle(a,c),f2=InCircle(b,c);
	P pa,pb;L l=L(a,b);
	if(f1&&f2)return Area(a,b,c.c);//两点都在圆内 
	if(PointToSegment(c.c,l)>c.r-eps)return SectorArea(c,a,b);//两点都在圆外且两点的线段与圆无交 
	vector<P> ret=LineCircleIntersect(l,c);
	if(SZ(ret)!=2)return 0;
	pa=ret[0];pb=ret[1];
	if(f1)return Area(a,pb,c.c)+SectorArea(c,pb,b);//一点在圆内，一点在圆外 
	if(f2)return Area(pa,b,c.c)+SectorArea(c,a,pa);
	return Area(pa,pb,c.c)+SectorArea(c,a,pa)+SectorArea(c,pb,b);//两点都在圆外且两点的线段与圆有交
} 
D CirclePolygonIntersect(const C&c,vector<P>a){//求圆与多边形交
	int n=SZ(a);D an=0;
	a.PB(a[0]);//最后添一个方便处理 
	FR(i,0,n-1){
		D d=CircleTriangleIntersect(c,a[i],a[i+1]);
		if(sgn(det(a[i]-c.c,a[i+1]-c.c))>0)an+=d;else an-=d;//叉积判断正负
	}
	return fabs(an);//返回正值 
}
