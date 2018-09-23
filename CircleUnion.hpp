#include "Circle.hpp"
D cal(C o,D l,D r){
	return o.r*.5*(o.r*(r-l)+(o.c.x*(sin(r)-sin(l))-o.c.y*(cos(r)-cos(l))));
}
D CircleUnionArea(vector<C>p){//返回圆并面积 
	int n=SZ(p),t;
	pair<D,D>q[N];
	bool vs[N];D v1,v2,Rv,an=0;
	FR(i,0,n-1)vs[i]=0;
	FR(i,0,n-1)FR(j,0,n-1)if(i!=j&&sgn(p[j].r-p[i].r)>=0&&!vs[j])
		if(CircleRelationship(p[i],p[j])<=2)vs[i]=1;
	FR(i,0,n-1)if(!vs[i]){
		t=0;
		FR(j,0,n-1)if(!vs[j]&&i!=j&&CircleRelationship(p[i],p[j])==3){
			vector<P> ret=CircleIntersect(p[i],p[j]);
			v1=(ret[1]-p[i].c).angle(),v2=(ret[0]-p[i].c).angle();
			if(sgn(v1)<0)v1+=pi*2;if(sgn(v2)<0)v2+=pi*2;
			if(sgn(v2-v1)>=0)q[++t]=MP(v1,v2);
			else q[++t]=MP(0,v2),q[++t]=MP(v1,pi*2);
		}
		sort(q+1,q+t+1);Rv=0;
		fr(j,t){
			if(q[j].X>Rv)an+=cal(p[i],Rv,q[j].X);
			Rv=max(Rv,q[j].Y);
		}
		an+=cal(p[i],Rv,pi*2);
	}
	return an;
}
D an[N];//an[k]返回被覆盖k次面积 
void CircleUnion(vector<C>p){
	int n=p.size(),V,j,k,t;
	pair<D,int>q[N];
	bool vs[N];D v1,v2,Rv;
	FR(i,0,n-1)vs[i]=0;
	FR(i,0,n-1)FR(j,0,n-1)if(i!=j&&sgn(p[j].r-p[i].r)>=0&&!vs[j])
		if(CircleRelationship(p[i],p[j])==0)vs[i]=1;
	FR(i,0,n-1)if(!vs[i]){
		for(k=j=0;j<n;j++)if(sgn(p[j].r-p[i].r)>=0&&CircleRelationship(p[i],p[j])<=2)k++;
		for(t=j=0;j<n;j++)if(!vs[j]&&i!=j&&CircleRelationship(p[i],p[j])==3){
			vector<P> ret=CircleIntersect(p[i],p[j]);
			v1=(ret[1]-p[i].c).angle(),v2=(ret[0]-p[i].c).angle();
			if(sgn(v1)<0)v1+=pi*2;if(sgn(v2)<0)v2+=pi*2;
			if(sgn(v2-v1)>=0)q[++t]=MP(v1,1),q[++t]=MP(v2,-1);
			else q[++t]=MP(0,1),q[++t]=MP(v2,-1),q[++t]=MP(v1,1),q[++t]=MP(pi*2,-1);
		}
		sort(q+1,q+t+1);Rv=0;
		fr(j,t){
			an[V+k]+=cal(p[i],q[j-1].first,q[j].first);
			V+=q[j].second;Rv=q[j].first;
		}
		an[V+k]+=cal(p[i],Rv,pi*2);
	}
	fr(i,n-1)an[i]-=an[i+1];
}
