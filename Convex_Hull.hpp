#ifndef CONVEX_HULL_HPP
#define CONVEX_HULL_HPP
#include "Point.hpp"
bool cmp(P a,P b,P s){return LeftTest(s,a,b)==1;}
vector<P> ConvexHull(vector<P>a){
	int n=SZ(a),t=0;vector<P>st;
	sort(ALL(a));
	for(auto o:a){
		for(;t>1&&LeftTest(st[t-2],o,st[t-1]);t--)st.pop_back();
		st.PB(o);t++;
	}
	int tmp=t;reverse(ALL(a));
	for(auto o:a){
		for(;t>tmp&&LeftTest(st[t-2],o,st[t-1]);t--)st.pop_back();
		st.PB(o);t++;
	}
	st.pop_back();
	return st;
}
#endif
