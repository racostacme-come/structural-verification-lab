#include "assembly.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

// Runtime checks deliberately remain enabled in Release builds (no assert).
void require(bool value, const char* what) {
    if (!value) throw std::runtime_error(what);
}
void near(double a, double b) {
    require(std::abs(a-b) <= 1e-12*std::max({1.0,std::abs(a),std::abs(b)}), "numeric invariant");
}
int main() {
    try {
        require(tm_abi_version() == 1, "ABI");
        double x[]{0,2,5}, e[]{10,20}, a[]{2,1}, alpha[]{.1,.2}, t[]{1,1};
        double d[3], o[2], f[3];
        auto assemble = [&]() { return tm_assemble(2,x,e,a,alpha,t,d,o,f); };
        require(assemble()==0, "assembly");
        near(d[0],10); near(d[1],10+20.0/3); near(d[2],20.0/3);
        near(o[0],-10); near(o[1],-20.0/3);
        near(f[0],-2); near(f[1],-2); near(f[2],4);
        near(d[0]+o[0],0); near(d[1]+o[0]+o[1],0); near(d[2]+o[1],0);
        near(f[0]+f[1]+f[2],0);
        // Independent exact compatibility: N=-2, u=(.3,.3,.6).
        near(d[0]*.3+o[0]*.3-f[0],2);
        near(o[0]*.3+d[1]*.3+o[1]*.6-f[1],0);
        near(o[1]*.3+d[2]*.6-f[2],-2);
        const double u[]{.4,-.7,.2};
        const double energy=d[0]*u[0]*u[0]+d[1]*u[1]*u[1]+d[2]*u[2]*u[2]
                            +2*o[0]*u[0]*u[1]+2*o[1]*u[1]*u[2];
        near(energy,10*1.1*1.1+(20.0/3)*.9*.9);
        // One-element free expansion and restraint, signed contraction.
        for (double signed_alpha : {-.1,0.0,.1}) {
            alpha[0]=signed_alpha;
            require(tm_assemble(1,x,e,a,alpha,t,d,o,f)==0,"one element");
            near(d[1]*(2*signed_alpha)-f[1],0);
            near(-f[0],20*signed_alpha);
        }
        require(tm_assemble(0,x,e,a,alpha,t,d,o,f)==1,"empty input");
        require(tm_assemble(2,nullptr,e,a,alpha,t,d,o,f)==1,"null input");
        x[1]=0; require(assemble()==2,"zero length"); x[1]=2;
        e[0]=-1; require(assemble()==2,"negative modulus"); e[0]=10;
        a[0]=0; require(assemble()==2,"zero area"); a[0]=2;
        t[0]=std::numeric_limits<double>::quiet_NaN();
        require(assemble()==2,"NaN temperature"); t[0]=1;
        alpha[0]=std::numeric_limits<double>::infinity();
        require(assemble()==2,"infinite expansion"); alpha[0]=.1;
        e[0]=std::numeric_limits<double>::max(); require(assemble()==3,"EA overflow");
        e[0]=std::numeric_limits<double>::min(); a[0]=e[0];
        require(assemble()==3,"stiffness underflow");
        std::cout << "Native thermoelastic assembly: analytical, energy, nullspace, signed patches and rejection checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
