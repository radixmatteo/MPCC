// Copyright 2019 Alexander Liniger

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#include "model.h"
namespace mpcc{
Model::Model()
:Ts_(1.0)
{
    std::cout << "default constructor, not everything is initialized properly" << std::endl;
}

Model::Model(double Ts,const PathToJson &path)
:Ts_(Ts),param_(Param(path.param_path))
{
}

double Model::getSlipAngleFront(const State &x) const
{
    // compute slip angels given current state
    return -std::atan2(x.vy+x.r*param_.lf,x.vx) + x.delta;
}

double Model::getSlipAngleRear(const State &x) const
{
    // compute slip angels given current state
    return -std::atan2(x.vy-x.r*param_.lr,x.vx);
}

TireForces Model::getForceFront(const State &x) const
{
    const double alpha_f = getSlipAngleFront(x);
    const double F_y = param_.Df * std::sin(param_.Cf * std::atan(param_.Bf * alpha_f ));
    const double F_x = 0.0;

    return {F_y,F_x};
}

TireForces Model::getForceRear(const State &x) const
{
    const double alpha_r = getSlipAngleRear(x);
    const double F_y = param_.Dr * std::sin(param_.Cr * std::atan(param_.Br * alpha_r ));
    const double F_x = x.Fm;// - param_.Cr0 - param_.Cr2*std::pow(x.vx,2.0);

    return {F_y,F_x};
}

double Model::getForceFriction(const State &x) const
{
    return -param_.Cr0 - param_.Cr2*std::pow(x.vx,2.0);
}

NormalForces Model::getForceNormal(const State &x) const
{
    // at this point aero forces could be modeled
    const double f_n_front = param_.lr/(param_.lf + param_.lr)*param_.m*param_.g;
    const double f_n_rear = param_.lf/(param_.lf + param_.lr)*param_.m*param_.g;
    return {f_n_front,f_n_rear};
}

StateVector Model::getF(const State &x,const Input &u) const
{
    const double n = x.n;
    const double mu = x.mu;
    const double vx = x.vx;
    const double vy = x.vy;
    const double r  = x.r;
    const double Fm = x.Fm;
    const double delta = x.delta;
    const double k = x.k;

    const double dFm = u.dFm;
    const double dDelta = u.dDelta;
    const double Mz = u.Mz;
    const double dk = u.dk;

    const TireForces tire_forces_front = getForceFront(x);
    const TireForces tire_forces_rear  = getForceRear(x);
    const double friction_force = getForceFriction(x);

    StateVector f;
    f(0) = vx * sin(mu) + vy * cos(mu);
    f(1) = r - k * (vx * cos(mu) - vy * sin (mu)/(1-n*k));
    f(2) = 1.0/param_.m*(Fm * (1+cos(delta)) - friction_force - tire_forces_front.F_y*std::sin(delta) + param_.m*vy*r);
    f(3) = 1.0/param_.m*(Fm * sin(delta) + tire_forces_rear.F_y + tire_forces_front.F_y*std::cos(delta) - param_.m*vx*r);
    f(4) = 1.0/param_.Iz*(Fm * sin(delta)*param_.lf + tire_forces_front.F_y*param_.lf*std::cos(delta) - tire_forces_rear.F_y*param_.lr + Mz);
    f(5) = dFm;
    f(6) = dDelta;

    return f;
}

LinModelMatrix Model::getModelJacobian(const State &x, const Input &u) const
{
    // compute jacobian of the model
    // state values
    const double n = x.n;
    const double mu = x.mu;
    const double vx = x.vx;
    const double vy = x.vy;
    const double r  = x.r;
    const double Fm = x.Fm;
    const double delta = x.delta;
    const double k = x.k;

    const double dFm = u.dFm;
    const double dDelta = u.dDelta;
    const double Mz = u.Mz;
    const double dk = u.dk;

//    LinModelMatrix lin_model_c;
    A_MPC A_c = A_MPC::Zero();
    B_MPC B_c = B_MPC::Zero();
    g_MPC g_c = g_MPC::Zero();

    const StateVector f = getF(x,u);

    // Derivatives of function
    // f1 = v_x*std::cos(phi) - v_y*std::sin(phi)
    const double df1_dmu  = vx*std::cos(mu) - vy*std::sin(mu);
    const double df1_dvx  = std::sin(mu);
    const double df1_dvy  = std::cos(mu);


    // f2 = v_y*std::cos(phi) + v_x*std::sin(phi);
    const double df2_dn = vy*pow(k,2.0)*sin(mu)*1.0/pow(n*k-1.0,2.0);
    const double df2_dmu  = k*(vx*sin(mu)-(vy*cos(mu))/(n*k-1.0));;
    const double df2_dvx  = -cos(mu)*k;
    const double df2_dvy  = -(k*sin(mu))/(n*k-1.0);
    const double df2_dr   = 1;
    const double df2_dFm  = 0;
    const double df2_ddelta = 0;
    const double df2_k = -vx*cos(mu)-(vy*sin(mu))/(k*n-1.0)+k*n*vy*sin(mu)*1.0/pow(k*n-1.0,2.0);


    // f3 = r;
    const double df3_dvx  = (param_.Cr2*vx*2.0-(param_.Bf*param_.Cf*param_.Df*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))*sin(delta)*(vy+param_.lf*r))/(((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx)))/param_.m;
    const double df3_dvy  = (param_.m*r+(param_.Bf*param_.Cf*param_.Df*vx*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))*sin(delta))/(((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx)))/param_.m;
    const double df3_dr   = (param_.m*vy+(param_.Bf*param_.Cf*param_.Df*param_.lf*vx*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))*sin(delta))/(((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx)))/param_.m;
    const double df3_dFm   = (cos(delta)+1.0)/param_.m;
    const double df3_ddelta   = -(Fm*sin(delta)+param_.Df*cos(delta)*sin(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))+(param_.Bf*param_.Cf*param_.Df*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))*sin(delta))/((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0))/param_.m;

    // f4 = 1/param_.m*(F_rx + F_fric - F_fy*std::sin(delta) + param_.m*v_y*r);
    const double df4_dvx     = -(param_.m*r-(param_.Bf*param_.Cf*param_.Df*cos(delta)*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))*(vy+param_.lf*r))/(((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx)))/param_.m;
    const double df4_dvy     = -(param_.Bf*param_.Cf*param_.Df*vx*cos(delta)*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx)))))/(param_.m*((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx));
    const double df4_dr      = -(param_.m*vx+(param_.Bf*param_.Cf*param_.Df*param_.lf*vx*cos(delta)*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx)))))/(((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx)))/param_.m;
    const double df4_dFm     = sin(delta)/param_.m;
    const double df4_ddelta  = (Fm*cos(delta)-param_.Df*sin(delta)*sin(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))+(param_.Bf*param_.Cf*param_.Df*cos(delta)*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx)))))/((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0))/param_.m;

    // f5 = 1/param_.m*(F_ry + F_fy*std::cos(delta) - param_.m*v_x*r);
    const double df5_dvx     = (param_.Bf*param_.Cf*param_.Df*param_.lf*cos(delta)*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))*(vy+param_.lf*r))/(param_.Iz*((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx));
    const double df5_dvy     = -(param_.Bf*param_.Cf*param_.Df*param_.lf*vx*cos(delta)*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx)))))/(param_.Iz*((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx));
    const double df5_dr      = -(param_.Bf*param_.Cf*param_.Df*(param_.lf*param_.lf)*vx*cos(delta)*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx)))))/(param_.Iz*((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0)*(pow(vy+param_.lf*r,2.0)+vx*vx));
    const double df5_dFm     = (param_.lf*sin(delta))/param_.Iz;
    const double df5_ddelta  = (Fm*param_.lf*cos(delta)-param_.Df*param_.lf*sin(delta)*sin(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx))))+(param_.Bf*param_.Cf*param_.Df*param_.lf*cos(delta)*cos(param_.Cf*atan(param_.Bf*(delta-atan2(vy+param_.lf*r,vx)))))/((param_.Bf*param_.Bf)*pow(delta-atan2(vy+param_.lf*r,vx),2.0)+1.0))/param_.Iz;
    const double df5_duMz     = 1.0/param_.Iz;

    // f6 = 1/param_.Iz*(F_fy*l_f*std::cos(delta)- F_ry*l_r)
    const double df6_duFm     = 1.0;

    // f7
    const double df7_dudelta     = 1.0;

    // f8
    const double df8_dk     = 1.0;


    // Jacobians
    // Column 0
    A_c(1,0) = df2_dn;
    // Column 1
    A_c(0,1) = df1_dmu;
    A_c(1,1) = df2_dmu;
    // Column 2
    A_c(0,2) = df1_dvx;
    A_c(1,2) = df2_dvx;
    A_c(2,2) = df3_dvx;
    A_c(3,2) = df4_dvx;
    A_c(4,2) = df5_dvx;
    // Column 3
    A_c(0,3) = df1_dvy;
    A_c(1,3) = df2_dvy;
    A_c(2,3) = df3_dvy;
    A_c(3,3) = df4_dvy;
    A_c(4,3) = df5_dvy;
    // Column 4
    A_c(1,4) = df2_dr;
    A_c(2,4) = df3_dr;
    A_c(3,4) = df4_dr;
    A_c(4,4) = df5_dr;
    // Column 5
    A_c(1,5) = df2_dFm;
    A_c(2,5) = df3_dFm;
    A_c(3,5) = df4_dFm;
    A_c(4,5) = df5_dFm;
    // Column 6
    A_c(1,6) = df2_ddelta;
    A_c(2,6) = df3_ddelta;
    A_c(3,6) = df4_ddelta;
    A_c(4,6) = df5_ddelta;
    // Column 7
    A_c(1,7) = df2_k;

    // Matrix B
    // Column 1
    B_c(4,2) = df5_duMz;
    // Column 2
    B_c(5,0) = df6_duFm;
    // Column 3
    B_c(6,1) = df7_dudelta;

    B_c(7,3) = df8_dk;

    //zero order term
    g_c = f - A_c*stateToVector(x) - B_c*inputToVector(u);

    return {A_c,B_c,g_c};
}

// Mixed RK4 - EXPM method
LinModelMatrix Model::discretizeModel(const LinModelMatrix &lin_model_c, const State &x, const Input &u,const State &x_next) const
{
    // disctetize the continuous time linear model \dot x = A x + B u + g using ZHO
    Eigen::Matrix<double,NX+NU,NX+NU> temp = Eigen::Matrix<double,NX+NU,NX+NU>::Zero();
    // building matrix necessary for expm
    // temp = Ts*[A,B,g;zeros]
    temp.block<NX,NX>(0,0) = lin_model_c.A;
    temp.block<NX,NU>(0,NX) = lin_model_c.B;
    temp = temp*Ts_;

    // take the matrix exponential of temp
    const Eigen::Matrix<double,NX+NU,NX+NU> temp_res = temp.exp();
    // extract dynamics out of big matrix
    // x_{k+1} = Ad x_k + Bd u_k
    //temp_res = [Ad,Bd;zeros]
    const A_MPC A_d = temp_res.block<NX,NX>(0,0);
    const B_MPC B_d = temp_res.block<NX,NU>(0,NX);

    // TODO: use correct RK4 instead of inline RK4
    const StateVector x_vec = stateToVector(x);

    const StateVector k1 = getF(vectorToState(x_vec),u);
    const StateVector k2 = getF(vectorToState(x_vec+Ts_/2.*k1),u);
    const StateVector k3 = getF(vectorToState(x_vec+Ts_/2.*k2),u);
    const StateVector k4 = getF(vectorToState(x_vec+Ts_*k3),u);
    // combining to give output
    const StateVector x_RK = x_vec + Ts_*(k1/6.+k2/3.+k3/3.+k4/6.);

    const g_MPC g_d =  -stateToVector(x_next) + x_RK;

    // return {A_d,B_d,g_d};
    return {A_d,B_d,g_d};
}

// // EF Method (not suited for used input lifting technique)
// LinModelMatrix Model::discretizeModel(const LinModelMatrix &lin_model_c, const State &x, const Input &u,const State &x_next) const
// {
//     // disctetize the continuous time linear model \dot x = A x + B u + g using ZHO
//     Eigen::Matrix<double,NX+NU,NX+NU> temp = Eigen::Matrix<double,NX+NU,NX+NU>::Zero();
//     // building matrix necessary for expm
//     // temp = Ts*[A,B,g;zeros]
//     temp.block<NX,NX>(0,0) = lin_model_c.A;
//     temp.block<NX,NU>(0,NX) = lin_model_c.B;
//     temp = temp*Ts_;

//     Eigen::Matrix<double,NX+NU,NX+NU> eye;
//     eye.setIdentity();
//     // take the matrix exponential of temp
//     const Eigen::Matrix<double,NX+NU,NX+NU> temp_res = eye + temp;
//     // extract dynamics out of big matrix
//     // x_{k+1} = Ad x_k + Bd u_k
//     //temp_res = [Ad,Bd;zeros]
//     const A_MPC A_d = temp_res.block<NX,NX>(0,0);
//     const B_MPC B_d = temp_res.block<NX,NU>(0,NX);

//     const StateVector x_vec = stateToVector(x);

//     const StateVector f = getF(vectorToState(x_vec),u);
//     // combining to give output
//     const StateVector x_EF = x_vec + Ts_*f;

//     const g_MPC g_d =  -stateToVector(x_next) + x_EF;
//     // return {A_d,B_d,g_d};
//     return {A_d,B_d,g_d};
// }

LinModelMatrix Model::getLinModel(const State &x, const Input &u, const State &x_next) const
{
    // compute linearized and discretized model
    const LinModelMatrix lin_model_c = getModelJacobian(x,u);
    // discretize the system
    return discretizeModel(lin_model_c,x,u,x_next);
}
}