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

#include "params.h"
namespace mpcc{
    
Param::Param(){
    std::cout << "Default initialization of model params" << std::endl;
}

Param::Param(std::string file){
    /////////////////////////////////////////////////////
    // Loading Model and Constraint Parameters //////////
    /////////////////////////////////////////////////////
    // std::cout << "model" << std::endl;

    std::ifstream iModel(file);
    json jsonModel;
    iModel >> jsonModel;
    Br 	= jsonModel["Br"];
    Cr 	= jsonModel["Cr"];
    Dr 	= jsonModel["Dr"];

    Bf 	= jsonModel["Bf"];
    Cf 	= jsonModel["Cf"];
    Df 	= jsonModel["Df"];

    m 	= jsonModel["m"];
    Iz 	= jsonModel["Iz"];
    lf 	= jsonModel["lf"];
    lr 	= jsonModel["lr"];

    car_l = jsonModel["car_l"];
    car_w = jsonModel["car_w"];
    
    g = jsonModel["g"];

    e_long = jsonModel["E_long"];
    e_eps = jsonModel["E_eps"];

    max_alpha = jsonModel["maxAlpha"];
    // initial warm start and trust region (model dependent)
    initial_velocity = jsonModel["initial_velocity"];

    vx_zero = jsonModel["vx_zero"];
}

CostParam::CostParam(){
    std::cout << "Default initialization of cost" << std::endl;
}

CostParam::CostParam(std::string file){
    /////////////////////////////////////////////////////
    // Loading Cost Parameters //////////////////////////
    /////////////////////////////////////////////////////
    // std::cout << "cost" << std::endl;

    std::ifstream iCost(file);
    json jsonCost;
    iCost >> jsonCost;

    q_c = jsonCost["qC"];
    q_l = jsonCost["qL"];
    q_vx = jsonCost["qVx"];
    q_n = jsonCost["qn"];

    q_mu = jsonCost["qMu"];

    q_r = jsonCost["qR"];

    q_beta = jsonCost["qBeta"];
    beta_kin_cost = 1;//jsonCost["betaKin"];

    r_delta = jsonCost["rDelta"];

    r_dDelta = jsonCost["rdDelta"];

    q_c_N_mult = jsonCost["qCNmult"];
    q_r_N_mult = jsonCost["qRNmult"];

    //sc_quad_track = jsonCost["sc_quad_track"];
    //sc_quad_tire= jsonCost["sc_quad_tire"];
    //sc_quad_alpha = jsonCost["sc_quad_alpha"];

    //sc_lin_track = jsonCost["sc_lin_track"];
    //sc_lin_tire = jsonCost["sc_lin_tire"];
    //sc_lin_alpha = jsonCost["sc_lin_alpha"];
}

BoundsParam::BoundsParam() {
    std::cout << "Default initialization of bounds" << std::endl;
}

BoundsParam::BoundsParam(std::string file) {

    /////////////////////////////////////////////////////
    // Loading Cost Parameters //////////////////////////
    /////////////////////////////////////////////////////
    // std::cout << "bounds" << std::endl;

    std::ifstream iBounds(file);
    json jsonBounds;
    iBounds >> jsonBounds;

    lower_state_bounds.n_l = jsonBounds["nl"];
    lower_state_bounds.mu_l = jsonBounds["mul"];
    lower_state_bounds.vx_l = jsonBounds["vxl"];
    lower_state_bounds.vy_l = jsonBounds["vyl"];
    lower_state_bounds.r_l = jsonBounds["rl"];
    lower_state_bounds.Fm_l = jsonBounds["Fml"];
    lower_state_bounds.delta_l = jsonBounds["deltal"];
    lower_state_bounds.k_l = jsonBounds["kl"];

    upper_state_bounds.n_u = jsonBounds["nu"];
    upper_state_bounds.mu_u = jsonBounds["muu"];
    upper_state_bounds.vx_u = jsonBounds["vxu"];
    upper_state_bounds.vy_u = jsonBounds["vyu"];
    upper_state_bounds.r_u = jsonBounds["ru"];
    upper_state_bounds.Fm_u = jsonBounds["Fmu"];
    upper_state_bounds.delta_u = jsonBounds["deltau"];
    upper_state_bounds.k_u = jsonBounds["ku"];

    lower_state_bounds.dFm_l = jsonBounds["dFml"];
    lower_state_bounds.dDelta_l = jsonBounds["dDeltal"];
    lower_state_bounds.Mz_l = jsonBounds["Mzl"];
    lower_state_bounds.dk_l = jsonBounds["dkl"];

    upper_input_bounds.dFm_u = jsonBounds["dFmu"];
    upper_input_bounds.dDelta_u = jsonBounds["dDeltau"];
    upper_input_bounds.Mz_u = jsonBounds["Mzu"];
    upper_input_bounds.dk_u = jsonBounds["dku"];
}

NormalizationParam::NormalizationParam(){
    std::cout << "Default initialization of normalization" << std::endl;
}

NormalizationParam::NormalizationParam(std::string file)
{
    /////////////////////////////////////////////////////
    // Loading Normalization Parameters /////////////////
    /////////////////////////////////////////////////////
    // std::cout << "norm" << std::endl;

    std::ifstream iNorm(file);
    json jsonNorm;
    iNorm >> jsonNorm;

    T_x.setIdentity();
    T_x(si_index.n,si_index.n) = jsonNorm["n"];
    T_x(si_index.mu,si_index.mu) = jsonNorm["mu"];
    T_x(si_index.vx,si_index.vx) = jsonNorm["vx"];
    T_x(si_index.vy,si_index.vy) = jsonNorm["vy"];
    T_x(si_index.r,si_index.r) = jsonNorm["r"];
    T_x(si_index.Fm,si_index.Fm) = jsonNorm["Fm"];
    T_x(si_index.delta,si_index.delta) = jsonNorm["delta"];
    T_x(si_index.k,si_index.k) = jsonNorm["k"];


    T_x_inv.setIdentity();
    for(int i = 0;i<NX;i++)
    {
        T_x_inv(i,i) = 1.0/T_x(i,i);
    }

    T_u.setIdentity();
    T_u(si_index.dFm,si_index.dFm) = jsonNorm["dFm"];
    T_u(si_index.dDelta,si_index.dDelta) = jsonNorm["dDelta"];
    T_u(si_index.Mz,si_index.Mz) = jsonNorm["dMz"];
    T_u(si_index.dk,si_index.dk) = jsonNorm["dk"];

    T_u_inv.setIdentity();
    for(int i = 0;i<NU;i++)
    {
        T_u_inv(i,i) = 1.0/T_u(i,i);
    }

    T_s.setIdentity();
    T_s_inv.setIdentity();
}

}
