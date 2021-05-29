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

#include "types.h"
namespace mpcc{

StateVector stateToVector(const State &x)
{
    StateVector xk;
    xk(0) = x.n;
    xk(1) = x.mu;
    xk(2) = x.vx;
    xk(3) = x.vy;
    xk(4) = x.r;
    xk(5) = x.Fm;
    xk(6) = x.delta;
    xk(7) = x.k;
    return xk;
}

InputVector inputToVector(const Input &u)
{
    InputVector uk = {u.dFm,u.dDelta,u.Mz,u.dk};
    return uk;
}

State vectorToState(const StateVector &xk)
{
    State x;
    x.n     = xk(0);
    x.mu     = xk(1);
    x.vx    = xk(2);
    x.vy    = xk(3);
    x.r     = xk(4);
    x.Fm     = xk(5);
    x.delta = xk(6);
    x.k    = xk(7);

    return x;
}

Input vectorToInput(const InputVector &uk)
{
    Input u;
    u.dFm     = uk(0);
    u.dDelta = uk(1);
    u.Mz    = uk(2);
    u.dk    = uk(3);

    return u;
}

State arrayToState(double *xk)
{
    State x;
    x.n     = xk[0];
    x.mu     = xk[1];
    x.vx    = xk[2];
    x.vy    = xk[3];
    x.r     = xk[4];
    x.Fm     = xk[5];
    x.delta = xk[6];
    x.k    = xk[7];

    return x;
}

Input arrayToInput(double *uk)
{
    Input u;
    u.dFm     = uk[0];
    u.dDelta = uk[1];
    u.Mz    = uk[2];
    u.dk    = uk[3];

    return u;
}

}