% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PPR1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [2x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x6]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-02 15:24
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PPR1A0_invdyn_para_pf_reg(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_reg: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_reg: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_reg: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR1A0_invdyn_para_pf_reg: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_reg: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR1A0_invdyn_para_pf_reg: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR1A0_invdyn_para_pf_reg: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR1A0_invdyn_para_pf_reg: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-02 15:24:32
% EndTime: 2019-05-02 15:24:32
% DurationCPUTime: 0.17s
% Computational Cost: add. (369->61), mult. (595->91), div. (0->0), fcn. (420->8), ass. (0->54)
t41 = xDDP(2) - g(2);
t42 = xDDP(1) - g(1);
t52 = xP(3);
t43 = sin(t52);
t44 = cos(t52);
t53 = koppelP(3,2);
t56 = koppelP(3,1);
t29 = t43 * t56 + t44 * t53;
t32 = -t43 * t53 + t44 * t56;
t45 = xDP(3) ^ 2;
t49 = xDDP(3);
t67 = t29 * t49 + t45 * t32 - t42;
t54 = koppelP(2,2);
t57 = koppelP(2,1);
t30 = t43 * t57 + t44 * t54;
t33 = -t43 * t54 + t44 * t57;
t66 = t30 * t49 + t45 * t33 - t42;
t55 = koppelP(1,2);
t58 = koppelP(1,1);
t31 = t43 * t58 + t44 * t55;
t34 = -t43 * t55 + t44 * t58;
t65 = t31 * t49 + t45 * t34 - t42;
t64 = -t45 * t29 + t32 * t49 + t41;
t63 = -t45 * t30 + t33 * t49 + t41;
t62 = -t45 * t31 + t34 * t49 + t41;
t46 = legFrame(3,3);
t35 = sin(t46);
t38 = cos(t46);
t23 = t35 * t56 - t38 * t53;
t24 = t35 * t53 + t38 * t56;
t47 = legFrame(2,3);
t36 = sin(t47);
t39 = cos(t47);
t25 = t36 * t57 - t39 * t54;
t26 = t36 * t54 + t39 * t57;
t48 = legFrame(1,3);
t37 = sin(t48);
t40 = cos(t48);
t27 = t37 * t58 - t40 * t55;
t28 = t37 * t55 + t40 * t58;
t7 = t67 * t35 + t64 * t38;
t8 = t66 * t36 + t63 * t39;
t9 = t65 * t37 + t62 * t40;
t61 = (t43 * t27 + t28 * t44) * t9 + (t43 * t25 + t26 * t44) * t8 + (t43 * t23 + t24 * t44) * t7;
t60 = t38 * t7 + t39 * t8 + t40 * t9;
t59 = -t35 * t7 - t36 * t8 - t37 * t9;
t22 = -t43 * t49 - t44 * t45;
t21 = -t43 * t45 + t44 * t49;
t20 = t43 * t41 + t44 * t42;
t19 = t44 * t41 - t43 * t42;
t12 = t62 * t37 - t65 * t40;
t11 = t63 * t36 - t66 * t39;
t10 = t64 * t35 - t67 * t38;
t1 = [t59, t38 * t10 + t39 * t11 + t40 * t12 + t59, 0, t22, -t21, -t43 * t19 + t44 * t20; t60, t35 * t10 + t36 * t11 + t37 * t12 + t60, 0, t21, t22, t44 * t19 + t43 * t20; t61 (t27 * t44 - t28 * t43) * t12 + (t25 * t44 - t26 * t43) * t11 + (t23 * t44 - t24 * t43) * t10 + t61, t49, t19, -t20, 0;];
tauX_reg  = t1;
