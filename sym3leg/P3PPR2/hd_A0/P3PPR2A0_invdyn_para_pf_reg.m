% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PPR2A0
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
% Datum: 2018-12-20 17:31
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX_reg = P3PPR2A0_invdyn_para_pf_reg(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_reg: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_reg: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_reg: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PPR2A0_invdyn_para_pf_reg: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_reg: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PPR2A0_invdyn_para_pf_reg: pkin has to be [3x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PPR2A0_invdyn_para_pf_reg: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PPR2A0_invdyn_para_pf_reg: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:31:34
% EndTime: 2018-12-20 17:31:34
% DurationCPUTime: 0.19s
% Computational Cost: add. (369->61), mult. (595->91), div. (0->0), fcn. (420->8), ass. (0->54)
t44 = xDDP(2) - g(2);
t45 = xDDP(1) - g(1);
t55 = xP(3);
t46 = sin(t55);
t47 = cos(t55);
t56 = koppelP(3,2);
t59 = koppelP(3,1);
t32 = t46 * t59 + t47 * t56;
t35 = -t46 * t56 + t47 * t59;
t48 = xDP(3) ^ 2;
t52 = xDDP(3);
t70 = -t32 * t52 - t48 * t35 + t45;
t57 = koppelP(2,2);
t60 = koppelP(2,1);
t33 = t46 * t60 + t47 * t57;
t36 = -t46 * t57 + t47 * t60;
t69 = -t33 * t52 - t48 * t36 + t45;
t58 = koppelP(1,2);
t61 = koppelP(1,1);
t34 = t46 * t61 + t47 * t58;
t37 = -t46 * t58 + t47 * t61;
t68 = -t34 * t52 - t48 * t37 + t45;
t67 = t48 * t32 - t35 * t52 - t44;
t66 = t48 * t33 - t36 * t52 - t44;
t65 = t48 * t34 - t37 * t52 - t44;
t49 = legFrame(3,3);
t38 = sin(t49);
t41 = cos(t49);
t13 = -t67 * t38 + t70 * t41;
t50 = legFrame(2,3);
t39 = sin(t50);
t42 = cos(t50);
t14 = -t66 * t39 + t69 * t42;
t51 = legFrame(1,3);
t40 = sin(t51);
t43 = cos(t51);
t15 = -t65 * t40 + t68 * t43;
t26 = t38 * t59 - t41 * t56;
t27 = t39 * t60 - t42 * t57;
t28 = t40 * t61 - t43 * t58;
t29 = t38 * t56 + t41 * t59;
t30 = t39 * t57 + t42 * t60;
t31 = t40 * t58 + t43 * t61;
t64 = (t28 * t47 - t31 * t46) * t15 + (t27 * t47 - t30 * t46) * t14 + (t26 * t47 - t29 * t46) * t13;
t63 = t38 * t13 + t39 * t14 + t40 * t15;
t62 = t41 * t13 + t42 * t14 + t43 * t15;
t25 = -t46 * t52 - t47 * t48;
t24 = -t46 * t48 + t47 * t52;
t23 = t46 * t44 + t47 * t45;
t22 = t47 * t44 - t46 * t45;
t12 = t68 * t40 + t65 * t43;
t11 = t69 * t39 + t66 * t42;
t10 = t70 * t38 + t67 * t41;
t1 = [t62, t38 * t10 + t39 * t11 + t40 * t12 + t62, 0, t25, -t24, -t46 * t22 + t47 * t23; t63, -t41 * t10 - t42 * t11 - t43 * t12 + t63, 0, t24, t25, t47 * t22 + t46 * t23; t64 (-t28 * t46 - t31 * t47) * t12 + (-t27 * t46 - t30 * t47) * t11 + (-t26 * t46 - t29 * t47) * t10 + t64, t52, t22, -t23, 0;];
tauX_reg  = t1;
