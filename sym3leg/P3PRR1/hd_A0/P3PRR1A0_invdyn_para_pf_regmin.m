% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRR1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x8]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:47
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRR1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3PRR1A0_invdyn_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRR1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3PRR1A0_invdyn_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRR1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRR1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:47:38
% EndTime: 2019-05-03 14:47:39
% DurationCPUTime: 0.65s
% Computational Cost: add. (1236->124), mult. (2422->239), div. (387->8), fcn. (1860->14), ass. (0->121)
t88 = xP(3);
t64 = sin(t88);
t65 = cos(t88);
t91 = koppelP(1,2);
t94 = koppelP(1,1);
t106 = t64 * t91 - t65 * t94;
t49 = t64 * t94 + t65 * t91;
t75 = legFrame(1,3);
t58 = sin(t75);
t61 = cos(t75);
t138 = t106 * t58 + t49 * t61;
t90 = koppelP(2,2);
t93 = koppelP(2,1);
t107 = t64 * t90 - t65 * t93;
t48 = t64 * t93 + t65 * t90;
t74 = legFrame(2,3);
t57 = sin(t74);
t60 = cos(t74);
t137 = t107 * t57 + t48 * t60;
t89 = koppelP(3,2);
t92 = koppelP(3,1);
t108 = t64 * t89 - t65 * t92;
t47 = t64 * t92 + t65 * t89;
t73 = legFrame(3,3);
t56 = sin(t73);
t59 = cos(t73);
t136 = t108 * t56 + t47 * t59;
t79 = sin(qJ(2,3));
t66 = 0.1e1 / t79;
t80 = sin(qJ(2,2));
t68 = 0.1e1 / t80;
t81 = sin(qJ(2,1));
t70 = 0.1e1 / t81;
t135 = 0.1e1 / t79 ^ 2;
t134 = 0.1e1 / t80 ^ 2;
t133 = 0.1e1 / t81 ^ 2;
t82 = cos(qJ(2,3));
t132 = ((-t108 * t59 + t56 * t47) * t79 - t136 * t82) * t66;
t83 = cos(qJ(2,2));
t131 = ((-t107 * t60 + t57 * t48) * t80 - t137 * t83) * t68;
t84 = cos(qJ(2,1));
t130 = ((-t106 * t61 + t58 * t49) * t81 - t138 * t84) * t70;
t85 = xDP(3);
t86 = xDP(2);
t87 = xDP(1);
t25 = t136 * t85 - t56 * t86 - t87 * t59;
t22 = t25 ^ 2;
t96 = 0.1e1 / pkin(2) ^ 2;
t129 = t22 * t96;
t26 = t137 * t85 - t57 * t86 - t87 * t60;
t23 = t26 ^ 2;
t128 = t23 * t96;
t27 = t138 * t85 - t58 * t86 - t87 * t61;
t24 = t27 ^ 2;
t127 = t24 * t96;
t126 = t138 * t70;
t125 = t136 * t66;
t124 = t137 * t68;
t39 = -t56 * t79 + t59 * t82;
t123 = t39 * t66;
t40 = t56 * t82 + t79 * t59;
t122 = t40 * t66;
t41 = -t57 * t80 + t60 * t83;
t121 = t41 * t68;
t42 = t57 * t83 + t80 * t60;
t120 = t42 * t68;
t43 = -t58 * t81 + t61 * t84;
t119 = t43 * t70;
t44 = t58 * t84 + t81 * t61;
t118 = t44 * t70;
t117 = t56 * t66;
t116 = t57 * t68;
t115 = t58 * t70;
t114 = t59 * t66;
t113 = t60 * t68;
t112 = t61 * t70;
t111 = t82 * t129;
t110 = t83 * t128;
t109 = t84 * t127;
t95 = 0.1e1 / pkin(2);
t78 = xDDP(1);
t77 = xDDP(2);
t76 = xDDP(3);
t72 = t85 ^ 2;
t71 = t70 * t133;
t69 = t68 * t134;
t67 = t66 * t135;
t63 = t78 - g(1);
t62 = t77 - g(2);
t55 = t61 * g(1) + t58 * g(2);
t54 = t60 * g(1) + t57 * g(2);
t53 = t59 * g(1) + t56 * g(2);
t46 = -t64 * t76 - t65 * t72;
t45 = -t64 * t72 + t65 * t76;
t38 = t64 * t62 + t65 * t63;
t37 = t65 * t62 - t64 * t63;
t36 = t106 * t72 - t49 * t76 + t78;
t35 = t107 * t72 - t48 * t76 + t78;
t34 = t108 * t72 - t47 * t76 + t78;
t33 = -t106 * t76 - t72 * t49 + t77;
t32 = -t107 * t76 - t72 * t48 + t77;
t31 = -t108 * t76 - t72 * t47 + t77;
t18 = -t71 * t109 + (-t33 * t58 - t36 * t61) * t95 * t70;
t17 = -t69 * t110 + (-t32 * t57 - t35 * t60) * t95 * t68;
t16 = -t67 * t111 + (-t31 * t56 - t34 * t59) * t95 * t66;
t15 = t24 * t95 * t71 + t58 * g(1) - t61 * g(2) + (t33 * t44 + t36 * t43) * t70;
t14 = t95 * t23 * t69 + t57 * g(1) - t60 * g(2) + (t32 * t42 + t35 * t41) * t68;
t13 = t95 * t22 * t67 + t56 * g(1) - t59 * g(2) + (t31 * t40 + t34 * t39) * t66;
t12 = -t15 * t81 + t55 * t84;
t11 = -t14 * t80 + t54 * t83;
t10 = -t13 * t79 + t53 * t82;
t9 = t15 * t84 + t55 * t81;
t8 = t14 * t83 + t54 * t80;
t7 = t13 * t82 + t53 * t79;
t6 = -t70 * t127 + t84 * t18;
t5 = -t68 * t128 + t83 * t17;
t4 = -t66 * t129 + t82 * t16;
t3 = -t133 * t109 - t18 * t81;
t2 = -t134 * t110 - t17 * t80;
t1 = -t135 * t111 - t16 * t79;
t19 = [t15 * t119 + t14 * t121 + t13 * t123, (-t18 * t112 - t17 * t113 - t16 * t114) * t95, t4 * t123 + t5 * t121 + t6 * t119 + (-t9 * t112 - t8 * t113 - t7 * t114) * t95, t1 * t123 + t2 * t121 + t3 * t119 + (-t10 * t114 - t11 * t113 - t12 * t112) * t95, 0, t46, -t45, -t64 * t37 + t65 * t38; t15 * t118 + t14 * t120 + t13 * t122, (-t18 * t115 - t17 * t116 - t16 * t117) * t95, t4 * t122 + t5 * t120 + t6 * t118 + (-t9 * t115 - t8 * t116 - t7 * t117) * t95, t1 * t122 + t2 * t120 + t3 * t118 + (-t10 * t117 - t11 * t116 - t12 * t115) * t95, 0, t45, t46, t65 * t37 + t64 * t38; t13 * t132 + t15 * t130 + t14 * t131, (t17 * t124 + t16 * t125 + t18 * t126) * t95, t4 * t132 + t5 * t131 + t6 * t130 + (t8 * t124 + t7 * t125 + t9 * t126) * t95, t1 * t132 + t2 * t131 + t3 * t130 + (t10 * t125 + t11 * t124 + t12 * t126) * t95, t76, t37, -t38, 0;];
tauX_reg  = t19;
