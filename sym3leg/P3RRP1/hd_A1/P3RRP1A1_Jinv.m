% Analytische Jacobi-Matrix für parallelen Roboter
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorpose und aktiven Gelenkkoordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% Jinv [3x3]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 18:10
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jinv = P3RRP1A1_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(4,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRP1A1_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRP1A1_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRP1A1_Jinv: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRP1A1_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRP1A1_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 18:09:49
% EndTime: 2018-12-20 18:09:52
% DurationCPUTime: 2.08s
% Computational Cost: add. (648->208), mult. (1356->362), div. (9->6), fcn. (459->20), ass. (0->149)
t106 = koppelP(1,1);
t103 = koppelP(1,2);
t140 = (qJ(3,1) * t103);
t127 = (pkin(2) * t140);
t100 = (qJ(3,1) ^ 2);
t107 = (pkin(2) ^ 2);
t84 = 1 + t107;
t128 = t100 - t84;
t175 = t128 * t106 + 2 * t127;
t139 = qJ(3,1) * t106;
t126 = pkin(2) * t139;
t29 = -t103 * t128 + 2 * t126;
t97 = xP(3);
t73 = sin(t97);
t74 = cos(t97);
t14 = -t175 * t73 + t29 * t74;
t176 = t175 * t74 + t73 * t29;
t90 = legFrame(1,3);
t69 = sin(t90);
t72 = cos(t90);
t177 = -t69 * t14 + t176 * t72;
t121 = 0.1e1 / 0.2e1 + t107 / 0.2e1;
t98 = qJ(3,3) ^ 2;
t52 = -t98 / 0.2e1 + t121;
t99 = qJ(3,2) ^ 2;
t53 = -t99 / 0.2e1 + t121;
t54 = -t100 / 0.2e1 + t121;
t101 = koppelP(3,2);
t104 = koppelP(3,1);
t135 = (qJ(3,3) * t104);
t122 = pkin(2) * t135;
t172 = t84 - t98;
t25 = t101 * t172 + (2 * t122);
t136 = qJ(3,3) * t101;
t123 = pkin(2) * t136;
t26 = t104 * t172 - 0.2e1 * t123;
t10 = t25 * t74 + t26 * t73;
t11 = -t73 * t25 + t26 * t74;
t88 = legFrame(3,3);
t67 = sin(t88);
t70 = cos(t88);
t174 = t67 * t10 + t11 * t70;
t102 = koppelP(2,2);
t105 = koppelP(2,1);
t137 = (qJ(3,2) * t105);
t124 = pkin(2) * t137;
t171 = t84 - t99;
t27 = t102 * t171 + (2 * t124);
t138 = qJ(3,2) * t102;
t125 = pkin(2) * t138;
t28 = t105 * t171 - 0.2e1 * t125;
t12 = t27 * t74 + t28 * t73;
t13 = -t73 * t27 + t28 * t74;
t89 = legFrame(2,3);
t68 = sin(t89);
t71 = cos(t89);
t173 = t68 * t12 + t13 * t71;
t132 = 2 * pkin(1);
t94 = cos(qJ(1,3));
t146 = qJ(3,3) * t94;
t91 = sin(qJ(1,3));
t147 = qJ(3,3) * t91;
t81 = qJ(1,3) + qJ(2,3);
t61 = sin(t81);
t64 = cos(t81);
t78 = t94 ^ 2;
t108 = pkin(1) ^ 2;
t85 = 2 + t108;
t170 = -(t61 * t146 - (pkin(1) * (t91 * t172 * t94 + (0.2e1 * t78 - 0.1e1) * qJ(3,3) * pkin(2)) * t61 + t147) * t64) * t132 - t98 * t85;
t95 = cos(qJ(1,2));
t149 = qJ(3,2) * t95;
t92 = sin(qJ(1,2));
t150 = qJ(3,2) * t92;
t82 = qJ(1,2) + qJ(2,2);
t62 = sin(t82);
t65 = cos(t82);
t79 = t95 ^ 2;
t169 = -(t62 * t149 - (pkin(1) * (t92 * t171 * t95 + (0.2e1 * t79 - 0.1e1) * qJ(3,2) * pkin(2)) * t62 + t150) * t65) * t132 - t99 * t85;
t96 = cos(qJ(1,1));
t152 = qJ(3,1) * t96;
t93 = sin(qJ(1,1));
t153 = qJ(3,1) * t93;
t83 = qJ(1,1) + qJ(2,1);
t63 = sin(t83);
t66 = cos(t83);
t80 = t96 ^ 2;
t168 = -(t63 * t152 - (pkin(1) * (-t93 * t128 * t96 + (0.2e1 * t80 - 0.1e1) * qJ(3,1) * pkin(2)) * t63 + t153) * t66) * t132 - (t100 * t85);
t167 = t61 * t64;
t166 = t62 * t65;
t165 = t63 * t66;
t154 = qJ(3,1) * t72;
t151 = qJ(3,2) * t71;
t148 = qJ(3,3) * t70;
t144 = t67 * qJ(3,3);
t143 = t68 * qJ(3,2);
t142 = t69 * qJ(3,1);
t133 = 0.1e1 / 0.4e1 + t107 / 0.4e1;
t131 = pkin(2) * t154;
t130 = pkin(2) * t151;
t129 = pkin(2) * t148;
t49 = pkin(2) * t144;
t120 = t172 * t70 + 0.2e1 * t49;
t50 = pkin(2) * t143;
t119 = t171 * t71 + 0.2e1 * t50;
t51 = pkin(2) * t142;
t118 = -t128 * t72 + 0.2e1 * t51;
t117 = pkin(2) * t93 * t152;
t116 = pkin(2) * t92 * t149;
t115 = pkin(2) * t91 * t146;
t111 = t128 * t80 + 0.2e1 * t117;
t110 = -t171 * t79 + 0.2e1 * t116;
t109 = -t172 * t78 + 0.2e1 * t115;
t57 = t66 ^ 2;
t56 = t65 ^ 2;
t55 = t64 ^ 2;
t48 = pkin(2) * t103 + t139;
t47 = pkin(2) * t106 - t140;
t46 = pkin(2) * t102 + t137;
t45 = pkin(2) * t105 - t138;
t44 = pkin(2) * t101 + t135;
t43 = pkin(2) * t104 - t136;
t42 = t107 * t103 + t103 + t126;
t41 = t107 * t106 + t106 - t127;
t40 = t107 * t102 + t102 + t124;
t39 = t107 * t105 + t105 - t125;
t38 = t107 * t101 + t101 + t122;
t37 = t107 * t104 + t104 - t123;
t36 = -t103 * t73 + t106 * t74;
t35 = -t102 * t73 + t105 * t74;
t34 = -t101 * t73 + t104 * t74;
t33 = t103 * t74 + t106 * t73;
t32 = t102 * t74 + t105 * t73;
t31 = t101 * t74 + t104 * t73;
t24 = -t128 * t69 - 0.2e1 * t131;
t23 = t171 * t68 - 0.2e1 * t130;
t22 = t172 * t67 - 0.2e1 * t129;
t21 = t24 * t96;
t20 = t23 * t95;
t19 = t22 * t94;
t9 = t118 * t96 + t24 * t93;
t8 = t119 * t95 + t23 * t92;
t7 = t120 * t94 + t22 * t91;
t6 = t14 * t72 + t176 * t69;
t5 = t12 * t71 - t13 * t68;
t4 = t10 * t70 - t11 * t67;
t3 = 0.1e1 / ((0.2e1 * (-t111 - t54) * t57 + t111) * t108 + t168);
t2 = 0.1e1 / ((0.2e1 * (-t110 - t53) * t56 + t110) * t108 + t169);
t1 = 0.1e1 / ((0.2e1 * (-t109 - t52) * t55 + t109) * t108 + t170);
t15 = [((t63 * t69 - t72 * t66) * qJ(3,1) + (-t21 * t57 - t9 * t165 + (t69 * t107 - t131 + t69) * t96 + (t118 * t57 - qJ(3,1) * (t69 * pkin(2) - t154)) * t93) * pkin(1)) * t3 ((-t72 * t63 - t69 * t66) * qJ(3,1) + (t9 * t57 - (t21 - 0.2e1 * (t54 * t72 + t51) * t93) * t165 - (t107 * t72 + t51 + t72) * t96 + (pkin(2) * t72 + t142) * t153) * pkin(1)) * t3 ((-(t33 * t72 - t69 * t36) * t66 + (t69 * t33 + t36 * t72) * t63) * qJ(3,1) + ((t177 * t96 + t6 * t93) * t57 - (-t177 * t93 + t6 * t96) * t165 - ((-t41 * t74 + t73 * t42) * t72 - (t73 * t41 + t42 * t74) * t69) * t96 - ((t47 * t74 - t73 * t48) * t72 + t69 * (t73 * t47 + t48 * t74)) * t153) * pkin(1)) / ((0.4e1 * (-t54 * t80 + t117 - t100 / 0.4e1 + t133) * t57 - t111) * t108 - t168); ((t62 * t68 - t71 * t65) * qJ(3,2) + (-t20 * t56 - t8 * t166 + (t68 * t107 - t130 + t68) * t95 + (t119 * t56 - qJ(3,2) * (t68 * pkin(2) - t151)) * t92) * pkin(1)) * t2 ((-t71 * t62 - t68 * t65) * qJ(3,2) + (t8 * t56 - (t20 - 0.2e1 * (t53 * t71 + t50) * t92) * t166 - (t107 * t71 + t50 + t71) * t95 + (pkin(2) * t71 + t143) * t150) * pkin(1)) * t2 ((-(t32 * t71 - t68 * t35) * t65 + (t68 * t32 + t35 * t71) * t62) * qJ(3,2) + ((-t173 * t95 + t5 * t92) * t56 - (t173 * t92 + t5 * t95) * t166 - ((-t39 * t74 + t73 * t40) * t71 - (t73 * t39 + t40 * t74) * t68) * t95 - ((t45 * t74 - t73 * t46) * t71 + t68 * (t73 * t45 + t46 * t74)) * t150) * pkin(1)) / ((0.4e1 * (-t53 * t79 + t116 - t99 / 0.4e1 + t133) * t56 - t110) * t108 - t169); ((t61 * t67 - t70 * t64) * qJ(3,3) + (-t19 * t55 - t7 * t167 + (t67 * t107 - t129 + t67) * t94 + (t120 * t55 - qJ(3,3) * (t67 * pkin(2) - t148)) * t91) * pkin(1)) * t1 ((-t70 * t61 - t67 * t64) * qJ(3,3) + (t7 * t55 - (t19 - 0.2e1 * (t52 * t70 + t49) * t91) * t167 - (t107 * t70 + t49 + t70) * t94 + (pkin(2) * t70 + t144) * t147) * pkin(1)) * t1 ((-(t31 * t70 - t67 * t34) * t64 + (t67 * t31 + t34 * t70) * t61) * qJ(3,3) + ((-t174 * t94 + t4 * t91) * t55 - (t174 * t91 + t4 * t94) * t167 - ((-t37 * t74 + t73 * t38) * t70 - (t73 * t37 + t38 * t74) * t67) * t94 - ((t43 * t74 - t73 * t44) * t70 + t67 * (t73 * t43 + t44 * t74)) * t147) * pkin(1)) / ((0.4e1 * (-t52 * t78 + t115 - t98 / 0.4e1 + t133) * t55 - t109) * t108 - t170);];
Jinv  = t15;
