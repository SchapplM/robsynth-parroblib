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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
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
% Datum: 2020-08-06 22:54
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V1G2P2A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G2P2A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G2P2A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G2P2A2_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G2P2A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G2P2A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 22:53:38
% EndTime: 2020-08-06 22:53:39
% DurationCPUTime: 1.14s
% Computational Cost: add. (276->165), mult. (792->351), div. (18->6), fcn. (699->26), ass. (0->137)
t63 = cos(qJ(1,3));
t61 = cos(qJ(3,3));
t38 = t61 ^ 2;
t149 = pkin(2) * t38;
t52 = sin(qJ(3,3));
t16 = pkin(5) * t52 + pkin(2);
t91 = -t16 + 0.2e1 * t149;
t157 = t63 * t91;
t66 = cos(qJ(1,2));
t64 = cos(qJ(3,2));
t41 = t64 ^ 2;
t148 = pkin(2) * t41;
t55 = sin(qJ(3,2));
t18 = pkin(5) * t55 + pkin(2);
t90 = -t18 + 0.2e1 * t148;
t156 = t66 * t90;
t69 = cos(qJ(1,1));
t67 = cos(qJ(3,1));
t44 = t67 ^ 2;
t147 = pkin(2) * t44;
t58 = sin(qJ(3,1));
t20 = pkin(5) * t58 + pkin(2);
t89 = -t20 + 0.2e1 * t147;
t155 = t69 * t89;
t47 = sin(pkin(3));
t57 = sin(qJ(1,2));
t124 = t47 * t57;
t56 = sin(qJ(2,2));
t140 = pkin(6) * t47;
t65 = cos(qJ(2,2));
t87 = (t65 + 0.1e1) * (t65 - 0.1e1) * t140;
t129 = t66 * pkin(2);
t99 = t56 * t129;
t154 = (-t56 * (pkin(1) * t124 - pkin(5) * t66) + t57 * t87) * t55 + t99;
t60 = sin(qJ(1,1));
t122 = t47 * t60;
t59 = sin(qJ(2,1));
t68 = cos(qJ(2,1));
t86 = (t68 + 0.1e1) * (t68 - 0.1e1) * t140;
t128 = t69 * pkin(2);
t98 = t59 * t128;
t153 = (-t59 * (pkin(1) * t122 - pkin(5) * t69) + t60 * t86) * t58 + t98;
t130 = t63 * pkin(2);
t53 = sin(qJ(2,3));
t100 = t53 * t130;
t54 = sin(qJ(1,3));
t126 = t47 * t54;
t62 = cos(qJ(2,3));
t88 = (t62 + 0.1e1) * (t62 - 0.1e1) * t140;
t152 = (-t53 * (pkin(1) * t126 - pkin(5) * t63) + t54 * t88) * t52 + t100;
t151 = pkin(1) * t47;
t48 = cos(pkin(3));
t150 = pkin(1) * t48;
t146 = pkin(2) * t52;
t145 = pkin(2) * t53;
t144 = pkin(2) * t55;
t143 = pkin(2) * t56;
t142 = pkin(2) * t58;
t141 = pkin(2) * t59;
t139 = pkin(6) * t53;
t138 = pkin(6) * t56;
t137 = pkin(6) * t59;
t136 = pkin(6) * t62;
t135 = pkin(6) * t65;
t134 = pkin(6) * t68;
t112 = t61 * t62;
t133 = 0.1e1 / ((-t61 * t145 + t136) * t150 + (-pkin(5) * t139 + (pkin(1) * t52 - pkin(5) * t112) * pkin(2)) * t47) / t61;
t111 = t64 * t65;
t132 = 0.1e1 / ((-t64 * t143 + t135) * t150 + (-pkin(5) * t138 + (pkin(1) * t55 - pkin(5) * t111) * pkin(2)) * t47) / t64;
t110 = t67 * t68;
t131 = 0.1e1 / ((-t67 * t141 + t134) * t150 + (-pkin(5) * t137 + (pkin(1) * t58 - pkin(5) * t110) * pkin(2)) * t47) / t67;
t127 = t47 * t53;
t125 = t47 * t56;
t123 = t47 * t59;
t121 = t52 * t62;
t120 = t52 * t63;
t119 = t54 * t62;
t118 = t55 * t65;
t117 = t55 * t66;
t116 = t57 * t65;
t115 = t58 * t68;
t114 = t58 * t69;
t113 = t60 * t68;
t109 = t47 * t136;
t108 = t47 * t135;
t107 = t47 * t134;
t106 = pkin(6) * t121;
t105 = pkin(6) * t118;
t104 = pkin(6) * t115;
t103 = t63 * t136;
t102 = t66 * t135;
t101 = t69 * t134;
t97 = t52 * t127;
t96 = t55 * t125;
t95 = t58 * t123;
t17 = pkin(1) + t139;
t94 = t17 * t121;
t19 = pkin(1) + t138;
t93 = t19 * t118;
t21 = pkin(1) + t137;
t92 = t21 * t115;
t85 = t52 * t103;
t84 = t55 * t102;
t83 = t38 * t100;
t82 = t41 * t99;
t81 = t44 * t98;
t51 = legFrame(1,2);
t27 = sin(t51);
t80 = t27 * t104;
t30 = cos(t51);
t79 = t30 * t104;
t40 = t62 ^ 2;
t13 = (t40 - 0.2e1) * t146 - pkin(5);
t43 = t65 ^ 2;
t14 = (t43 - 0.2e1) * t144 - pkin(5);
t46 = t68 ^ 2;
t15 = (t46 - 0.2e1) * t142 - pkin(5);
t78 = pkin(2) * t62 * t97;
t77 = pkin(2) * t65 * t96;
t76 = pkin(2) * t68 * t95;
t75 = t54 * t78;
t74 = t57 * t77;
t73 = t60 * t76;
t50 = legFrame(2,2);
t49 = legFrame(3,2);
t37 = t48 ^ 2;
t29 = cos(t50);
t28 = cos(t49);
t26 = sin(t50);
t25 = sin(t49);
t12 = -pkin(5) + (t46 - 0.1e1) * t142;
t11 = -pkin(5) + (t43 - 0.1e1) * t144;
t10 = -pkin(5) + (t40 - 0.1e1) * t146;
t6 = t15 * t69 * t47 + t21 * t60;
t5 = t14 * t66 * t47 + t19 * t57;
t4 = t13 * t63 * t47 + t17 * t54;
t1 = [(((-t30 * t101 - t15 * t27) * t67 + (t155 * t30 - t80) * t59) * t37 + ((t30 * t113 + 0.2e1 * t27 * t123) * t147 + (-t27 * t107 + t30 * t6) * t67 + (-t27 * t20 + t69 * t79) * t123) * t48 - t30 * t81 + (t12 * t27 - t30 * t73) * t67 + t153 * t30 + t27 * t92) * t131, (((t27 * t101 - t15 * t30) * t67 + (-t155 * t27 - t79) * t59) * t37 + (-(t27 * t113 - 0.2e1 * t30 * t123) * t147 + (-t30 * t107 - t27 * t6) * t67 - (t30 * t20 + t69 * t80) * t123) * t48 + t27 * t81 + (t12 * t30 + t27 * t73) * t67 - t153 * t27 + t30 * t92) * t131, ((pkin(6) * t110 - t89 * t59) * t60 * t37 + ((-t15 * t122 + t21 * t69) * t67 + (-t60 * pkin(6) * t95 + t44 * t128) * t68) * t48 + t44 * t60 * t141 - t67 * t69 * t76 + (-t114 * t151 - t60 * t20) * t59 + t86 * t114) * t131; (((-t29 * t102 - t14 * t26) * t64 + (-t26 * t105 + t156 * t29) * t56) * t37 + ((t29 * t116 + 0.2e1 * t26 * t125) * t148 + (-t26 * t108 + t29 * t5) * t64 + (-t18 * t26 + t29 * t84) * t125) * t48 - t29 * t82 + (t11 * t26 - t29 * t74) * t64 + t154 * t29 + t26 * t93) * t132, (((t26 * t102 - t14 * t29) * t64 + (-t29 * t105 - t26 * t156) * t56) * t37 + (-(t26 * t116 - 0.2e1 * t29 * t125) * t148 + (-t29 * t108 - t26 * t5) * t64 - (t18 * t29 + t26 * t84) * t125) * t48 + t26 * t82 + (t11 * t29 + t26 * t74) * t64 - t154 * t26 + t29 * t93) * t132, ((pkin(6) * t111 - t90 * t56) * t57 * t37 + ((-t14 * t124 + t19 * t66) * t64 + (-t57 * pkin(6) * t96 + t41 * t129) * t65) * t48 + t41 * t57 * t143 - t64 * t66 * t77 + (-t117 * t151 - t18 * t57) * t56 + t87 * t117) * t132; (((-t28 * t103 - t13 * t25) * t61 + (-t25 * t106 + t157 * t28) * t53) * t37 + ((t28 * t119 + 0.2e1 * t25 * t127) * t149 + (-t25 * t109 + t28 * t4) * t61 + (-t16 * t25 + t28 * t85) * t127) * t48 - t28 * t83 + (t10 * t25 - t28 * t75) * t61 + t152 * t28 + t25 * t94) * t133, (((t25 * t103 - t13 * t28) * t61 + (-t28 * t106 - t157 * t25) * t53) * t37 + (-(t25 * t119 - 0.2e1 * t28 * t127) * t149 + (-t28 * t109 - t25 * t4) * t61 - (t16 * t28 + t25 * t85) * t127) * t48 + t25 * t83 + (t10 * t28 + t25 * t75) * t61 - t152 * t25 + t28 * t94) * t133, ((pkin(6) * t112 - t91 * t53) * t54 * t37 + ((-t13 * t126 + t17 * t63) * t61 + (-t54 * pkin(6) * t97 + t38 * t130) * t62) * t48 + t38 * t54 * t145 - t61 * t63 * t78 + (-t120 * t151 - t16 * t54) * t53 + t88 * t120) * t133;];
Jinv  = t1;
