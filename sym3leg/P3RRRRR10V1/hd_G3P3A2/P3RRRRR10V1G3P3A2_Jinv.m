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
% Datum: 2020-08-06 23:33
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jinv = P3RRRRR10V1G3P3A2_Jinv(xP, qJ, pkin, koppelP, ...
legFrame)


%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(6,1),zeros(3,3),zeros(3,3)}
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G3P3A2_Jinv: qJ has to be [3x3] (double)');
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G3P3A2_Jinv: xP has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G3P3A2_Jinv: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G3P3A2_Jinv: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G3P3A2_Jinv: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From Jinv_para_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 23:32:54
% EndTime: 2020-08-06 23:32:55
% DurationCPUTime: 1.15s
% Computational Cost: add. (276->165), mult. (792->335), div. (18->6), fcn. (699->26), ass. (0->134)
t54 = sin(qJ(1,3));
t61 = cos(qJ(3,3));
t38 = t61 ^ 2;
t149 = pkin(2) * t38;
t52 = sin(qJ(3,3));
t16 = t52 * pkin(5) + pkin(2);
t93 = -t16 + 0.2e1 * t149;
t157 = t54 * t93;
t57 = sin(qJ(1,2));
t64 = cos(qJ(3,2));
t41 = t64 ^ 2;
t148 = pkin(2) * t41;
t55 = sin(qJ(3,2));
t18 = t55 * pkin(5) + pkin(2);
t92 = -t18 + 0.2e1 * t148;
t156 = t57 * t92;
t60 = sin(qJ(1,1));
t67 = cos(qJ(3,1));
t44 = t67 ^ 2;
t147 = pkin(2) * t44;
t58 = sin(qJ(3,1));
t20 = t58 * pkin(5) + pkin(2);
t91 = -t20 + 0.2e1 * t147;
t155 = t60 * t91;
t59 = sin(qJ(2,1));
t140 = pkin(2) * t59;
t100 = t60 * t140;
t47 = sin(pkin(3));
t69 = cos(qJ(1,1));
t124 = t47 * t69;
t139 = pkin(6) * t47;
t68 = cos(qJ(2,1));
t88 = (t68 + 0.1e1) * (t68 - 0.1e1) * t139;
t154 = (-(pkin(1) * t124 + t60 * pkin(5)) * t59 + t69 * t88) * t58 - t100;
t56 = sin(qJ(2,2));
t142 = pkin(2) * t56;
t101 = t57 * t142;
t66 = cos(qJ(1,2));
t125 = t47 * t66;
t65 = cos(qJ(2,2));
t89 = (t65 + 0.1e1) * (t65 - 0.1e1) * t139;
t153 = (-(pkin(1) * t125 + t57 * pkin(5)) * t56 + t66 * t89) * t55 - t101;
t53 = sin(qJ(2,3));
t144 = pkin(2) * t53;
t102 = t54 * t144;
t63 = cos(qJ(1,3));
t126 = t47 * t63;
t62 = cos(qJ(2,3));
t90 = (t62 + 0.1e1) * (t62 - 0.1e1) * t139;
t152 = (-(pkin(1) * t126 + t54 * pkin(5)) * t53 + t63 * t90) * t52 - t102;
t151 = pkin(1) * t47;
t48 = cos(pkin(3));
t150 = pkin(1) * t48;
t146 = pkin(2) * t47;
t145 = pkin(2) * t52;
t143 = pkin(2) * t55;
t141 = pkin(2) * t58;
t117 = t61 * t62;
t132 = t62 * pkin(6);
t135 = t53 * pkin(6);
t138 = 0.1e1 / ((-t61 * t144 + t132) * t150 + (-pkin(5) * t135 + (pkin(1) * t52 - pkin(5) * t117) * pkin(2)) * t47) / t61;
t115 = t64 * t65;
t131 = t65 * pkin(6);
t134 = t56 * pkin(6);
t137 = 0.1e1 / ((-t64 * t142 + t131) * t150 + (-pkin(5) * t134 + (pkin(1) * t55 - pkin(5) * t115) * pkin(2)) * t47) / t64;
t113 = t67 * t68;
t130 = t68 * pkin(6);
t133 = t59 * pkin(6);
t136 = 0.1e1 / ((-t67 * t140 + t130) * t150 + (-pkin(5) * t133 + (pkin(1) * t58 - pkin(5) * t113) * pkin(2)) * t47) / t67;
t129 = t47 * t53;
t128 = t47 * t56;
t127 = t47 * t59;
t122 = t52 * t62;
t120 = t55 * t65;
t118 = t58 * t68;
t116 = t62 * t63;
t114 = t65 * t66;
t112 = t68 * t69;
t111 = t47 * t132;
t110 = t47 * t131;
t109 = t47 * t130;
t108 = pkin(6) * t122;
t107 = t54 * t132;
t106 = pkin(6) * t120;
t105 = t57 * t131;
t104 = pkin(6) * t118;
t103 = t60 * t130;
t99 = t53 * t122;
t98 = t56 * t120;
t97 = t59 * t118;
t17 = pkin(1) + t135;
t96 = t17 * t122;
t19 = pkin(1) + t134;
t95 = t19 * t120;
t21 = pkin(1) + t133;
t94 = t21 * t118;
t87 = t38 * t102;
t86 = t41 * t101;
t85 = t44 * t100;
t49 = legFrame(3,2);
t25 = sin(t49);
t84 = t25 * t107;
t28 = cos(t49);
t83 = t28 * t107;
t50 = legFrame(2,2);
t26 = sin(t50);
t82 = t26 * t105;
t29 = cos(t50);
t81 = t29 * t105;
t51 = legFrame(1,2);
t27 = sin(t51);
t80 = t27 * t103;
t30 = cos(t51);
t79 = t30 * t103;
t40 = t62 ^ 2;
t13 = (t40 - 0.2e1) * t145 - pkin(5);
t43 = t65 ^ 2;
t14 = (t43 - 0.2e1) * t143 - pkin(5);
t46 = t68 ^ 2;
t15 = (t46 - 0.2e1) * t141 - pkin(5);
t78 = t99 * t146;
t77 = t98 * t146;
t76 = t97 * t146;
t75 = t63 * t78;
t74 = t66 * t77;
t73 = t69 * t76;
t37 = t48 ^ 2;
t12 = -pkin(5) + (t46 - 0.1e1) * t141;
t11 = -pkin(5) + (t43 - 0.1e1) * t143;
t10 = -pkin(5) + (t40 - 0.1e1) * t145;
t6 = t15 * t60 * t47 - t21 * t69;
t5 = t14 * t57 * t47 - t19 * t66;
t4 = t13 * t54 * t47 - t17 * t63;
t1 = [(((-t15 * t27 + t79) * t67 + (-t27 * t104 - t30 * t155) * t59) * t37 + ((t30 * t112 + 0.2e1 * t27 * t127) * t147 + (-t27 * t109 - t6 * t30) * t67 - (t27 * t20 + t58 * t79) * t127) * t48 + t30 * t85 + (t12 * t27 - t30 * t73) * t67 + t154 * t30 + t27 * t94) * t136, (((-t30 * t15 - t80) * t67 + (-t30 * t104 + t27 * t155) * t59) * t37 + (-(t27 * t112 - 0.2e1 * t30 * t127) * t147 + (-t30 * t109 + t6 * t27) * t67 + (-t30 * t20 + t58 * t80) * t127) * t48 - t27 * t85 + (t30 * t12 + t27 * t73) * t67 - t154 * t27 + t30 * t94) * t136, ((-pkin(6) * t97 - t15 * t67) * t124 * t48 + ((pkin(6) * t113 - t91 * t59) * t37 + t44 * t140 - t59 * t20) * t69 + ((-t68 * t147 - t21 * t67) * t48 + t67 * t76 + (t151 * t59 - t88) * t58) * t60) * t136; (((-t14 * t26 + t81) * t64 + (-t26 * t106 - t29 * t156) * t56) * t37 + ((t29 * t114 + 0.2e1 * t26 * t128) * t148 + (-t26 * t110 - t5 * t29) * t64 - (t26 * t18 + t55 * t81) * t128) * t48 + t29 * t86 + (t11 * t26 - t29 * t74) * t64 + t153 * t29 + t26 * t95) * t137, (((-t29 * t14 - t82) * t64 + (-t29 * t106 + t26 * t156) * t56) * t37 + (-(t26 * t114 - 0.2e1 * t29 * t128) * t148 + (-t29 * t110 + t5 * t26) * t64 + (-t29 * t18 + t55 * t82) * t128) * t48 - t26 * t86 + (t29 * t11 + t26 * t74) * t64 - t153 * t26 + t29 * t95) * t137, ((-pkin(6) * t98 - t14 * t64) * t125 * t48 + ((pkin(6) * t115 - t92 * t56) * t37 + t41 * t142 - t56 * t18) * t66 + ((-t65 * t148 - t19 * t64) * t48 + t64 * t77 + (t151 * t56 - t89) * t55) * t57) * t137; (((-t13 * t25 + t83) * t61 + (-t25 * t108 - t28 * t157) * t53) * t37 + ((t28 * t116 + 0.2e1 * t25 * t129) * t149 + (-t25 * t111 - t4 * t28) * t61 - (t25 * t16 + t52 * t83) * t129) * t48 + t28 * t87 + (t10 * t25 - t28 * t75) * t61 + t152 * t28 + t25 * t96) * t138, (((-t28 * t13 - t84) * t61 + (-t28 * t108 + t25 * t157) * t53) * t37 + (-(t25 * t116 - 0.2e1 * t28 * t129) * t149 + (-t28 * t111 + t4 * t25) * t61 + (-t28 * t16 + t52 * t84) * t129) * t48 - t25 * t87 + (t28 * t10 + t25 * t75) * t61 - t152 * t25 + t28 * t96) * t138, ((-pkin(6) * t99 - t13 * t61) * t126 * t48 + ((pkin(6) * t117 - t93 * t53) * t37 + t38 * t144 - t53 * t16) * t63 + ((-t62 * t149 - t17 * t61) * t48 + t61 * t78 + (t151 * t53 - t90) * t52) * t54) * t138;];
Jinv  = t1;
