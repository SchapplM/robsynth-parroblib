% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRRR1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% xDDP [3x1]
%   Generalized platform accelerations
% qJ [3x3]
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x12]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRRR1G3A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:24
% EndTime: 2020-03-09 21:02:25
% DurationCPUTime: 1.65s
% Computational Cost: add. (1934->230), mult. (5658->491), div. (2718->17), fcn. (5772->18), ass. (0->193)
t105 = xDP(2);
t106 = xDP(1);
t104 = xDP(3);
t99 = cos(qJ(2,3));
t173 = t104 * t99;
t86 = legFrame(3,2);
t58 = sin(t86);
t61 = cos(t86);
t92 = sin(qJ(3,3));
t98 = cos(qJ(3,3));
t40 = -t92 * t173 + (t105 * t58 - t106 * t61) * t98;
t37 = t40 ^ 2;
t93 = sin(qJ(2,3));
t66 = 0.1e1 / t93 ^ 2;
t202 = t37 * t66;
t100 = cos(qJ(3,2));
t101 = cos(qJ(2,2));
t166 = t101 * t104;
t87 = legFrame(2,2);
t59 = sin(t87);
t62 = cos(t87);
t94 = sin(qJ(3,2));
t41 = -t94 * t166 + (t105 * t59 - t106 * t62) * t100;
t38 = t41 ^ 2;
t95 = sin(qJ(2,2));
t69 = 0.1e1 / t95 ^ 2;
t201 = t38 * t69;
t102 = cos(qJ(3,1));
t103 = cos(qJ(2,1));
t165 = t103 * t104;
t88 = legFrame(1,2);
t60 = sin(t88);
t63 = cos(t88);
t96 = sin(qJ(3,1));
t42 = -t96 * t165 + (t105 * t60 - t106 * t63) * t102;
t39 = t42 ^ 2;
t97 = sin(qJ(2,1));
t72 = 0.1e1 / t97 ^ 2;
t200 = t39 * t72;
t107 = 1 / pkin(2);
t65 = 0.1e1 / t93;
t68 = 0.1e1 / t95;
t71 = 0.1e1 / t97;
t73 = 0.1e1 / t98;
t77 = 0.1e1 / t100;
t81 = 0.1e1 / t102;
t108 = 0.1e1 / pkin(2) ^ 2;
t112 = t100 ^ 2;
t78 = 0.1e1 / t112;
t115 = t102 ^ 2;
t82 = 0.1e1 / t115;
t124 = t98 ^ 2;
t74 = 0.1e1 / t124;
t207 = 2 * t107;
t76 = 0.1e1 / t124 ^ 2;
t170 = t108 * t76;
t144 = t66 * t170;
t184 = t74 * t99;
t149 = t92 * t184;
t164 = t104 * t107;
t188 = t65 * t99;
t75 = t73 * t74;
t89 = xDDP(3);
t90 = xDDP(2);
t91 = xDDP(1);
t19 = -(-t104 * t92 * t93 + t40 * t188) * t40 * t144 + (-t89 * t149 + (t58 * t90 - t61 * t91 - (-t40 * t92 + t173) * t75 * t164) * t73) * t65 * t107;
t206 = t19 * t98;
t85 = t104 ^ 2;
t172 = t107 * t85;
t185 = t73 * t92;
t49 = -t58 * t99 + t61 * t93;
t50 = t58 * t93 + t61 * t99;
t28 = -t58 * g(1) - t61 * g(2) + (t89 * t185 + t49 * t90 + t50 * t91 + (t107 * t202 + t172) * t75) * t65;
t55 = g(1) * t61 - g(2) * t58;
t22 = t28 * t99 + t55 * t93;
t205 = t22 * t65;
t183 = t77 * t94;
t51 = -t101 * t59 + t62 * t95;
t52 = t101 * t62 + t59 * t95;
t79 = t77 * t78;
t29 = -t59 * g(1) - t62 * g(2) + (t89 * t183 + t51 * t90 + t52 * t91 + (t107 * t201 + t172) * t79) * t68;
t56 = g(1) * t62 - g(2) * t59;
t23 = t101 * t29 + t56 * t95;
t204 = t23 * t68;
t182 = t81 * t96;
t53 = -t103 * t60 + t63 * t97;
t54 = t103 * t63 + t60 * t97;
t83 = t81 * t82;
t30 = -t60 * g(1) - t63 * g(2) + (t89 * t182 + t53 * t90 + t54 * t91 + (t107 * t200 + t172) * t83) * t71;
t57 = g(1) * t63 - g(2) * t60;
t24 = t103 * t30 + t57 * t97;
t203 = t24 * t71;
t167 = t108 * t85;
t141 = t92 * t167;
t171 = t107 * t89;
t46 = t75 * t141 + t73 * t171;
t199 = t46 * t92;
t198 = t46 * t98;
t140 = t94 * t167;
t47 = t79 * t140 + t77 * t171;
t197 = t47 * t94;
t139 = t96 * t167;
t48 = t83 * t139 + t81 * t171;
t196 = t48 * t96;
t195 = t49 * t65;
t194 = t50 * t65;
t193 = t51 * t68;
t192 = t52 * t68;
t191 = t53 * t71;
t190 = t54 * t71;
t189 = t65 * t73;
t187 = t68 * t77;
t186 = t71 * t81;
t80 = 0.1e1 / t112 ^ 2;
t169 = t108 * t80;
t143 = t69 * t169;
t178 = t101 * t78;
t147 = t94 * t178;
t179 = t101 * t68;
t20 = -(-t104 * t94 * t95 + t41 * t179) * t41 * t143 + (-t89 * t147 + (t59 * t90 - t62 * t91 - (-t41 * t94 + t166) * t79 * t164) * t77) * t68 * t107;
t181 = t100 * t20;
t180 = t100 * t47;
t84 = 0.1e1 / t115 ^ 2;
t168 = t108 * t84;
t142 = t72 * t168;
t174 = t103 * t82;
t145 = t96 * t174;
t175 = t103 * t71;
t21 = -(-t104 * t96 * t97 + t42 * t175) * t42 * t142 + (-t89 * t145 + (t60 * t90 - t63 * t91 - (-t42 * t96 + t165) * t83 * t164) * t81) * t71 * t107;
t177 = t102 * t21;
t176 = t102 * t48;
t163 = t104 * t108;
t162 = t76 * t202;
t161 = t80 * t201;
t160 = t84 * t200;
t159 = t58 * t189;
t158 = t59 * t187;
t157 = t60 * t186;
t156 = t61 * t189;
t155 = t62 * t187;
t154 = t63 * t186;
t153 = t65 * t185;
t152 = t65 * t184;
t151 = t68 * t183;
t150 = t71 * t182;
t148 = t68 * t178;
t146 = t71 * t174;
t138 = t22 * t153;
t137 = t23 * t151;
t136 = t24 * t150;
t135 = t65 * t149;
t134 = t68 * t147;
t133 = t71 * t145;
t132 = t40 * t65 * t163;
t131 = t41 * t68 * t163;
t130 = t42 * t71 * t163;
t129 = 0.2e1 * t74 * t132;
t128 = 0.2e1 * t78 * t131;
t127 = 0.2e1 * t82 * t130;
t70 = t96 ^ 2;
t67 = t94 ^ 2;
t64 = t92 ^ 2;
t45 = t81 * t167 + t196;
t44 = t77 * t167 + t197;
t43 = t73 * t167 + t199;
t36 = -t82 * t139 + t176;
t35 = -t78 * t140 + t180;
t34 = -t74 * t141 + t198;
t33 = (t82 * t85 + t160) * t108;
t32 = (t78 * t85 + t161) * t108;
t31 = (t74 * t85 + t162) * t108;
t27 = t103 * t57 - t30 * t97;
t26 = t101 * t56 - t29 * t95;
t25 = -t28 * t93 + t55 * t99;
t18 = -t39 * t71 * t168 + t103 * t21;
t17 = -t38 * t68 * t169 + t101 * t20;
t16 = -t37 * t65 * t170 + t19 * t99;
t15 = -t103 * t39 * t142 - t21 * t97;
t14 = -t101 * t38 * t143 - t20 * t95;
t13 = -t37 * t99 * t144 - t19 * t93;
t12 = t96 * t127 + t21 * t70;
t11 = t94 * t128 + t20 * t67;
t10 = t92 * t129 + t19 * t64;
t9 = t96 * t177 + (0.2e1 * t81 - t83) * t130;
t8 = t94 * t181 + (0.2e1 * t77 - t79) * t131;
t7 = t92 * t206 + (0.2e1 * t73 - t75) * t132;
t6 = (t33 * t96 - t176) * t97 - t103 * (t21 * t96 + t127);
t5 = (t32 * t94 - t180) * t95 - t101 * (t20 * t94 + t128);
t4 = (t31 * t92 - t198) * t93 - t99 * (t19 * t92 + t129);
t3 = (-t102 * t33 - t196) * t97 + (-0.2e1 * t83 * t96 * t130 + t177) * t103;
t2 = (-t100 * t32 - t197) * t95 + (-0.2e1 * t79 * t94 * t131 + t181) * t101;
t1 = (-t31 * t98 - t199) * t93 + (-0.2e1 * t75 * t92 * t132 + t206) * t99;
t109 = [t30 * t190 + t29 * t192 + t28 * t194, (-t21 * t154 - t20 * t155 - t19 * t156) * t107, t16 * t194 + t17 * t192 + t18 * t190 + (-t24 * t154 - t23 * t155 - t22 * t156) * t107, t13 * t194 + t14 * t192 + t15 * t190 + (-t27 * t154 - t26 * t155 - t25 * t156) * t107, (-t10 * t156 - t11 * t155 - t12 * t154) * t107, (-t9 * t154 - t8 * t155 - t7 * t156) * t207, (-t45 * t154 - t44 * t155 - t43 * t156) * t107, (-t36 * t154 - t35 * t155 - t34 * t156) * t107, 0, t1 * t194 + t2 * t192 + t3 * t190 + (-t63 * t203 - t62 * t204 - t61 * t205) * t107, t4 * t194 + t5 * t192 + t6 * t190 + (t63 * t136 + t62 * t137 + t61 * t138) * t107, t91 - g(1); t30 * t191 + t29 * t193 + t28 * t195, (t21 * t157 + t20 * t158 + t19 * t159) * t107, t16 * t195 + t17 * t193 + t18 * t191 + (t24 * t157 + t23 * t158 + t22 * t159) * t107, t13 * t195 + t14 * t193 + t15 * t191 + (t27 * t157 + t26 * t158 + t25 * t159) * t107, (t10 * t159 + t11 * t158 + t12 * t157) * t107, (t9 * t157 + t8 * t158 + t7 * t159) * t207, (t45 * t157 + t44 * t158 + t43 * t159) * t107, (t36 * t157 + t35 * t158 + t34 * t159) * t107, 0, t1 * t195 + t2 * t193 + t3 * t191 + (t60 * t203 + t59 * t204 + t58 * t205) * t107, t4 * t195 + t5 * t193 + t6 * t191 + (-t60 * t136 - t59 * t137 - t58 * t138) * t107, t90 - g(2); t30 * t150 + t29 * t151 + t28 * t153, (-t21 * t133 - t20 * t134 - t19 * t135) * t107, t16 * t153 + t17 * t151 + t18 * t150 + (-t24 * t133 - t23 * t134 - t22 * t135) * t107, t13 * t153 + t14 * t151 + t15 * t150 + (-t27 * t133 - t26 * t134 - t25 * t135) * t107, ((-t96 * t160 - t94 * t161 - t92 * t162) * t108 - t10 * t135 - t11 * t134 - t12 * t133) * t107, (-0.2e1 * t8 * t134 - 0.2e1 * t9 * t133 - 0.2e1 * t7 * t135 + (t81 * (-0.2e1 * t82 + t84) * t200 + t77 * (-0.2e1 * t78 + t80) * t201 + t73 * (-0.2e1 * t74 + t76) * t202) * t108) * t107, ((-t45 * t146 + t21 * t81) * t96 + (-t44 * t148 + t20 * t77) * t94 + (-t43 * t152 + t19 * t73) * t92) * t107, (-t36 * t133 - t35 * t134 - t34 * t135 + t19 + t20 + t21) * t107, (t46 * t73 + t47 * t77 + t48 * t81) * t107, t1 * t153 + t2 * t151 + t3 * t150 + ((-g(3) * t102 + (-t24 * t175 + t27) * t96) * t81 + (-g(3) * t100 + (-t23 * t179 + t26) * t94) * t77 + (-g(3) * t98 + (-t22 * t188 + t25) * t92) * t73) * t107, t4 * t153 + t5 * t151 + t6 * t150 + (t70 * t24 * t146 + t81 * (g(3) * t96 + t102 * t27) + t67 * t23 * t148 + t77 * (g(3) * t94 + t100 * t26) + t64 * t22 * t152 + t73 * (g(3) * t92 + t25 * t98)) * t107, t89 - g(3);];
tauX_reg  = t109;
