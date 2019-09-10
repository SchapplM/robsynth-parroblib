% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RPR1G1P1A0
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
%   pkin=[a2,a3,d1,d3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x10]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPR1G1P1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_regmin: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:09
% EndTime: 2019-05-03 14:58:11
% DurationCPUTime: 1.86s
% Computational Cost: add. (12813->280), mult. (19981->467), div. (648->9), fcn. (9609->14), ass. (0->200)
t145 = sin(qJ(1,3));
t148 = cos(qJ(1,3));
t139 = legFrame(3,3);
t127 = sin(t139);
t130 = cos(t139);
t95 = t127 * g(1) - t130 * g(2);
t98 = t130 * g(1) + t127 * g(2);
t58 = t98 * t145 + t95 * t148;
t146 = sin(qJ(1,2));
t149 = cos(qJ(1,2));
t140 = legFrame(2,3);
t128 = sin(t140);
t131 = cos(t140);
t96 = t128 * g(1) - t131 * g(2);
t99 = t131 * g(1) + t128 * g(2);
t59 = t99 * t146 + t96 * t149;
t141 = legFrame(1,3);
t129 = sin(t141);
t132 = cos(t141);
t100 = t132 * g(1) + t129 * g(2);
t147 = sin(qJ(1,1));
t150 = cos(qJ(1,1));
t97 = t129 * g(1) - t132 * g(2);
t60 = t100 * t147 + t97 * t150;
t63 = t100 * t150 - t97 * t147;
t62 = -t96 * t146 + t99 * t149;
t61 = -t95 * t145 + t98 * t148;
t166 = 0.1e1 / qJ(2,1);
t162 = 0.1e1 / qJ(2,2);
t158 = 0.1e1 / qJ(2,3);
t226 = pkin(1) * g(2);
t154 = pkin(1) + pkin(2);
t157 = qJ(2,3) ^ 2;
t160 = t158 / t157;
t169 = koppelP(3,2);
t172 = koppelP(3,1);
t113 = -qJ(2,3) * t169 + t154 * t172;
t116 = qJ(2,3) * t172 + t154 * t169;
t152 = xDP(2);
t125 = t152 * t154;
t156 = xP(3);
t135 = sin(t156);
t136 = cos(t156);
t151 = xDP(3);
t153 = xDP(1);
t40 = qJ(2,3) * t153 + t125 + (t113 * t136 - t116 * t135) * t151;
t126 = t153 * t154;
t43 = qJ(2,3) * t152 - t126 + (t113 * t135 + t116 * t136) * t151;
t25 = (t40 * t145 - t43 * t148) * t130 + (t43 * t145 + t40 * t148) * t127;
t89 = t135 * t172 + t136 * t169;
t70 = t89 * t151 - t153;
t92 = -t135 * t169 + t136 * t172;
t73 = t92 * t151 + t152;
t34 = (t145 * t73 - t70 * t148) * t130 + t127 * (t145 * t70 + t73 * t148);
t222 = t34 * t25;
t189 = t160 * t222;
t159 = 0.1e1 / qJ(2,3) ^ 2;
t213 = t159 * t34;
t137 = t151 ^ 2;
t142 = xDDP(3);
t143 = xDDP(2);
t46 = -t137 * t89 + t142 * t92 + t143;
t144 = xDDP(1);
t49 = -t137 * t92 - t142 * t89 + t144;
t78 = t127 * t148 + t130 * t145;
t79 = -t127 * t145 + t130 * t148;
t13 = -t189 + (t79 * t49 + t78 * t46 - (-t154 * t34 + t25) * t213) * t158;
t225 = pkin(1) * t13;
t161 = qJ(2,2) ^ 2;
t164 = t162 / t161;
t170 = koppelP(2,2);
t173 = koppelP(2,1);
t114 = -qJ(2,2) * t170 + t154 * t173;
t117 = qJ(2,2) * t173 + t154 * t170;
t41 = qJ(2,2) * t153 + t125 + (t114 * t136 - t117 * t135) * t151;
t44 = qJ(2,2) * t152 - t126 + (t114 * t135 + t117 * t136) * t151;
t26 = (t41 * t146 - t44 * t149) * t131 + (t44 * t146 + t41 * t149) * t128;
t90 = t135 * t173 + t136 * t170;
t71 = t90 * t151 - t153;
t93 = -t135 * t170 + t136 * t173;
t74 = t93 * t151 + t152;
t35 = (t146 * t74 - t71 * t149) * t131 + t128 * (t146 * t71 + t74 * t149);
t221 = t35 * t26;
t188 = t164 * t221;
t163 = 0.1e1 / qJ(2,2) ^ 2;
t208 = t163 * t35;
t47 = -t137 * t90 + t142 * t93 + t143;
t50 = -t137 * t93 - t142 * t90 + t144;
t80 = t128 * t149 + t131 * t146;
t81 = -t128 * t146 + t131 * t149;
t14 = -t188 + (t81 * t50 + t80 * t47 - (-t154 * t35 + t26) * t208) * t162;
t224 = pkin(1) * t14;
t165 = qJ(2,1) ^ 2;
t168 = t166 / t165;
t171 = koppelP(1,2);
t174 = koppelP(1,1);
t115 = -qJ(2,1) * t171 + t154 * t174;
t118 = qJ(2,1) * t174 + t154 * t171;
t42 = qJ(2,1) * t153 + t125 + (t115 * t136 - t118 * t135) * t151;
t45 = qJ(2,1) * t152 - t126 + (t115 * t135 + t118 * t136) * t151;
t27 = (t42 * t147 - t45 * t150) * t132 + (t45 * t147 + t42 * t150) * t129;
t91 = t135 * t174 + t136 * t171;
t72 = t91 * t151 - t153;
t94 = -t135 * t171 + t136 * t174;
t75 = t94 * t151 + t152;
t36 = (t147 * t75 - t72 * t150) * t132 + t129 * (t147 * t72 + t75 * t150);
t220 = t36 * t27;
t187 = t168 * t220;
t167 = 0.1e1 / qJ(2,1) ^ 2;
t203 = t167 * t36;
t48 = -t137 * t91 + t142 * t94 + t143;
t51 = -t137 * t94 - t142 * t91 + t144;
t82 = t129 * t150 + t132 * t147;
t83 = -t129 * t147 + t132 * t150;
t15 = -t187 + (t83 * t51 + t82 * t48 - (-t154 * t36 + t27) * t203) * t166;
t223 = pkin(1) * t15;
t186 = -pkin(2) ^ 2 + (-0.2e1 * pkin(2) - pkin(1)) * pkin(1);
t219 = ((-t157 + t186) * t34 + t25 * t154) * t158 * t213 + t154 * t189;
t218 = ((-t161 + t186) * t35 + t26 * t154) * t162 * t208 + t154 * t188;
t217 = ((-t165 + t186) * t36 + t27 * t154) * t166 * t203 + t154 * t187;
t101 = t145 * t172 - t148 * t169;
t102 = t145 * t169 + t148 * t172;
t37 = (t101 * t136 - t135 * t102) * t130 + t127 * (t135 * t101 + t102 * t136);
t216 = t158 * t37;
t215 = t158 * t78;
t214 = t158 * t79;
t31 = t34 ^ 2;
t212 = t160 * t31;
t103 = t146 * t173 - t149 * t170;
t104 = t146 * t170 + t149 * t173;
t38 = (t103 * t136 - t135 * t104) * t131 + t128 * (t135 * t103 + t104 * t136);
t211 = t162 * t38;
t210 = t162 * t80;
t209 = t162 * t81;
t32 = t35 ^ 2;
t207 = t164 * t32;
t105 = t147 * t174 - t150 * t171;
t106 = t147 * t171 + t150 * t174;
t39 = (t105 * t136 - t135 * t106) * t132 + t129 * (t135 * t105 + t106 * t136);
t206 = t166 * t39;
t205 = t166 * t82;
t204 = t166 * t83;
t33 = t36 ^ 2;
t202 = t168 * t33;
t192 = 0.2e1 * t222;
t191 = 0.2e1 * t221;
t190 = 0.2e1 * t220;
t185 = t219 + t225;
t184 = t218 + t224;
t183 = t217 + t223;
t107 = -t148 * qJ(2,3) + t145 * t154;
t110 = t145 * qJ(2,3) + t154 * t148;
t52 = t107 * t130 + t127 * t110;
t55 = -t127 * t107 + t110 * t130;
t182 = -t46 * t52 - t49 * t55;
t108 = -t149 * qJ(2,2) + t146 * t154;
t111 = t146 * qJ(2,2) + t154 * t149;
t53 = t108 * t131 + t128 * t111;
t56 = -t128 * t108 + t111 * t131;
t181 = -t47 * t53 - t50 * t56;
t109 = -t150 * qJ(2,1) + t147 * t154;
t112 = t147 * qJ(2,1) + t154 * t150;
t54 = t109 * t132 + t129 * t112;
t57 = -t129 * t109 + t112 * t132;
t180 = -t48 * t54 - t51 * t57;
t155 = pkin(1) * g(1);
t134 = t144 - g(1);
t133 = t143 - g(2);
t124 = g(1) * qJ(2,1) + t226;
t123 = -g(2) * qJ(2,1) + t155;
t122 = g(1) * qJ(2,2) + t226;
t121 = -g(2) * qJ(2,2) + t155;
t120 = g(1) * qJ(2,3) + t226;
t119 = -g(2) * qJ(2,3) + t155;
t88 = -t135 * t142 - t136 * t137;
t87 = -t135 * t137 + t136 * t142;
t77 = t135 * t133 + t136 * t134;
t76 = t136 * t133 - t135 * t134;
t69 = t147 * t115 - t118 * t150;
t68 = t146 * t114 - t117 * t149;
t67 = t145 * t113 - t116 * t148;
t66 = t115 * t150 + t118 * t147;
t65 = t114 * t149 + t117 * t146;
t64 = t113 * t148 + t116 * t145;
t30 = (-t135 * t66 + t69 * t136) * t132 + (t69 * t135 + t66 * t136) * t129;
t29 = (-t135 * t65 + t68 * t136) * t131 + (t68 * t135 + t65 * t136) * t128;
t28 = (-t135 * t64 + t67 * t136) * t130 + (t67 * t135 + t64 * t136) * t127;
t12 = 0.2e1 * qJ(2,1) * t15 + t167 * t190 - t63;
t11 = 0.2e1 * qJ(2,2) * t14 + t163 * t191 - t62;
t10 = 0.2e1 * qJ(2,3) * t13 + t159 * t192 - t61;
t9 = t180 * t166 + t217 + 0.2e1 * t223 + t60;
t8 = t181 * t162 + t218 + 0.2e1 * t224 + t59;
t7 = t182 * t158 + t219 + 0.2e1 * t225 + t58;
t6 = (-t180 - t33) * t166 - t183 - t60;
t5 = (-t181 - t32) * t162 - t184 - t59;
t4 = (-t182 - t31) * t158 - t185 - t58;
t3 = (t147 * t123 - t124 * t150) * t132 + (t123 * t150 + t147 * t124) * t129 + t15 * t165 + t183 * pkin(1) + (t180 * pkin(1) + t190) * t166;
t2 = (t146 * t121 - t122 * t149) * t131 + (t121 * t149 + t146 * t122) * t128 + t14 * t161 + t184 * pkin(1) + (t181 * pkin(1) + t191) * t162;
t1 = (t145 * t119 - t120 * t148) * t130 + (t119 * t148 + t145 * t120) * t127 + t13 * t157 + t185 * pkin(1) + (t182 * pkin(1) + t192) * t158;
t16 = [t13 * t214 + t14 * t209 + t15 * t204, t60 * t204 + t59 * t209 + t58 * t214, t63 * t204 + t62 * t209 + t61 * t214, (-t15 * t57 + t83 * t9) * t166 + (-t14 * t56 + t8 * t81) * t162 + (-t13 * t55 + t7 * t79) * t158, t10 * t214 + t11 * t209 + t12 * t204 - t57 * t202 - t56 * t207 - t55 * t212, (t3 * t83 + t57 * t6) * t166 + (t2 * t81 + t5 * t56) * t162 + (t1 * t79 + t4 * t55) * t158, 0, t88, -t87, -t135 * t76 + t136 * t77; t13 * t215 + t14 * t210 + t15 * t205, t60 * t205 + t59 * t210 + t58 * t215, t63 * t205 + t62 * t210 + t61 * t215, (-t15 * t54 + t82 * t9) * t166 + (-t14 * t53 + t8 * t80) * t162 + (-t13 * t52 + t7 * t78) * t158, t10 * t215 + t11 * t210 + t12 * t205 - t54 * t202 - t53 * t207 - t52 * t212, (t3 * t82 + t54 * t6) * t166 + (t2 * t80 + t5 * t53) * t162 + (t1 * t78 + t4 * t52) * t158, 0, t87, t88, t135 * t77 + t136 * t76; t13 * t216 + t14 * t211 + t15 * t206, t60 * t206 + t59 * t211 + t58 * t216, t63 * t206 + t62 * t211 + t61 * t216, (-t15 * t30 + t39 * t9) * t166 + (-t14 * t29 + t38 * t8) * t162 + (-t13 * t28 + t37 * t7) * t158, t10 * t216 + t11 * t211 + t12 * t206 - t30 * t202 - t29 * t207 - t28 * t212, (t3 * t39 + t30 * t6) * t166 + (t2 * t38 + t29 * t5) * t162 + (t1 * t37 + t28 * t4) * t158, t142, t76, -t77, 0;];
tauX_reg  = t16;
