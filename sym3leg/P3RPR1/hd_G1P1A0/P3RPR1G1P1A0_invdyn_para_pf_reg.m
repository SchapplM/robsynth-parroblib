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
% Datum: 2018-12-20 17:54
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX_reg = P3RPR1G1P1A0_invdyn_para_pf_reg(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_reg: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_reg: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_reg: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_reg: qJ has to be [2x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_reg: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_reg: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_reg: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_invdyn_para_pf_reg: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:54:36
% EndTime: 2018-12-20 17:54:38
% DurationCPUTime: 1.85s
% Computational Cost: add. (12813->280), mult. (19981->467), div. (648->9), fcn. (9609->14), ass. (0->200)
t145 = sin(qJ(1,3));
t148 = cos(qJ(1,3));
t139 = legFrame(3,3);
t127 = sin(t139);
t130 = cos(t139);
t95 = g(1) * t127 - g(2) * t130;
t98 = g(1) * t130 + g(2) * t127;
t58 = t98 * t145 + t148 * t95;
t146 = sin(qJ(1,2));
t149 = cos(qJ(1,2));
t140 = legFrame(2,3);
t128 = sin(t140);
t131 = cos(t140);
t96 = g(1) * t128 - g(2) * t131;
t99 = g(1) * t131 + g(2) * t128;
t59 = t99 * t146 + t149 * t96;
t141 = legFrame(1,3);
t129 = sin(t141);
t132 = cos(t141);
t100 = g(1) * t132 + g(2) * t129;
t147 = sin(qJ(1,1));
t150 = cos(qJ(1,1));
t97 = g(1) * t129 - g(2) * t132;
t60 = t100 * t147 + t150 * t97;
t63 = t100 * t150 - t97 * t147;
t62 = -t96 * t146 + t149 * t99;
t61 = -t95 * t145 + t148 * t98;
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
t25 = (t145 * t40 - t148 * t43) * t130 + (t145 * t43 + t148 * t40) * t127;
t89 = t135 * t172 + t136 * t169;
t70 = t89 * t151 - t153;
t92 = -t135 * t169 + t136 * t172;
t73 = t92 * t151 + t152;
t34 = (t145 * t73 - t148 * t70) * t130 + (t145 * t70 + t148 * t73) * t127;
t222 = t34 * t25;
t189 = t160 * t222;
t159 = 0.1e1 / qJ(2,3) ^ 2;
t208 = t159 * t34;
t137 = t151 ^ 2;
t142 = xDDP(3);
t143 = xDDP(2);
t46 = -t137 * t89 + t142 * t92 + t143;
t144 = xDDP(1);
t49 = -t137 * t92 - t142 * t89 + t144;
t78 = t127 * t148 + t130 * t145;
t79 = -t127 * t145 + t130 * t148;
t13 = -t189 + (t79 * t49 + t78 * t46 - (-t154 * t34 + t25) * t208) * t158;
t225 = pkin(1) * t13;
t161 = qJ(2,2) ^ 2;
t164 = t162 / t161;
t170 = koppelP(2,2);
t173 = koppelP(2,1);
t114 = -qJ(2,2) * t170 + t154 * t173;
t117 = qJ(2,2) * t173 + t154 * t170;
t41 = qJ(2,2) * t153 + t125 + (t114 * t136 - t117 * t135) * t151;
t44 = qJ(2,2) * t152 - t126 + (t114 * t135 + t117 * t136) * t151;
t26 = (t146 * t41 - t149 * t44) * t131 + (t146 * t44 + t149 * t41) * t128;
t90 = t135 * t173 + t136 * t170;
t71 = t90 * t151 - t153;
t93 = -t135 * t170 + t136 * t173;
t74 = t93 * t151 + t152;
t35 = (t146 * t74 - t149 * t71) * t131 + (t146 * t71 + t149 * t74) * t128;
t221 = t35 * t26;
t188 = t164 * t221;
t163 = 0.1e1 / qJ(2,2) ^ 2;
t203 = t163 * t35;
t47 = -t137 * t90 + t142 * t93 + t143;
t50 = -t137 * t93 - t142 * t90 + t144;
t80 = t128 * t149 + t131 * t146;
t81 = -t128 * t146 + t131 * t149;
t14 = -t188 + (t81 * t50 + t80 * t47 - (-t154 * t35 + t26) * t203) * t162;
t224 = pkin(1) * t14;
t165 = qJ(2,1) ^ 2;
t168 = t166 / t165;
t171 = koppelP(1,2);
t174 = koppelP(1,1);
t115 = -qJ(2,1) * t171 + t154 * t174;
t118 = qJ(2,1) * t174 + t154 * t171;
t42 = qJ(2,1) * t153 + t125 + (t115 * t136 - t118 * t135) * t151;
t45 = qJ(2,1) * t152 - t126 + (t115 * t135 + t118 * t136) * t151;
t27 = (t147 * t42 - t150 * t45) * t132 + (t147 * t45 + t150 * t42) * t129;
t91 = t135 * t174 + t136 * t171;
t72 = t91 * t151 - t153;
t94 = -t135 * t171 + t136 * t174;
t75 = t94 * t151 + t152;
t36 = (t147 * t75 - t150 * t72) * t132 + (t147 * t72 + t150 * t75) * t129;
t220 = t36 * t27;
t187 = t168 * t220;
t167 = 0.1e1 / qJ(2,1) ^ 2;
t198 = t167 * t36;
t48 = -t137 * t91 + t142 * t94 + t143;
t51 = -t137 * t94 - t142 * t91 + t144;
t82 = t129 * t150 + t132 * t147;
t83 = -t129 * t147 + t132 * t150;
t15 = -t187 + (t83 * t51 + t82 * t48 - (-t154 * t36 + t27) * t198) * t166;
t223 = pkin(1) * t15;
t186 = -pkin(2) ^ 2 + (-0.2e1 * pkin(2) - pkin(1)) * pkin(1);
t219 = ((-t157 + t186) * t34 + t25 * t154) * t158 * t208 + t154 * t189;
t218 = ((-t161 + t186) * t35 + t26 * t154) * t162 * t203 + t154 * t188;
t217 = ((-t165 + t186) * t36 + t27 * t154) * t166 * t198 + t154 * t187;
t101 = t145 * t172 - t148 * t169;
t102 = t145 * t169 + t148 * t172;
t37 = (t101 * t136 - t102 * t135) * t130 + (t101 * t135 + t102 * t136) * t127;
t211 = t158 * t37;
t210 = t158 * t78;
t209 = t158 * t79;
t31 = t34 ^ 2;
t207 = t160 * t31;
t103 = t146 * t173 - t149 * t170;
t104 = t146 * t170 + t149 * t173;
t38 = (t103 * t136 - t104 * t135) * t131 + (t103 * t135 + t104 * t136) * t128;
t206 = t162 * t38;
t205 = t162 * t80;
t204 = t162 * t81;
t32 = t35 ^ 2;
t202 = t164 * t32;
t105 = t147 * t174 - t150 * t171;
t106 = t147 * t171 + t150 * t174;
t39 = (t105 * t136 - t106 * t135) * t132 + (t105 * t135 + t106 * t136) * t129;
t201 = t166 * t39;
t200 = t166 * t82;
t199 = t166 * t83;
t33 = t36 ^ 2;
t197 = t168 * t33;
t192 = 0.2e1 * t222;
t191 = 0.2e1 * t221;
t190 = 0.2e1 * t220;
t185 = t219 + t225;
t184 = t218 + t224;
t183 = t217 + t223;
t107 = -qJ(2,3) * t148 + t145 * t154;
t110 = qJ(2,3) * t145 + t148 * t154;
t52 = t107 * t130 + t110 * t127;
t55 = -t107 * t127 + t110 * t130;
t182 = -t46 * t52 - t49 * t55;
t108 = -qJ(2,2) * t149 + t146 * t154;
t111 = qJ(2,2) * t146 + t149 * t154;
t53 = t108 * t131 + t111 * t128;
t56 = -t108 * t128 + t111 * t131;
t181 = -t47 * t53 - t50 * t56;
t109 = -qJ(2,1) * t150 + t147 * t154;
t112 = qJ(2,1) * t147 + t150 * t154;
t54 = t109 * t132 + t112 * t129;
t57 = -t109 * t129 + t112 * t132;
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
t77 = t133 * t135 + t134 * t136;
t76 = t133 * t136 - t134 * t135;
t69 = t115 * t147 - t118 * t150;
t68 = t114 * t146 - t117 * t149;
t67 = t113 * t145 - t116 * t148;
t66 = t115 * t150 + t118 * t147;
t65 = t114 * t149 + t117 * t146;
t64 = t113 * t148 + t116 * t145;
t30 = (-t135 * t66 + t136 * t69) * t132 + (t135 * t69 + t136 * t66) * t129;
t29 = (-t135 * t65 + t136 * t68) * t131 + (t135 * t68 + t136 * t65) * t128;
t28 = (-t135 * t64 + t136 * t67) * t130 + (t135 * t67 + t136 * t64) * t127;
t12 = 0.2e1 * qJ(2,1) * t15 + t167 * t190 - t63;
t11 = 0.2e1 * qJ(2,2) * t14 + t163 * t191 - t62;
t10 = 0.2e1 * qJ(2,3) * t13 + t159 * t192 - t61;
t9 = t180 * t166 + t217 + 0.2e1 * t223 + t60;
t8 = t181 * t162 + t218 + 0.2e1 * t224 + t59;
t7 = t182 * t158 + t219 + 0.2e1 * t225 + t58;
t6 = (-t180 - t33) * t166 - t183 - t60;
t5 = (-t181 - t32) * t162 - t184 - t59;
t4 = (-t182 - t31) * t158 - t185 - t58;
t3 = (t123 * t147 - t124 * t150) * t132 + (t123 * t150 + t124 * t147) * t129 + t15 * t165 + t183 * pkin(1) + (t180 * pkin(1) + t190) * t166;
t2 = (t121 * t146 - t122 * t149) * t131 + (t121 * t149 + t122 * t146) * t128 + t14 * t161 + t184 * pkin(1) + (t181 * pkin(1) + t191) * t162;
t1 = (t119 * t145 - t120 * t148) * t130 + (t119 * t148 + t120 * t145) * t127 + t13 * t157 + t185 * pkin(1) + (t182 * pkin(1) + t192) * t158;
t16 = [t13 * t209 + t14 * t204 + t15 * t199, t60 * t199 + t59 * t204 + t58 * t209, t63 * t199 + t62 * t204 + t61 * t209 (-t15 * t57 + t83 * t9) * t166 + (-t14 * t56 + t8 * t81) * t162 + (-t13 * t55 + t7 * t79) * t158, t10 * t209 + t11 * t204 + t12 * t199 - t57 * t197 - t56 * t202 - t55 * t207 (t3 * t83 + t57 * t6) * t166 + (t2 * t81 + t5 * t56) * t162 + (t1 * t79 + t4 * t55) * t158, 0, t88, -t87, -t135 * t76 + t136 * t77; t13 * t210 + t14 * t205 + t15 * t200, t60 * t200 + t59 * t205 + t58 * t210, t63 * t200 + t62 * t205 + t61 * t210 (-t15 * t54 + t82 * t9) * t166 + (-t14 * t53 + t8 * t80) * t162 + (-t13 * t52 + t7 * t78) * t158, t10 * t210 + t11 * t205 + t12 * t200 - t54 * t197 - t53 * t202 - t52 * t207 (t3 * t82 + t54 * t6) * t166 + (t2 * t80 + t5 * t53) * t162 + (t1 * t78 + t4 * t52) * t158, 0, t87, t88, t135 * t77 + t136 * t76; t13 * t211 + t14 * t206 + t15 * t201, t60 * t201 + t59 * t206 + t58 * t211, t63 * t201 + t62 * t206 + t61 * t211 (-t15 * t30 + t39 * t9) * t166 + (-t14 * t29 + t38 * t8) * t162 + (-t13 * t28 + t37 * t7) * t158, t10 * t211 + t11 * t206 + t12 * t201 - t30 * t197 - t29 * t202 - t28 * t207 (t3 * t39 + t30 * t6) * t166 + (t2 * t38 + t29 * t5) * t162 + (t1 * t37 + t28 * t4) * t158, t142, t76, -t77, 0;];
tauX_reg  = t16;
