% Calculate minimal parameter regressor of inverse dynamics forces for
% P3PRRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-09 21:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3PRRR1G1P1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G1P1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:14:52
% EndTime: 2020-03-09 21:14:54
% DurationCPUTime: 1.54s
% Computational Cost: add. (16642->186), mult. (11084->356), div. (1944->9), fcn. (11400->30), ass. (0->183)
t145 = pkin(7) + qJ(2,1);
t135 = qJ(3,1) + t145;
t123 = sin(t135);
t126 = cos(t135);
t129 = sin(t145);
t132 = cos(t145);
t75 = t123 * t132 - t129 * t126;
t223 = 0.1e1 / t75;
t148 = legFrame(1,3);
t138 = sin(t148);
t141 = cos(t148);
t83 = t141 * t123 + t138 * t126;
t200 = t223 * t83;
t84 = -t123 * t138 + t141 * t126;
t199 = t223 * t84;
t144 = pkin(7) + qJ(2,2);
t134 = qJ(3,2) + t144;
t122 = sin(t134);
t125 = cos(t134);
t128 = sin(t144);
t131 = cos(t144);
t74 = t122 * t131 - t128 * t125;
t224 = 0.1e1 / t74;
t147 = legFrame(2,3);
t137 = sin(t147);
t140 = cos(t147);
t81 = t140 * t122 + t137 * t125;
t202 = t224 * t81;
t82 = -t122 * t137 + t140 * t125;
t201 = t224 * t82;
t143 = pkin(7) + qJ(2,3);
t133 = qJ(3,3) + t143;
t121 = sin(t133);
t124 = cos(t133);
t127 = sin(t143);
t130 = cos(t143);
t73 = t121 * t130 - t127 * t124;
t225 = 0.1e1 / t73;
t146 = legFrame(3,3);
t136 = sin(t146);
t139 = cos(t146);
t79 = t139 * t121 + t136 * t124;
t204 = t225 * t79;
t80 = -t121 * t136 + t139 * t124;
t203 = t225 * t80;
t227 = -2 * pkin(2);
t226 = 2 * pkin(2);
t162 = 1 / pkin(2) ^ 2;
t165 = t121 * t127 + t124 * t130;
t157 = xDP(2);
t158 = xDP(1);
t85 = t136 * t157 + t139 * t158;
t86 = t136 * t158 - t139 * t157;
t40 = t86 * t121 - t85 * t124;
t219 = (-(-t127 * t86 + t85 * t130) * pkin(2) + t40 * pkin(3)) * t225;
t161 = 1 / pkin(2);
t160 = 0.1e1 / pkin(3);
t180 = t160 * t219;
t216 = t40 * t225;
t31 = (t180 - t216) * t161;
t186 = t31 * t219;
t150 = xDDP(1);
t191 = t150 * t161;
t149 = xDDP(2);
t192 = t149 * t161;
t68 = 0.1e1 / t73 ^ 2;
t215 = t40 * t68;
t19 = t192 * t204 + t191 * t203 + (t225 * t186 - (pkin(3) * t31 - t165 * t216) * t215) * t162;
t222 = pkin(2) * t19;
t164 = t122 * t128 + t125 * t131;
t87 = t137 * t157 + t140 * t158;
t88 = t137 * t158 - t140 * t157;
t41 = t88 * t122 - t87 * t125;
t218 = (-(-t128 * t88 + t87 * t131) * pkin(2) + t41 * pkin(3)) * t224;
t179 = t160 * t218;
t214 = t41 * t224;
t32 = (t179 - t214) * t161;
t185 = t32 * t218;
t70 = 0.1e1 / t74 ^ 2;
t213 = t41 * t70;
t20 = t192 * t202 + t191 * t201 + (t224 * t185 - (pkin(3) * t32 - t164 * t214) * t213) * t162;
t221 = pkin(2) * t20;
t163 = t123 * t129 + t126 * t132;
t89 = t138 * t157 + t141 * t158;
t90 = t138 * t158 - t141 * t157;
t42 = t90 * t123 - t89 * t126;
t217 = (-(-t129 * t90 + t89 * t132) * pkin(2) + t42 * pkin(3)) * t223;
t178 = t160 * t217;
t212 = t42 * t223;
t33 = (t178 - t212) * t161;
t184 = t33 * t217;
t72 = 0.1e1 / t75 ^ 2;
t211 = t42 * t72;
t21 = t192 * t200 + t191 * t199 + (t223 * t184 - (t33 * pkin(3) - t163 * t212) * t211) * t162;
t220 = pkin(2) * t21;
t49 = pkin(2) * (t127 * t139 + t136 * t130) + t79 * pkin(3);
t210 = t49 * t225;
t50 = pkin(2) * (t128 * t140 + t137 * t131) + t81 * pkin(3);
t209 = t50 * t224;
t51 = pkin(2) * (t129 * t141 + t138 * t132) + t83 * pkin(3);
t208 = t51 * t223;
t52 = -pkin(2) * (t127 * t136 - t139 * t130) + t80 * pkin(3);
t207 = t52 * t225;
t53 = -pkin(2) * (t128 * t137 - t140 * t131) + t82 * pkin(3);
t206 = t53 * t224;
t54 = -pkin(2) * (t129 * t138 - t141 * t132) + t84 * pkin(3);
t205 = t54 * t223;
t198 = t49 * t149;
t197 = t50 * t149;
t196 = t51 * t149;
t195 = t52 * t150;
t194 = t53 * t150;
t193 = t54 * t150;
t118 = t146 + t133;
t112 = sin(t118);
t115 = cos(t118);
t190 = g(1) * t115 + g(2) * t112;
t119 = t147 + t134;
t113 = sin(t119);
t116 = cos(t119);
t189 = g(1) * t116 + g(2) * t113;
t120 = t148 + t135;
t114 = sin(t120);
t117 = cos(t120);
t188 = g(1) * t117 + g(2) * t114;
t187 = 0.2e1 * pkin(3);
t183 = t161 * t40 ^ 2 * t68;
t182 = t161 * t41 ^ 2 * t70;
t181 = t161 * t42 ^ 2 * t72;
t177 = g(1) * t112 - g(2) * t115;
t176 = g(1) * t113 - g(2) * t116;
t175 = g(1) * t114 - g(2) * t117;
t159 = pkin(3) ^ 2;
t28 = (-t216 + t180 / 0.2e1) * t161;
t174 = (-t31 * t159 + (-t165 * t28 * t187 + t216) * pkin(2)) * t162 * t215;
t29 = (-t214 + t179 / 0.2e1) * t161;
t173 = (-t32 * t159 + (-t164 * t29 * t187 + t214) * pkin(2)) * t162 * t213;
t30 = (-t212 + t178 / 0.2e1) * t161;
t172 = (-t33 * t159 + (-t163 * t30 * t187 + t212) * pkin(2)) * t162 * t211;
t171 = (t165 * pkin(2) + pkin(3)) * t162 * t186;
t170 = (t164 * pkin(2) + pkin(3)) * t162 * t185;
t169 = (t163 * pkin(2) + pkin(3)) * t162 * t184;
t168 = t161 * t28 * t180;
t167 = t161 * t29 * t179;
t166 = t161 * t30 * t178;
t156 = cos(qJ(3,1));
t155 = cos(qJ(3,2));
t154 = cos(qJ(3,3));
t153 = sin(qJ(3,1));
t152 = sin(qJ(3,2));
t151 = sin(qJ(3,3));
t142 = (xDDP(3) - g(3));
t96 = t141 * g(1) + t138 * g(2);
t95 = t140 * g(1) + t137 * g(2);
t94 = t139 * g(1) + t136 * g(2);
t93 = t138 * g(1) - t141 * g(2);
t92 = t137 * g(1) - t140 * g(2);
t91 = t136 * g(1) - t139 * g(2);
t60 = -t93 * t129 + t96 * t132;
t59 = -t92 * t128 + t95 * t131;
t58 = -t91 * t127 + t94 * t130;
t57 = t96 * t129 + t93 * t132;
t56 = t95 * t128 + t92 * t131;
t55 = t94 * t127 + t91 * t130;
t18 = t152 * t182 + t155 * t221 + t176;
t17 = t151 * t183 + t154 * t222 + t177;
t16 = -t153 * t220 + t156 * t181 + t188;
t15 = -t152 * t221 + t155 * t182 + t189;
t14 = -t151 * t222 + t154 * t183 + t190;
t13 = t153 * t181 + t156 * t220 + t175;
t12 = (-t172 - (t169 + (t193 + t196) * t161) * t223) * t160 + t21;
t11 = (-t173 - (t170 + (t194 + t197) * t161) * t224) * t160 + t20;
t10 = (-t174 - (t171 + (t195 + t198) * t161) * t225) * t160 + t19;
t9 = (-t172 / 0.2e1 - (t169 / 0.2e1 + (t193 / 0.2e1 + t196 / 0.2e1) * t161) * t223) * t160 + t21;
t8 = (-t173 / 0.2e1 - (t170 / 0.2e1 + (t194 / 0.2e1 + t197 / 0.2e1) * t161) * t224) * t160 + t20;
t7 = (-t174 / 0.2e1 - (t171 / 0.2e1 + (t195 / 0.2e1 + t198 / 0.2e1) * t161) * t225) * t160 + t19;
t6 = (t153 * t9 + t156 * t166) * t227 + t188;
t5 = (t152 * t8 + t155 * t167) * t227 + t189;
t4 = (t151 * t7 + t154 * t168) * t227 + t190;
t3 = (-t153 * t166 + t9 * t156) * t226 + t175;
t2 = (-t152 * t167 + t8 * t155) * t226 + t176;
t1 = (-t151 * t168 + t7 * t154) * t226 + t177;
t22 = [0, (t19 * t203 + t21 * t199 + t20 * t201) * t161, (t57 * t199 + t56 * t201 + t55 * t203) * t161, (t60 * t199 + t59 * t201 + t58 * t203) * t161, (t10 * t203 + t11 * t201 + t12 * t199 + (-t10 * t207 - t11 * t206 - t12 * t205) * t160) * t161, (t1 * t203 + t2 * t201 + t3 * t199 + (-t13 * t205 - t17 * t207 - t18 * t206) * t160) * t161, (t4 * t203 + t5 * t201 + t6 * t199 + (-t14 * t207 - t15 * t206 - t16 * t205) * t160) * t161, t150 - g(1); 0, (t19 * t204 + t20 * t202 + t21 * t200) * t161, (t57 * t200 + t56 * t202 + t55 * t204) * t161, (t60 * t200 + t59 * t202 + t58 * t204) * t161, (t10 * t204 + t11 * t202 + t12 * t200 + (-t10 * t210 - t11 * t209 - t12 * t208) * t160) * t161, (t1 * t204 + t2 * t202 + t3 * t200 + (-t13 * t208 - t17 * t210 - t18 * t209) * t160) * t161, (t4 * t204 + t5 * t202 + t6 * t200 + (-t14 * t210 - t15 * t209 - t16 * t208) * t160) * t161, t149 - g(2); 3 * t142, 0, 0, 0, 0, 0, 0, t142;];
tauX_reg  = t22;
