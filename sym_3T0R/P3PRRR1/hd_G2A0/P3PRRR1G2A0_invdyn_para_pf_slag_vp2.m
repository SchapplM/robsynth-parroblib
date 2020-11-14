% Calculate vector of inverse dynamics forces for parallel robot
% P3PRRR1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:18
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3PRRR1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:18:19
% EndTime: 2020-03-09 21:18:22
% DurationCPUTime: 2.23s
% Computational Cost: add. (15642->227), mult. (7887->435), div. (3366->5), fcn. (8976->24), ass. (0->205)
t140 = 0.1e1 / pkin(3);
t119 = pkin(7) + qJ(2,1);
t103 = sin(t119);
t109 = qJ(3,1) + t119;
t93 = sin(t109);
t72 = pkin(2) * t103 + pkin(3) * t93;
t203 = t140 * t72;
t129 = sin(qJ(3,1));
t233 = mrSges(3,2) * pkin(2);
t100 = t129 * t233;
t141 = pkin(2) ^ 2;
t149 = m(3) * t141 + Ifges(2,3) + Ifges(3,3);
t132 = cos(qJ(3,1));
t234 = mrSges(3,1) * pkin(2);
t182 = t132 * t234;
t69 = -0.2e1 * t100 + t149 + 0.2e1 * t182;
t81 = Ifges(3,3) - t100 + t182;
t240 = t81 * t203 - t69 * t93;
t118 = pkin(7) + qJ(2,2);
t102 = sin(t118);
t108 = qJ(3,2) + t118;
t92 = sin(t108);
t71 = pkin(2) * t102 + pkin(3) * t92;
t204 = t140 * t71;
t131 = cos(qJ(3,2));
t183 = t131 * t234;
t128 = sin(qJ(3,2));
t99 = t128 * t233;
t68 = t149 - 0.2e1 * t99 + 0.2e1 * t183;
t80 = Ifges(3,3) - t99 + t183;
t239 = t80 * t204 - t68 * t92;
t117 = pkin(7) + qJ(2,3);
t101 = sin(t117);
t107 = qJ(3,3) + t117;
t91 = sin(t107);
t70 = pkin(2) * t101 + pkin(3) * t91;
t205 = t140 * t70;
t130 = cos(qJ(3,3));
t184 = t130 * t234;
t127 = sin(qJ(3,3));
t98 = t127 * t233;
t67 = t149 - 0.2e1 * t98 + 0.2e1 * t184;
t79 = Ifges(3,3) - t98 + t184;
t238 = t79 * t205 - t67 * t91;
t104 = cos(t117);
t94 = cos(t107);
t237 = 0.1e1 / (t101 * t94 - t104 * t91);
t105 = cos(t118);
t95 = cos(t108);
t236 = 0.1e1 / (t102 * t95 - t105 * t92);
t106 = cos(t119);
t96 = cos(t109);
t235 = 0.1e1 / (t103 * t96 - t106 * t93);
t232 = pkin(2) * (mrSges(3,1) * t127 + mrSges(3,2) * t130);
t231 = pkin(2) * (mrSges(3,1) * t128 + mrSges(3,2) * t131);
t230 = pkin(2) * (mrSges(3,1) * t129 + mrSges(3,2) * t132);
t142 = 0.1e1 / pkin(2);
t185 = t140 * t142;
t156 = t237 * t185;
t121 = legFrame(3,2);
t113 = cos(t121);
t138 = xDP(1);
t193 = t113 * t138;
t110 = sin(t121);
t137 = xDP(2);
t196 = t110 * t137;
t136 = xDP(3);
t73 = pkin(2) * t104 + pkin(3) * t94;
t199 = t73 * t136;
t19 = (t199 + (t193 - t196) * t70) * t156;
t202 = t142 * t237;
t159 = t110 * t202;
t223 = t237 * t91;
t173 = t113 * t223;
t186 = t138 * t142;
t187 = t136 * t142;
t222 = t237 * t94;
t22 = t91 * t137 * t159 - t173 * t186 - t187 * t222;
t16 = t19 + t22;
t229 = t16 * t19;
t155 = t236 * t185;
t122 = legFrame(2,2);
t114 = cos(t122);
t191 = t114 * t138;
t111 = sin(t122);
t195 = t111 * t137;
t74 = pkin(2) * t105 + pkin(3) * t95;
t198 = t74 * t136;
t20 = (t198 + (t191 - t195) * t71) * t155;
t201 = t142 * t236;
t158 = t111 * t201;
t221 = t236 * t92;
t171 = t114 * t221;
t220 = t236 * t95;
t23 = t92 * t137 * t158 - t171 * t186 - t187 * t220;
t17 = t20 + t23;
t228 = t17 * t20;
t154 = t235 * t185;
t123 = legFrame(1,2);
t115 = cos(t123);
t189 = t115 * t138;
t112 = sin(t123);
t194 = t112 * t137;
t75 = pkin(2) * t106 + pkin(3) * t96;
t197 = t75 * t136;
t21 = (t197 + (t189 - t194) * t72) * t154;
t200 = t142 * t235;
t157 = t112 * t200;
t219 = t235 * t93;
t169 = t115 * t219;
t218 = t235 * t96;
t24 = t93 * t137 * t157 - t169 * t186 - t187 * t218;
t18 = t21 + t24;
t227 = t18 * t21;
t226 = t237 * t73;
t225 = t236 * t74;
t224 = t235 * t75;
t214 = t79 * t91;
t213 = t80 * t92;
t212 = t81 * t93;
t211 = t110 * t237;
t210 = t111 * t236;
t209 = t112 * t235;
t208 = t140 * t237;
t207 = t140 * t236;
t206 = t140 * t235;
t192 = t113 * t142;
t190 = t114 * t142;
t188 = t115 * t142;
t120 = m(1) + m(2) + m(3);
t181 = (t110 * t113 + t111 * t114 + t112 * t115) * t120;
t180 = t70 * t211;
t179 = t91 * t211;
t178 = t71 * t210;
t177 = t92 * t210;
t176 = t72 * t209;
t175 = t93 * t209;
t174 = t113 * t237 * t70;
t172 = t114 * t236 * t71;
t170 = t115 * t235 * t72;
t168 = t237 * t205;
t167 = t73 * t208;
t166 = t236 * t204;
t165 = t74 * t207;
t164 = t235 * t203;
t163 = t75 * t206;
t153 = 0.2e1 * pkin(2) * pkin(3);
t133 = mrSges(3,2) * g(3);
t135 = mrSges(3,1) * g(3);
t85 = g(1) * t113 - g(2) * t110;
t152 = -t94 * (mrSges(3,1) * t85 - t133) + (mrSges(3,2) * t85 + t135) * t91;
t86 = g(1) * t114 - g(2) * t111;
t151 = -t95 * (mrSges(3,1) * t86 - t133) + (mrSges(3,2) * t86 + t135) * t92;
t87 = g(1) * t115 - g(2) * t112;
t150 = -t96 * (mrSges(3,1) * t87 - t133) + (mrSges(3,2) * t87 + t135) * t93;
t148 = t101 * t91 + t104 * t94;
t147 = t102 * t92 + t105 * t95;
t146 = t103 * t93 + t106 * t96;
t145 = t148 * pkin(2);
t144 = t147 * pkin(2);
t143 = t146 * pkin(2);
t139 = pkin(3) ^ 2;
t134 = mrSges(2,2) * g(3);
t126 = xDDP(1);
t125 = xDDP(2);
t124 = xDDP(3);
t116 = m(3) * pkin(2) + mrSges(2,1);
t97 = t116 * g(3);
t84 = g(1) * t112 + g(2) * t115;
t83 = g(1) * t111 + g(2) * t114;
t82 = g(1) * t110 + g(2) * t113;
t42 = (Ifges(3,3) * t163 - t81 * t218) * t142;
t41 = (Ifges(3,3) * t165 - t80 * t220) * t142;
t40 = (Ifges(3,3) * t167 - t79 * t222) * t142;
t39 = (-Ifges(3,3) * t203 + t212) * t157;
t38 = (-Ifges(3,3) * t204 + t213) * t158;
t37 = (-Ifges(3,3) * t205 + t214) * t159;
t36 = (Ifges(3,3) * t164 - t212 * t235) * t188;
t35 = (Ifges(3,3) * t166 - t213 * t236) * t190;
t34 = (Ifges(3,3) * t168 - t214 * t237) * t192;
t33 = (t81 * t163 - t69 * t218) * t142;
t32 = (t80 * t165 - t68 * t220) * t142;
t31 = (t79 * t167 - t67 * t222) * t142;
t30 = t240 * t157;
t29 = t239 * t158;
t28 = t238 * t159;
t27 = t240 * t235 * t188;
t26 = t239 * t236 * t190;
t25 = t238 * t237 * t192;
t15 = (t197 / 0.2e1 + (t189 / 0.2e1 - t194 / 0.2e1) * t72) * t154 + t24;
t14 = (t198 / 0.2e1 + (t191 / 0.2e1 - t195 / 0.2e1) * t71) * t155 + t23;
t13 = (t199 / 0.2e1 + (t193 / 0.2e1 - t196 / 0.2e1) * t70) * t156 + t22;
t12 = (pkin(3) * t229 + (t16 * pkin(3) + t22 * t145) * t22) * t202;
t11 = (pkin(3) * t227 + (pkin(3) * t18 + t24 * t143) * t24) * t200;
t10 = (pkin(3) * t228 + (t17 * pkin(3) + t23 * t144) * t23) * t201;
t9 = ((-t146 * t15 * t153 - t139 * t18 - t141 * t24) * t24 * t206 - (pkin(3) + t143) * t235 * t227) * t142;
t8 = ((-t147 * t14 * t153 - t139 * t17 - t141 * t23) * t23 * t207 - (pkin(3) + t144) * t236 * t228) * t142;
t7 = ((-t148 * t13 * t153 - t139 * t16 - t141 * t22) * t22 * t208 - (pkin(3) + t145) * t237 * t229) * t142;
t6 = t24 ^ 2 * t230 - Ifges(3,3) * t9 - t11 * t81 + t150;
t5 = t23 ^ 2 * t231 - Ifges(3,3) * t8 - t10 * t80 + t151;
t4 = t22 ^ 2 * t232 - Ifges(3,3) * t7 - t12 * t79 + t152;
t3 = -t69 * t11 - t81 * t9 - 0.2e1 * t21 * t15 * t230 + (-t116 * t87 + t134) * t106 + t103 * (mrSges(2,2) * t87 + t97) + t150;
t2 = -t68 * t10 - t80 * t8 - 0.2e1 * t20 * t14 * t231 + (-t116 * t86 + t134) * t105 + t102 * (mrSges(2,2) * t86 + t97) + t151;
t1 = -t67 * t12 - t79 * t7 - 0.2e1 * t19 * t13 * t232 + (-t116 * t85 + t134) * t104 + t101 * (mrSges(2,2) * t85 + t97) + t152;
t43 = [(-g(1) + t126) * m(4) + t181 * t125 + (-t110 * t82 - t111 * t83 - t112 * t84 + (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) * t126) * t120 + ((t25 * t179 + t26 * t177 + t27 * t175 + (-t36 * t176 - t35 * t178 - t34 * t180) * t140) * t125 + (-t25 * t222 - t26 * t220 - t27 * t218 + (t36 * t224 + t35 * t225 + t34 * t226) * t140) * t124 + ((t36 * t164 - t27 * t219) * t126 - t3 * t219 + t6 * t164) * t115 + ((t35 * t166 - t26 * t221) * t126 - t2 * t221 + t5 * t166) * t114 + ((t34 * t168 - t25 * t223) * t126 - t1 * t223 + t4 * t168) * t113) * t142; (-g(2) + t125) * m(4) + t181 * t126 + (-t113 * t82 - t114 * t83 - t115 * t84 + (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) * t125) * t120 + ((t28 * t173 + t29 * t171 + t30 * t169 + (t39 * t170 + t38 * t172 + t37 * t174) * t140) * t126 + (t28 * t222 + t29 * t220 + t30 * t218 + (t39 * t224 + t38 * t225 + t37 * t226) * t140) * t124 + ((-t39 * t203 - t30 * t93) * t125 + t93 * t3 - t6 * t203) * t209 + ((-t38 * t204 - t29 * t92) * t125 + t92 * t2 - t5 * t204) * t210 + ((-t37 * t205 - t28 * t91) * t125 + t91 * t1 - t4 * t205) * t211) * t142; (-g(3) + t124) * m(4) + (-t1 * t222 - t2 * t220 - t3 * t218 + (t6 * t224 + t5 * t225 + t4 * t226) * t140 + (-t31 * t173 - t32 * t171 - t33 * t169 + (t42 * t170 + t41 * t172 + t40 * t174) * t140) * t126 + (t31 * t179 + t32 * t177 + t33 * t175 + (-t42 * t176 - t41 * t178 - t40 * t180) * t140) * t125 + (-t31 * t222 - t32 * t220 - t33 * t218 + (t42 * t224 + t41 * t225 + t40 * t226) * t140) * t124) * t142;];
tauX  = t43;
