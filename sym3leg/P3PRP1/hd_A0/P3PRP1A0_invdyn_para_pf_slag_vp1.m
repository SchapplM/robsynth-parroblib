% Calculate vector of inverse dynamics forces for parallel robot
% P3PRP1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% m [4x1]
%   mass of all robot links (including platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tauX [3x1]
%   forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)
%   in platform coordinates xP

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:35
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX = P3PRP1A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:34:45
% EndTime: 2018-12-20 17:34:49
% DurationCPUTime: 3.64s
% Computational Cost: add. (26231->420), mult. (49631->639), div. (2214->3), fcn. (23525->14), ass. (0->249)
t193 = (pkin(2) ^ 2);
t157 = 1 + t193;
t183 = (qJ(3,1) ^ 2);
t141 = -t183 + t157;
t172 = cos(qJ(2,1));
t155 = t172 ^ 2;
t169 = sin(qJ(2,1));
t231 = t169 * t172;
t220 = pkin(2) * t231;
t267 = 2 * qJ(3,1);
t111 = 0.1e1 / (t141 * t155 + t220 * t267 - t157 - t183);
t178 = xP(3);
t149 = sin(t178);
t150 = cos(t178);
t186 = koppelP(1,2);
t189 = koppelP(1,1);
t128 = -t149 * t186 + t150 * t189;
t175 = xDP(3);
t176 = xDP(2);
t114 = t128 * t175 + t176;
t163 = legFrame(1,3);
t144 = sin(t163);
t147 = cos(t163);
t245 = qJ(3,1) * t144;
t226 = pkin(2) * t245;
t204 = t147 * t183 - t226;
t244 = qJ(3,1) * t147;
t221 = pkin(2) * t244;
t205 = -t144 * t183 - t221;
t90 = t204 * t172 - t169 * (t144 - t205);
t249 = t90 * t114;
t125 = t149 * t189 + t150 * t186;
t177 = xDP(1);
t117 = -t125 * t175 + t177;
t87 = t205 * t172 + t169 * (-t147 - t204);
t252 = t87 * t117;
t274 = (t249 + t252) * t111;
t182 = (qJ(3,2) ^ 2);
t140 = -t182 + t157;
t171 = cos(qJ(2,2));
t154 = t171 ^ 2;
t168 = sin(qJ(2,2));
t232 = t168 * t171;
t218 = pkin(2) * t232;
t266 = 2 * qJ(3,2);
t110 = 0.1e1 / (t140 * t154 + t218 * t266 - t157 - t182);
t185 = koppelP(2,2);
t188 = koppelP(2,1);
t127 = -t149 * t185 + t150 * t188;
t113 = t127 * t175 + t176;
t162 = legFrame(2,3);
t143 = sin(t162);
t146 = cos(t162);
t241 = qJ(3,2) * t143;
t225 = pkin(2) * t241;
t201 = t146 * t182 - t225;
t240 = qJ(3,2) * t146;
t222 = pkin(2) * t240;
t202 = -t143 * t182 - t222;
t89 = t201 * t171 - t168 * (t143 - t202);
t250 = t89 * t113;
t124 = t149 * t188 + t150 * t185;
t116 = -t124 * t175 + t177;
t86 = t202 * t171 + t168 * (-t146 - t201);
t253 = t86 * t116;
t273 = (t250 + t253) * t110;
t181 = (qJ(3,3) ^ 2);
t139 = -t181 + t157;
t170 = cos(qJ(2,3));
t153 = t170 ^ 2;
t167 = sin(qJ(2,3));
t233 = t167 * t170;
t219 = pkin(2) * t233;
t265 = 2 * qJ(3,3);
t109 = 0.1e1 / (t139 * t153 + t219 * t265 - t157 - t181);
t184 = koppelP(3,2);
t187 = koppelP(3,1);
t126 = -t149 * t184 + t150 * t187;
t112 = t126 * t175 + t176;
t161 = legFrame(3,3);
t142 = sin(t161);
t145 = cos(t161);
t237 = qJ(3,3) * t142;
t224 = pkin(2) * t237;
t198 = t145 * t181 - t224;
t236 = qJ(3,3) * t145;
t223 = pkin(2) * t236;
t199 = -t142 * t181 - t223;
t88 = t198 * t170 - t167 * (t142 - t199);
t251 = t88 * t112;
t123 = t149 * t187 + t150 * t184;
t115 = -t123 * t175 + t177;
t85 = t199 * t170 + t167 * (-t145 - t198);
t254 = t85 * t115;
t272 = (t251 + t254) * t109;
t271 = -2 * m(3);
t270 = 0.2e1 * t272;
t269 = 0.2e1 * t273;
t268 = 0.2e1 * t274;
t264 = m(3) * t109;
t263 = m(3) * t110;
t262 = m(3) * t111;
t261 = m(3) * t170;
t260 = m(3) * t171;
t259 = m(3) * t172;
t173 = pkin(2) + rSges(3,1);
t258 = m(3) * t173;
t234 = qJ(3,3) * t170;
t227 = -0.2e1 * t234;
t97 = t142 * t227 + t167 * (pkin(2) * t142 - t236);
t98 = t145 * t227 + t167 * (pkin(2) * t145 + t237);
t67 = (t112 * t97 + t115 * t98) * t109;
t238 = qJ(3,2) * t171;
t228 = -0.2e1 * t238;
t100 = t146 * t228 + t168 * (pkin(2) * t146 + t241);
t99 = t143 * t228 + t168 * (pkin(2) * t143 - t240);
t68 = (t100 * t116 + t113 * t99) * t110;
t242 = qJ(3,1) * t172;
t229 = -0.2e1 * t242;
t101 = t144 * t229 + t169 * (pkin(2) * t144 - t244);
t102 = t147 * t229 + t169 * (pkin(2) * t147 + t245);
t69 = (t101 * t114 + t102 * t117) * t111;
t257 = t109 * t67;
t256 = t110 * t68;
t255 = t111 * t69;
t192 = pkin(2) * t193;
t58 = t181 * t67;
t214 = pkin(2) * t270 - t193 * t67 - t58;
t230 = -t193 - t157;
t10 = (t214 * t234 + ((t192 + (1 + t181) * pkin(2)) * t67 - t272 + t230 * t272) * t167) * t257;
t129 = -g(1) * t142 + g(2) * t145;
t248 = t10 + t129;
t59 = t182 * t68;
t213 = pkin(2) * t269 - t193 * t68 - t59;
t11 = (t213 * t238 + ((t192 + (1 + t182) * pkin(2)) * t68 - t273 + t230 * t273) * t168) * t256;
t130 = -g(1) * t143 + g(2) * t146;
t247 = t11 + t130;
t60 = t183 * t69;
t212 = pkin(2) * t268 - t193 * t69 - t60;
t12 = (t212 * t242 + ((t192 + (1 + t183) * pkin(2)) * t69 - t274 + t230 * t274) * t169) * t255;
t131 = -g(1) * t144 + g(2) * t147;
t246 = t12 + t131;
t243 = qJ(3,1) * t155;
t239 = qJ(3,2) * t154;
t235 = qJ(3,3) * t153;
t158 = rSges(3,3) + qJ(3,3);
t174 = (m(2) * rSges(2,2));
t217 = -m(3) * t158 + t174;
t159 = rSges(3,3) + qJ(3,2);
t216 = -m(3) * t159 + t174;
t160 = rSges(3,3) + qJ(3,1);
t215 = -m(3) * t160 + t174;
t211 = (rSges(3,3) ^ 2) + t193 + (0.2e1 * pkin(2) + rSges(3,1)) * rSges(3,1);
t210 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,2) + Icges(2,3);
t206 = -t141 * t144 + 0.2e1 * t221;
t203 = -t140 * t143 + 0.2e1 * t222;
t200 = -t139 * t142 + 0.2e1 * t223;
t180 = rSges(4,1);
t179 = rSges(4,2);
t166 = xDDP(1);
t165 = xDDP(2);
t164 = xDDP(3);
t156 = t175 ^ 2;
t151 = m(1) + m(2) + m(3);
t138 = (m(2) * rSges(2,1)) + t258;
t134 = g(1) * t147 + g(2) * t144;
t133 = g(1) * t146 + g(2) * t143;
t132 = g(1) * t145 + g(2) * t142;
t122 = -t149 * t179 + t150 * t180;
t121 = t149 * t180 + t150 * t179;
t120 = t141 * t147 + 0.2e1 * t226;
t119 = t140 * t146 + 0.2e1 * t225;
t118 = t139 * t145 + 0.2e1 * t224;
t108 = t138 * t172 - t169 * t215;
t107 = t138 * t171 - t168 * t216;
t106 = t138 * t170 - t167 * t217;
t105 = ((rSges(3,3) * t267) + t183 + t211) * m(3) + t210;
t104 = ((rSges(3,3) * t266) + t182 + t211) * m(3) + t210;
t103 = ((rSges(3,3) * t265) + t181 + t211) * m(3) + t210;
t96 = -t125 * t164 - t128 * t156 + t166;
t95 = -t124 * t164 - t127 * t156 + t166;
t94 = -t123 * t164 - t126 * t156 + t166;
t93 = -t125 * t156 + t128 * t164 + t165;
t92 = -t124 * t156 + t127 * t164 + t165;
t91 = -t123 * t156 + t126 * t164 + t165;
t84 = -t120 * t231 + t193 * t144 + t206 * t155 + t144 - t221;
t83 = t120 * t155 - t147 * t193 + t206 * t231 - t147 - t226;
t82 = -t119 * t232 + t193 * t143 + t203 * t154 + t143 - t222;
t81 = t119 * t154 - t146 * t193 + t203 * t232 - t146 - t225;
t80 = -t118 * t233 + t193 * t142 + t200 * t153 + t142 - t223;
t79 = t118 * t153 - t145 * t193 + t200 * t233 - t145 - t224;
t72 = (t101 * t128 - t102 * t125) * t111;
t71 = (-t100 * t124 + t127 * t99) * t110;
t70 = (-t123 * t98 + t126 * t97) * t109;
t66 = pkin(2) * t69;
t65 = pkin(2) * t68;
t64 = pkin(2) * t67;
t57 = (-t125 * t87 + t128 * t90) * t111;
t56 = (-t124 * t86 + t127 * t89) * t110;
t55 = (-t123 * t85 + t126 * t88) * t109;
t51 = (-t125 * t84 + t128 * t83) * t111;
t50 = (-t124 * t82 + t127 * t81) * t110;
t49 = (-t123 * t80 + t126 * t79) * t109;
t48 = (-t102 * t173 - t172 * t84 + t87) * t262;
t47 = (-t101 * t173 - t172 * t83 + t90) * t262;
t46 = (-t100 * t173 - t171 * t82 + t86) * t263;
t45 = (-t171 * t81 - t173 * t99 + t89) * t263;
t44 = (-t170 * t80 - t173 * t98 + t85) * t264;
t43 = (-t170 * t79 - t173 * t97 + t88) * t264;
t42 = (t102 * t108 + t151 * t84 - t87 * t259) * t111;
t41 = (t101 * t108 + t151 * t83 - t90 * t259) * t111;
t40 = (t100 * t107 + t151 * t82 - t86 * t260) * t110;
t39 = (t107 * t99 + t151 * t81 - t89 * t260) * t110;
t38 = (t106 * t98 + t151 * t80 - t85 * t261) * t109;
t37 = (t106 * t97 + t151 * t79 - t88 * t261) * t109;
t36 = (t102 * t105 + t108 * t84 - t87 * t258) * t111;
t35 = (t101 * t105 + t108 * t83 - t90 * t258) * t111;
t34 = (t100 * t104 + t107 * t82 - t86 * t258) * t110;
t33 = (t104 * t99 + t107 * t81 - t89 * t258) * t110;
t32 = (t103 * t98 + t106 * t80 - t85 * t258) * t109;
t31 = (t103 * t97 + t106 * t79 - t88 * t258) * t109;
t30 = t66 - t274;
t29 = t65 - t273;
t28 = t64 - t272;
t27 = (-t172 * t51 - t173 * t72 + t57) * m(3);
t26 = (-t171 * t50 - t173 * t71 + t56) * m(3);
t25 = (-t170 * t49 - t173 * t70 + t55) * m(3);
t24 = t108 * t72 + t151 * t51 - t57 * t259;
t23 = t107 * t71 + t151 * t50 - t56 * t260;
t22 = t106 * t70 + t151 * t49 - t55 * t261;
t21 = t105 * t72 + t108 * t51 - t57 * t258;
t20 = t104 * t71 + t107 * t50 - t56 * t258;
t19 = t103 * t70 + t106 * t49 - t55 * t258;
t18 = ((0.2e1 * (t66 + (-t252 / 0.2e1 - t249 / 0.2e1) * t111) * t243 - (pkin(2) * t30 - t60) * t231 - qJ(3,1) * t274) * t111 + (-qJ(3,1) + t220 - t243) * t111 * t274) * t69;
t17 = ((0.2e1 * (t65 + (-t253 / 0.2e1 - t250 / 0.2e1) * t110) * t239 - (pkin(2) * t29 - t59) * t232 - qJ(3,2) * t273) * t110 + (-qJ(3,2) + t218 - t239) * t110 * t273) * t68;
t16 = ((0.2e1 * (t64 + (-t254 / 0.2e1 - t251 / 0.2e1) * t109) * t235 - (pkin(2) * t28 - t58) * t233 - qJ(3,3) * t272) * t109 + (-qJ(3,3) + t219 - t235) * t109 * t272) * t67;
t15 = ((t30 - t274) * t231 + (-t155 * t69 - t212 + t69) * qJ(3,1)) * t255;
t14 = ((t29 - t273) * t232 + (-t154 * t68 - t213 + t68) * qJ(3,2)) * t256;
t13 = ((t28 - t272) * t233 + (-t153 * t67 - t214 + t67) * qJ(3,3)) * t257;
t9 = (-t69 ^ 2 * t160 - t134 * t169 + t246 * t172 + t173 * t18 - t15) * m(3);
t8 = (-t68 ^ 2 * t159 - t133 * t168 + t173 * t17 + t247 * t171 - t14) * m(3);
t7 = (-t67 ^ 2 * t158 - t132 * t167 + t173 * t16 + t248 * t170 - t13) * m(3);
t6 = -t105 * t18 - t108 * t12 + ((-rSges(2,1) * t131 + rSges(2,2) * t134) * t172 + t169 * (rSges(2,1) * t134 + rSges(2,2) * t131)) * m(2) + (t173 * t15 + t69 * t160 * t268 + (-t131 * t173 - t134 * t160) * t172 + t169 * (-t131 * t160 + t134 * t173)) * m(3);
t5 = -t104 * t17 - t107 * t11 + ((-rSges(2,1) * t130 + rSges(2,2) * t133) * t171 + t168 * (rSges(2,1) * t133 + rSges(2,2) * t130)) * m(2) + (t173 * t14 + t68 * t159 * t269 + (-t130 * t173 - t133 * t159) * t171 + t168 * (-t130 * t159 + t133 * t173)) * m(3);
t4 = -t106 * t10 - t103 * t16 + ((-rSges(2,1) * t129 + rSges(2,2) * t132) * t170 + t167 * (rSges(2,1) * t132 + rSges(2,2) * t129)) * m(2) + (t173 * t13 + t67 * t158 * t270 + (-t129 * t173 - t132 * t158) * t170 + t167 * (-t129 * t158 + t132 * t173)) * m(3);
t3 = -t108 * t18 + t15 * t259 - ((t138 * t69 + t271 * t274) * t169 + t172 * t215 * t69) * t69 - t246 * t151;
t2 = -t107 * t17 + t14 * t260 - ((t138 * t68 + t271 * t273) * t168 + t171 * t216 * t68) * t68 - t247 * t151;
t1 = -t106 * t16 + t13 * t261 - ((t138 * t67 + t271 * t272) * t167 + t170 * t217 * t67) * t67 - t248 * t151;
t52 = [(-t121 * t164 - t122 * t156 - g(1) + t166) * m(4) + ((t102 * t36 + t42 * t84 + t48 * t87) * t96 + (t101 * t36 + t42 * t83 + t48 * t90) * t93 + t84 * t3 + t102 * t6 + t87 * t9) * t111 + ((t100 * t34 + t40 * t82 + t46 * t86) * t95 + (t34 * t99 + t40 * t81 + t46 * t89) * t92 + t82 * t2 + t100 * t5 + t86 * t8) * t110 + ((t32 * t98 + t38 * t80 + t44 * t85) * t94 + (t32 * t97 + t38 * t79 + t44 * t88) * t91 + t80 * t1 + t98 * t4 + t85 * t7) * t109; (-t121 * t156 + t122 * t164 - g(2) + t165) * m(4) + ((t102 * t35 + t41 * t84 + t47 * t87) * t96 + (t101 * t35 + t41 * t83 + t47 * t90) * t93 + t83 * t3 + t101 * t6 + t90 * t9) * t111 + ((t100 * t33 + t39 * t82 + t45 * t86) * t95 + (t33 * t99 + t39 * t81 + t45 * t89) * t92 + t81 * t2 + t99 * t5 + t89 * t8) * t110 + ((t31 * t98 + t37 * t80 + t43 * t85) * t94 + (t31 * t97 + t37 * t79 + t43 * t88) * t91 + t79 * t1 + t97 * t4 + t88 * t7) * t109; Icges(4,3) * t164 + t49 * t1 + t50 * t2 + t51 * t3 + t70 * t4 + t71 * t5 + t55 * t7 + t56 * t8 + t57 * t9 + t72 * t6 + (-t121 * t166 + t122 * t165 + (t179 ^ 2 + t180 ^ 2) * t164 + (g(1) * t180 + g(2) * t179) * t149 + (g(1) * t179 - g(2) * t180) * t150) * m(4) + ((t102 * t21 + t24 * t84 + t27 * t87) * t96 + (t101 * t21 + t24 * t83 + t27 * t90) * t93) * t111 + ((t100 * t20 + t23 * t82 + t26 * t86) * t95 + (t20 * t99 + t23 * t81 + t26 * t89) * t92) * t110 + ((t19 * t98 + t22 * t80 + t25 * t85) * t94 + (t19 * t97 + t22 * t79 + t25 * t88) * t91) * t109;];
tauX  = t52;
