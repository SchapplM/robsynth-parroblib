% Calculate vector of inverse dynamics forces for parallel robot
% P3PRP2A0
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
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauX = P3PRP2A0_invdyn_para_pf_slag_vp1(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_invdyn_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:38:47
% EndTime: 2018-12-20 17:38:50
% DurationCPUTime: 3.80s
% Computational Cost: add. (26231->420), mult. (49631->642), div. (2214->3), fcn. (23525->14), ass. (0->246)
t192 = (qJ(3,1) ^ 2);
t202 = (pkin(2) ^ 2);
t265 = -t202 - 1;
t150 = -t192 - t265;
t181 = cos(qJ(2,1));
t164 = t181 ^ 2;
t178 = sin(qJ(2,1));
t231 = t178 * t181;
t221 = pkin(2) * t231;
t268 = 2 * qJ(3,1);
t111 = 0.1e1 / (t150 * t164 + t221 * t268 - t192 + t265);
t187 = xP(3);
t158 = sin(t187);
t159 = cos(t187);
t195 = koppelP(1,2);
t198 = koppelP(1,1);
t131 = -t158 * t195 + t159 * t198;
t184 = xDP(3);
t185 = xDP(2);
t114 = t131 * t184 + t185;
t172 = legFrame(1,3);
t156 = cos(t172);
t153 = sin(t172);
t244 = qJ(3,1) * t156;
t226 = pkin(2) * t244;
t209 = t153 * t192 - t226;
t245 = qJ(3,1) * t153;
t144 = pkin(2) * t245;
t228 = t156 * t192 + t144;
t90 = t209 * t181 - t178 * (t156 + t228);
t249 = t90 * t114;
t128 = t158 * t198 + t159 * t195;
t186 = xDP(1);
t117 = -t128 * t184 + t186;
t87 = t228 * t181 - t178 * (-t153 - t209);
t252 = t87 * t117;
t275 = (t249 + t252) * t111;
t191 = (qJ(3,2) ^ 2);
t149 = -t191 - t265;
t180 = cos(qJ(2,2));
t163 = t180 ^ 2;
t177 = sin(qJ(2,2));
t232 = t177 * t180;
t222 = pkin(2) * t232;
t267 = 2 * qJ(3,2);
t110 = 0.1e1 / (t149 * t163 + t222 * t267 - t191 + t265);
t194 = koppelP(2,2);
t197 = koppelP(2,1);
t130 = -t158 * t194 + t159 * t197;
t113 = t130 * t184 + t185;
t171 = legFrame(2,3);
t155 = cos(t171);
t152 = sin(t171);
t240 = qJ(3,2) * t155;
t225 = pkin(2) * t240;
t208 = t152 * t191 - t225;
t241 = qJ(3,2) * t152;
t143 = pkin(2) * t241;
t229 = t155 * t191 + t143;
t89 = t208 * t180 - t177 * (t155 + t229);
t250 = t89 * t113;
t127 = t158 * t197 + t159 * t194;
t116 = -t127 * t184 + t186;
t86 = t229 * t180 - t177 * (-t152 - t208);
t253 = t86 * t116;
t274 = (t250 + t253) * t110;
t190 = (qJ(3,3) ^ 2);
t148 = -t190 - t265;
t179 = cos(qJ(2,3));
t162 = t179 ^ 2;
t176 = sin(qJ(2,3));
t233 = t176 * t179;
t223 = pkin(2) * t233;
t266 = 2 * qJ(3,3);
t109 = 0.1e1 / (t148 * t162 + t223 * t266 - t190 + t265);
t193 = koppelP(3,2);
t196 = koppelP(3,1);
t129 = -t158 * t193 + t159 * t196;
t112 = t129 * t184 + t185;
t170 = legFrame(3,3);
t154 = cos(t170);
t151 = sin(t170);
t236 = qJ(3,3) * t154;
t224 = pkin(2) * t236;
t207 = t151 * t190 - t224;
t237 = qJ(3,3) * t151;
t142 = pkin(2) * t237;
t230 = t154 * t190 + t142;
t88 = t207 * t179 - t176 * (t154 + t230);
t251 = t88 * t112;
t126 = t158 * t196 + t159 * t193;
t115 = -t126 * t184 + t186;
t85 = t230 * t179 - t176 * (-t151 - t207);
t254 = t85 * t115;
t273 = (t251 + t254) * t109;
t272 = -2 * m(3);
t271 = 0.2e1 * t273;
t270 = 0.2e1 * t274;
t269 = 0.2e1 * t275;
t264 = m(3) * t109;
t263 = m(3) * t110;
t262 = m(3) * t111;
t261 = m(3) * t179;
t260 = m(3) * t180;
t259 = m(3) * t181;
t182 = pkin(2) + rSges(3,1);
t258 = m(3) * t182;
t234 = qJ(3,3) * t179;
t97 = 0.2e1 * t151 * t234 - t176 * (pkin(2) * t151 + t236);
t98 = -0.2e1 * t154 * t234 + t176 * (pkin(2) * t154 - t237);
t67 = (t112 * t98 + t115 * t97) * t109;
t238 = qJ(3,2) * t180;
t100 = -0.2e1 * t155 * t238 + t177 * (pkin(2) * t155 - t241);
t99 = 0.2e1 * t152 * t238 - t177 * (pkin(2) * t152 + t240);
t68 = (t100 * t113 + t116 * t99) * t110;
t242 = qJ(3,1) * t181;
t101 = 0.2e1 * t153 * t242 - t178 * (pkin(2) * t153 + t244);
t102 = -0.2e1 * t156 * t242 + t178 * (pkin(2) * t156 - t245);
t69 = (t101 * t117 + t102 * t114) * t111;
t257 = t109 * t67;
t256 = t110 * t68;
t255 = t111 * t69;
t201 = pkin(2) * t202;
t58 = t190 * t67;
t217 = pkin(2) * t271 - t202 * t67 - t58;
t227 = -t202 + t265;
t10 = (t217 * t234 + ((t201 + (1 + t190) * pkin(2)) * t67 - t273 + t227 * t273) * t176) * t257;
t135 = g(1) * t154 + g(2) * t151;
t248 = t10 + t135;
t59 = t191 * t68;
t216 = pkin(2) * t270 - t202 * t68 - t59;
t11 = (t216 * t238 + ((t201 + (1 + t191) * pkin(2)) * t68 - t274 + t227 * t274) * t177) * t256;
t136 = g(1) * t155 + g(2) * t152;
t247 = t11 + t136;
t60 = t192 * t69;
t215 = pkin(2) * t269 - t202 * t69 - t60;
t12 = (t215 * t242 + ((t201 + (1 + t192) * pkin(2)) * t69 - t275 + t227 * t275) * t178) * t255;
t137 = g(1) * t156 + g(2) * t153;
t246 = t12 + t137;
t243 = qJ(3,1) * t164;
t239 = qJ(3,2) * t163;
t235 = qJ(3,3) * t162;
t167 = rSges(3,3) + qJ(3,3);
t183 = (m(2) * rSges(2,2));
t220 = -m(3) * t167 + t183;
t168 = rSges(3,3) + qJ(3,2);
t219 = -m(3) * t168 + t183;
t169 = rSges(3,3) + qJ(3,1);
t218 = -m(3) * t169 + t183;
t214 = (rSges(3,3) ^ 2) + t202 + (0.2e1 * pkin(2) + rSges(3,1)) * rSges(3,1);
t213 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,2) + Icges(2,3);
t189 = rSges(4,1);
t188 = rSges(4,2);
t175 = xDDP(1);
t174 = xDDP(2);
t173 = xDDP(3);
t165 = t184 ^ 2;
t160 = m(1) + m(2) + m(3);
t141 = (m(2) * rSges(2,1)) + t258;
t134 = -g(1) * t153 + g(2) * t156;
t133 = -g(1) * t152 + g(2) * t155;
t132 = -g(1) * t151 + g(2) * t154;
t125 = -t158 * t188 + t159 * t189;
t124 = t158 * t189 + t159 * t188;
t123 = t150 * t153 + 0.2e1 * t226;
t122 = t149 * t152 + 0.2e1 * t225;
t121 = t148 * t151 + 0.2e1 * t224;
t120 = t150 * t156 - 0.2e1 * t144;
t119 = t149 * t155 - 0.2e1 * t143;
t118 = t148 * t154 - 0.2e1 * t142;
t108 = t141 * t181 - t178 * t218;
t107 = t141 * t180 - t177 * t219;
t106 = t141 * t179 - t176 * t220;
t105 = ((rSges(3,3) * t268) + t192 + t214) * m(3) + t213;
t104 = ((rSges(3,3) * t267) + t191 + t214) * m(3) + t213;
t103 = ((rSges(3,3) * t266) + t190 + t214) * m(3) + t213;
t96 = -t128 * t173 - t131 * t165 + t175;
t95 = -t127 * t173 - t130 * t165 + t175;
t94 = -t126 * t173 - t129 * t165 + t175;
t93 = -t128 * t165 + t131 * t173 + t174;
t92 = -t127 * t165 + t130 * t173 + t174;
t91 = -t126 * t165 + t129 * t173 + t174;
t84 = -t120 * t231 + t123 * t164 - t153 * t202 - t153 - t226;
t83 = -t119 * t232 + t122 * t163 - t152 * t202 - t152 - t225;
t82 = -t118 * t233 + t121 * t162 - t151 * t202 - t151 - t224;
t81 = t120 * t164 + t123 * t231 - t156 * t202 + t144 - t156;
t80 = t119 * t163 + t122 * t232 - t155 * t202 + t143 - t155;
t79 = t118 * t162 + t121 * t233 - t154 * t202 + t142 - t154;
t72 = (-t101 * t128 + t102 * t131) * t111;
t71 = (t100 * t130 - t127 * t99) * t110;
t70 = (-t126 * t97 + t129 * t98) * t109;
t66 = pkin(2) * t69;
t65 = pkin(2) * t68;
t64 = pkin(2) * t67;
t57 = (-t128 * t87 + t131 * t90) * t111;
t56 = (-t127 * t86 + t130 * t89) * t110;
t55 = (-t126 * t85 + t129 * t88) * t109;
t51 = (-t128 * t81 + t131 * t84) * t111;
t50 = (-t127 * t80 + t130 * t83) * t110;
t49 = (-t126 * t79 + t129 * t82) * t109;
t48 = (-t102 * t182 - t181 * t84 + t90) * t262;
t47 = (-t100 * t182 - t180 * t83 + t89) * t263;
t46 = (-t179 * t82 - t182 * t98 + t88) * t264;
t45 = (-t101 * t182 - t181 * t81 + t87) * t262;
t44 = (-t180 * t80 - t182 * t99 + t86) * t263;
t43 = (-t179 * t79 - t182 * t97 + t85) * t264;
t42 = (t102 * t108 + t160 * t84 - t259 * t90) * t111;
t41 = (t100 * t107 + t160 * t83 - t260 * t89) * t110;
t40 = (t106 * t98 + t160 * t82 - t261 * t88) * t109;
t39 = (t101 * t108 + t160 * t81 - t259 * t87) * t111;
t38 = (t107 * t99 + t160 * t80 - t260 * t86) * t110;
t37 = (t106 * t97 + t160 * t79 - t261 * t85) * t109;
t36 = (t102 * t105 + t108 * t84 - t258 * t90) * t111;
t35 = (t100 * t104 + t107 * t83 - t258 * t89) * t110;
t34 = (t103 * t98 + t106 * t82 - t258 * t88) * t109;
t33 = (t101 * t105 + t108 * t81 - t258 * t87) * t111;
t32 = (t104 * t99 + t107 * t80 - t258 * t86) * t110;
t31 = (t103 * t97 + t106 * t79 - t258 * t85) * t109;
t30 = t66 - t275;
t29 = t65 - t274;
t28 = t64 - t273;
t27 = (-t181 * t51 - t182 * t72 + t57) * m(3);
t26 = (-t180 * t50 - t182 * t71 + t56) * m(3);
t25 = (-t179 * t49 - t182 * t70 + t55) * m(3);
t24 = t108 * t72 + t160 * t51 - t259 * t57;
t23 = t107 * t71 + t160 * t50 - t260 * t56;
t22 = t106 * t70 + t160 * t49 - t261 * t55;
t21 = t105 * t72 + t108 * t51 - t258 * t57;
t20 = t104 * t71 + t107 * t50 - t258 * t56;
t19 = t103 * t70 + t106 * t49 - t258 * t55;
t18 = ((0.2e1 * (t66 + (-t252 / 0.2e1 - t249 / 0.2e1) * t111) * t243 - (pkin(2) * t30 - t60) * t231 - qJ(3,1) * t275) * t111 + (-qJ(3,1) + t221 - t243) * t111 * t275) * t69;
t17 = ((0.2e1 * (t65 + (-t253 / 0.2e1 - t250 / 0.2e1) * t110) * t239 - (pkin(2) * t29 - t59) * t232 - qJ(3,2) * t274) * t110 + (-qJ(3,2) + t222 - t239) * t110 * t274) * t68;
t16 = ((0.2e1 * (t64 + (-t254 / 0.2e1 - t251 / 0.2e1) * t109) * t235 - (pkin(2) * t28 - t58) * t233 - qJ(3,3) * t273) * t109 + (-qJ(3,3) + t223 - t235) * t109 * t273) * t67;
t15 = ((t30 - t275) * t231 + (-t164 * t69 - t215 + t69) * qJ(3,1)) * t255;
t14 = ((t29 - t274) * t232 + (-t163 * t68 - t216 + t68) * qJ(3,2)) * t256;
t13 = ((t28 - t273) * t233 + (-t162 * t67 - t217 + t67) * qJ(3,3)) * t257;
t9 = (-t169 * t69 ^ 2 - t134 * t178 + t18 * t182 + t181 * t246 - t15) * m(3);
t8 = (-t68 ^ 2 * t168 - t133 * t177 + t182 * t17 + t180 * t247 - t14) * m(3);
t7 = (-t167 * t67 ^ 2 - t132 * t176 + t16 * t182 + t179 * t248 - t13) * m(3);
t6 = -t105 * t18 - t108 * t12 + (-(rSges(2,1) * t137 - rSges(2,2) * t134) * t181 + (rSges(2,1) * t134 + rSges(2,2) * t137) * t178) * m(2) + (t182 * t15 + t69 * t169 * t269 + (-t134 * t169 - t137 * t182) * t181 + (t134 * t182 - t137 * t169) * t178) * m(3);
t5 = -t104 * t17 - t107 * t11 + (-(rSges(2,1) * t136 - rSges(2,2) * t133) * t180 + (rSges(2,1) * t133 + rSges(2,2) * t136) * t177) * m(2) + (t182 * t14 + t68 * t168 * t270 + (-t133 * t168 - t136 * t182) * t180 + (t133 * t182 - t136 * t168) * t177) * m(3);
t4 = -t106 * t10 - t103 * t16 + (-(rSges(2,1) * t135 - rSges(2,2) * t132) * t179 + (rSges(2,1) * t132 + rSges(2,2) * t135) * t176) * m(2) + (t182 * t13 + t67 * t167 * t271 + (-t132 * t167 - t135 * t182) * t179 + (t132 * t182 - t135 * t167) * t176) * m(3);
t3 = -t108 * t18 + t15 * t259 - t69 * ((t141 * t69 + t272 * t275) * t178 + t69 * t181 * t218) - t246 * t160;
t2 = -t107 * t17 + t14 * t260 - t68 * ((t141 * t68 + t272 * t274) * t177 + t68 * t180 * t219) - t247 * t160;
t1 = -t106 * t16 + t13 * t261 - t67 * ((t141 * t67 + t272 * t273) * t176 + t67 * t179 * t220) - t248 * t160;
t52 = [(-t124 * t173 - t125 * t165 - g(1) + t175) * m(4) + ((t101 * t33 + t39 * t81 + t45 * t87) * t96 + (t102 * t33 + t39 * t84 + t45 * t90) * t93 + t81 * t3 + t101 * t6 + t87 * t9) * t111 + ((t32 * t99 + t38 * t80 + t44 * t86) * t95 + (t100 * t32 + t38 * t83 + t44 * t89) * t92 + t80 * t2 + t99 * t5 + t86 * t8) * t110 + ((t31 * t97 + t37 * t79 + t43 * t85) * t94 + (t31 * t98 + t37 * t82 + t43 * t88) * t91 + t79 * t1 + t97 * t4 + t85 * t7) * t109; (-t124 * t165 + t125 * t173 - g(2) + t174) * m(4) + ((t101 * t36 + t42 * t81 + t48 * t87) * t96 + (t102 * t36 + t42 * t84 + t48 * t90) * t93 + t84 * t3 + t102 * t6 + t90 * t9) * t111 + ((t35 * t99 + t41 * t80 + t47 * t86) * t95 + (t100 * t35 + t41 * t83 + t47 * t89) * t92 + t83 * t2 + t100 * t5 + t89 * t8) * t110 + ((t34 * t97 + t40 * t79 + t46 * t85) * t94 + (t34 * t98 + t40 * t82 + t46 * t88) * t91 + t82 * t1 + t98 * t4 + t88 * t7) * t109; Icges(4,3) * t173 + t49 * t1 + t50 * t2 + t51 * t3 + t70 * t4 + t71 * t5 + t55 * t7 + t56 * t8 + t57 * t9 + t72 * t6 + (-t124 * t175 + t125 * t174 + (t188 ^ 2 + t189 ^ 2) * t173 + (g(1) * t189 + g(2) * t188) * t158 + (g(1) * t188 - g(2) * t189) * t159) * m(4) + ((t101 * t21 + t24 * t81 + t27 * t87) * t96 + (t102 * t21 + t24 * t84 + t27 * t90) * t93) * t111 + ((t20 * t99 + t23 * t80 + t26 * t86) * t95 + (t100 * t20 + t23 * t83 + t26 * t89) * t92) * t110 + ((t19 * t97 + t22 * t79 + t25 * t85) * t94 + (t19 * t98 + t22 * t82 + t25 * t88) * t91) * t109;];
tauX  = t52;
