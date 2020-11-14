% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:38
% EndTime: 2020-08-06 19:05:44
% DurationCPUTime: 6.19s
% Computational Cost: add. (30168->480), mult. (46056->788), div. (6867->6), fcn. (38043->18), ass. (0->306)
t205 = pkin(1) + pkin(2);
t194 = cos(qJ(2,3));
t202 = xDP(3);
t203 = xDP(2);
t204 = xDP(1);
t247 = t194 * t205;
t188 = sin(qJ(2,3));
t294 = qJ(3,3) * t188;
t142 = t247 + t294;
t189 = sin(qJ(1,3));
t195 = cos(qJ(1,3));
t224 = -pkin(4) * t189 + t142 * t195;
t133 = 0.1e1 / t142;
t207 = 0.1e1 / qJ(3,3);
t280 = t133 * t207;
t259 = t188 * t189;
t333 = pkin(4) * t195;
t130 = qJ(3,3) * t259 + t333;
t182 = legFrame(3,2);
t160 = sin(t182);
t163 = cos(t182);
t172 = t194 ^ 2;
t256 = t189 * t205;
t257 = t188 * t205;
t274 = t160 * qJ(3,3);
t82 = (t163 * t256 - t274) * t172 + (t130 * t163 + t160 * t257) * t194 + t274;
t268 = t163 * qJ(3,3);
t85 = (-t160 * t256 - t268) * t172 + (-t130 * t160 + t163 * t257) * t194 + t268;
t46 = (t194 * t202 * t224 + t203 * t85 + t204 * t82) * t280;
t313 = t205 * t46;
t196 = cos(qJ(2,2));
t245 = t196 * t205;
t190 = sin(qJ(2,2));
t296 = qJ(3,2) * t190;
t143 = t245 + t296;
t191 = sin(qJ(1,2));
t197 = cos(qJ(1,2));
t223 = -pkin(4) * t191 + t143 * t197;
t134 = 0.1e1 / t143;
t209 = 0.1e1 / qJ(3,2);
t279 = t134 * t209;
t255 = t190 * t191;
t332 = pkin(4) * t197;
t131 = qJ(3,2) * t255 + t332;
t183 = legFrame(2,2);
t161 = sin(t183);
t164 = cos(t183);
t173 = t196 ^ 2;
t252 = t191 * t205;
t253 = t190 * t205;
t272 = t161 * qJ(3,2);
t83 = (t164 * t252 - t272) * t173 + (t131 * t164 + t161 * t253) * t196 + t272;
t266 = t164 * qJ(3,2);
t86 = (-t161 * t252 - t266) * t173 + (-t131 * t161 + t164 * t253) * t196 + t266;
t47 = (t196 * t202 * t223 + t203 * t86 + t204 * t83) * t279;
t312 = t205 * t47;
t198 = cos(qJ(2,1));
t243 = t198 * t205;
t192 = sin(qJ(2,1));
t298 = qJ(3,1) * t192;
t144 = t243 + t298;
t193 = sin(qJ(1,1));
t199 = cos(qJ(1,1));
t222 = -pkin(4) * t193 + t144 * t199;
t135 = 0.1e1 / t144;
t211 = 0.1e1 / qJ(3,1);
t278 = t135 * t211;
t251 = t192 * t193;
t331 = pkin(4) * t199;
t132 = qJ(3,1) * t251 + t331;
t184 = legFrame(1,2);
t162 = sin(t184);
t165 = cos(t184);
t174 = t198 ^ 2;
t248 = t193 * t205;
t249 = t192 * t205;
t270 = t162 * qJ(3,1);
t84 = (t165 * t248 - t270) * t174 + (t132 * t165 + t162 * t249) * t198 + t270;
t264 = t165 * qJ(3,1);
t87 = (-t162 * t248 - t264) * t174 + (-t132 * t162 + t165 * t249) * t198 + t264;
t48 = (t198 * t202 * t222 + t203 * t87 + t204 * t84) * t278;
t311 = t205 * t48;
t289 = t224 * t207;
t106 = t142 * t189 + t333;
t293 = qJ(3,3) * t194;
t216 = -t257 + t293;
t98 = -t106 * t160 - t163 * t216;
t307 = t207 * t98;
t97 = t106 * t163 - t160 * t216;
t308 = t207 * t97;
t64 = t202 * t289 + t203 * t307 + t204 * t308;
t339 = 0.2e1 * t64;
t288 = t223 * t209;
t107 = t143 * t191 + t332;
t295 = qJ(3,2) * t196;
t217 = -t253 + t295;
t100 = -t107 * t161 - t164 * t217;
t292 = t100 * t209;
t99 = t107 * t164 - t161 * t217;
t304 = t209 * t99;
t65 = t202 * t288 + t203 * t292 + t204 * t304;
t338 = 0.2e1 * t65;
t287 = t222 * t211;
t108 = t144 * t193 + t331;
t297 = qJ(3,1) * t198;
t218 = -t249 + t297;
t102 = -t108 * t162 - t165 * t218;
t290 = t102 * t211;
t101 = t108 * t165 - t162 * t218;
t291 = t101 * t211;
t66 = t202 * t287 + t203 * t290 + t204 * t291;
t337 = 0.2e1 * t66;
t79 = (-t189 * t202 + (-t160 * t203 + t163 * t204) * t195) * t133;
t336 = pkin(4) * t79;
t80 = (-t191 * t202 + (-t161 * t203 + t164 * t204) * t197) * t134;
t335 = pkin(4) * t80;
t81 = (-t193 * t202 + (-t162 * t203 + t165 * t204) * t199) * t135;
t334 = pkin(4) * t81;
t330 = Ifges(2,1) + Ifges(3,1);
t329 = Ifges(2,6) - Ifges(3,6);
t328 = mrSges(3,2) * t188;
t327 = mrSges(3,2) * t190;
t326 = mrSges(3,2) * t192;
t325 = mrSges(3,3) * qJ(3,1);
t324 = mrSges(3,3) * qJ(3,2);
t323 = mrSges(3,3) * qJ(3,3);
t322 = t133 * t79;
t321 = t134 * t80;
t320 = t135 * t81;
t319 = t172 * t79;
t318 = t173 * t80;
t317 = t174 * t81;
t22 = -t64 + t313;
t316 = t188 * t22;
t23 = -t65 + t312;
t315 = t190 * t23;
t24 = -t66 + t311;
t314 = t192 * t24;
t310 = t207 * t82;
t309 = t207 * t85;
t306 = t209 * t83;
t305 = t209 * t86;
t303 = t211 * t84;
t302 = t211 * t87;
t157 = m(3) * qJ(3,3) + mrSges(3,3);
t301 = t64 * t157;
t158 = m(3) * qJ(3,2) + mrSges(3,3);
t300 = t65 * t158;
t159 = m(3) * qJ(3,1) + mrSges(3,3);
t299 = t66 * t159;
t148 = mrSges(3,2) * qJ(3,3) + t329;
t154 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t115 = t148 * t194 - t188 * t154;
t286 = t115 * t207;
t149 = mrSges(3,2) * qJ(3,2) + t329;
t116 = t149 * t196 - t190 * t154;
t285 = t116 * t209;
t150 = qJ(3,1) * mrSges(3,2) + t329;
t117 = t150 * t198 - t192 * t154;
t284 = t117 * t211;
t168 = 0.2e1 * t323;
t206 = qJ(3,3) ^ 2;
t214 = pkin(1) ^ 2;
t241 = 0.2e1 * mrSges(3,1) * pkin(1);
t225 = Ifges(3,2) + Ifges(2,3) + t241;
t121 = (t206 + t214) * m(3) + t168 + t225;
t283 = t121 * t207;
t169 = 0.2e1 * t324;
t208 = qJ(3,2) ^ 2;
t122 = (t208 + t214) * m(3) + t169 + t225;
t282 = t122 * t209;
t170 = 0.2e1 * t325;
t210 = qJ(3,1) ^ 2;
t123 = (t210 + t214) * m(3) + t170 + t225;
t281 = t123 * t211;
t166 = m(3) * pkin(1) + mrSges(3,1);
t156 = mrSges(2,1) + t166;
t277 = t156 * t194;
t276 = t156 * t196;
t275 = t156 * t198;
t273 = t160 * t195;
t271 = t161 * t197;
t269 = t162 * t199;
t267 = t163 * t195;
t265 = t164 * t197;
t263 = t165 * t199;
t262 = t166 * t207;
t261 = t166 * t209;
t260 = t166 * t211;
t258 = t188 * t194;
t254 = t190 * t196;
t250 = t192 * t198;
t246 = t194 * t207;
t244 = t196 * t209;
t242 = t198 * t211;
t240 = 0.2e1 * t258;
t239 = 0.2e1 * t254;
t238 = 0.2e1 * t250;
t237 = Ifges(1,3) + t330;
t236 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t235 = t195 * t328;
t234 = t207 * t328;
t233 = t197 * t327;
t232 = t209 * t327;
t231 = t199 * t326;
t230 = t211 * t326;
t229 = t224 * t246;
t228 = t223 * t244;
t227 = t222 * t242;
t226 = -t214 + (-0.2e1 * pkin(1) - pkin(2)) * pkin(2);
t127 = g(1) * t163 - g(2) * t160;
t221 = g(3) * t195 + t127 * t189;
t128 = g(1) * t164 - g(2) * t161;
t220 = g(3) * t197 + t128 * t191;
t129 = g(1) * t165 - g(2) * t162;
t219 = g(3) * t199 + t129 * t193;
t215 = m(3) * t214 + Ifges(2,2) + Ifges(3,3) + t241 - t330;
t212 = pkin(4) ^ 2;
t200 = g(3) * mrSges(1,1);
t187 = xDDP(1);
t186 = xDDP(2);
t185 = xDDP(3);
t180 = m(3) * t210;
t178 = m(3) * t208;
t176 = m(3) * t206;
t167 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t155 = t167 * g(3);
t153 = -mrSges(2,2) + t159;
t152 = -mrSges(2,2) + t158;
t151 = -mrSges(2,2) + t157;
t138 = qJ(3,1) * t166 + t236;
t137 = qJ(3,2) * t166 + t236;
t136 = qJ(3,3) * t166 + t236;
t126 = g(1) * t162 + g(2) * t165;
t125 = g(1) * t161 + g(2) * t164;
t124 = g(1) * t160 + g(2) * t163;
t120 = -t180 + t215 - 0.2e1 * t325;
t119 = -t178 + t215 - 0.2e1 * t324;
t118 = -t176 + t215 - 0.2e1 * t323;
t96 = t120 * t174 + t138 * t238 + t170 + t180 + t237;
t95 = t119 * t173 + t137 * t239 + t169 + t178 + t237;
t94 = t118 * t172 + t136 * t240 + t168 + t176 + t237;
t78 = t81 ^ 2;
t77 = t80 ^ 2;
t76 = t79 ^ 2;
t75 = t192 * t334;
t74 = t190 * t335;
t73 = t188 * t336;
t72 = m(3) * t287 + (-mrSges(3,2) * t251 - t166 * t227) * t135;
t71 = m(3) * t288 + (-mrSges(3,2) * t255 - t166 * t228) * t134;
t70 = m(3) * t289 + (-mrSges(3,2) * t259 - t166 * t229) * t133;
t69 = -t222 * t260 + (-t117 * t193 + t123 * t227) * t135;
t68 = -t223 * t261 + (-t116 * t191 + t122 * t228) * t134;
t67 = -t224 * t262 + (-t115 * t189 + t121 * t229) * t133;
t63 = t222 * t230 + (t117 * t227 - t193 * t96) * t135;
t62 = t223 * t232 + (t116 * t228 - t191 * t95) * t134;
t61 = t224 * t234 + (t115 * t229 - t189 * t94) * t133;
t60 = m(3) * t291 + (t165 * t231 - t84 * t260) * t135;
t59 = m(3) * t304 + (t164 * t233 - t83 * t261) * t134;
t58 = m(3) * t308 + (t163 * t235 - t82 * t262) * t133;
t57 = m(3) * t290 + (-t162 * t231 - t87 * t260) * t135;
t56 = m(3) * t292 + (-t161 * t233 - t86 * t261) * t134;
t55 = m(3) * t307 + (-t160 * t235 - t85 * t262) * t133;
t54 = -t101 * t260 + (t117 * t263 + t84 * t281) * t135;
t53 = -t99 * t261 + (t116 * t265 + t83 * t282) * t134;
t52 = -t97 * t262 + (t115 * t267 + t82 * t283) * t133;
t51 = -t102 * t260 + (-t117 * t269 + t87 * t281) * t135;
t50 = -t100 * t261 + (-t116 * t271 + t86 * t282) * t134;
t49 = -t98 * t262 + (-t115 * t273 + t85 * t283) * t133;
t45 = t48 ^ 2;
t44 = t47 ^ 2;
t43 = t46 ^ 2;
t36 = t138 * t48;
t35 = t137 * t47;
t34 = t136 * t46;
t33 = t101 * t230 + (t96 * t263 + t84 * t284) * t135;
t32 = t99 * t232 + (t95 * t265 + t83 * t285) * t134;
t31 = t97 * t234 + (t94 * t267 + t82 * t286) * t133;
t30 = t102 * t230 + (-t96 * t269 + t87 * t284) * t135;
t29 = t100 * t232 + (-t95 * t271 + t86 * t285) * t134;
t28 = t98 * t234 + (-t94 * t273 + t85 * t286) * t133;
t27 = t75 + t311;
t26 = t74 + t312;
t25 = t73 + t313;
t21 = (-t48 * t297 + t314) * pkin(4) + ((qJ(3,1) + t205) * (-qJ(3,1) + t205) * t174 + t205 * qJ(3,1) * t238 + t210 + t212) * t81;
t20 = (-t47 * t295 + t315) * pkin(4) + ((qJ(3,2) + t205) * (-qJ(3,2) + t205) * t173 + t205 * qJ(3,2) * t239 + t208 + t212) * t80;
t19 = (-t46 * t293 + t316) * pkin(4) + ((qJ(3,3) + t205) * (-qJ(3,3) + t205) * t172 + t205 * qJ(3,3) * t240 + t206 + t212) * t79;
t18 = (-t334 + (t218 + t297) * t48 + (t337 - t311) * t192) * t320;
t17 = (-t335 + (t217 + t295) * t47 + (t338 - t312) * t190) * t321;
t16 = (-t336 + (t216 + t293) * t46 + (t339 - t313) * t188) * t322;
t15 = (((-t210 + t226) * t48 + t205 * t66) * t48 + t27 * t66 + (t218 * t48 * pkin(4) - t21) * t81) * t211;
t14 = (((-t208 + t226) * t47 + t205 * t65) * t47 + t26 * t65 + (t217 * t47 * pkin(4) - t20) * t80) * t209;
t13 = (((-t206 + t226) * t46 + t205 * t64) * t46 + t25 * t64 + (t216 * t46 * pkin(4) - t19) * t79) * t207;
t12 = -t21 * t242 * t320 + ((-(t75 + t24) * t243 + (pkin(4) * t317 - t314) * qJ(3,1)) * t48 + (t198 * t27 + t48 * t298) * t66) * t278;
t11 = -t20 * t244 * t321 + ((-(t74 + t23) * t245 + (pkin(4) * t318 - t315) * qJ(3,2)) * t47 + (t196 * t26 + t47 * t296) * t65) * t279;
t10 = -t19 * t246 * t322 + ((-(t73 + t22) * t247 + (pkin(4) * t319 - t316) * qJ(3,3)) * t46 + (t194 * t25 + t46 * t294) * t64) * t280;
t9 = t166 * t12 + (t174 * t78 - t45 - t78) * t159 + (t198 * t126 - t15) * m(3) + (-t78 * t166 * t198 - m(3) * t219 - mrSges(3,2) * t18) * t192;
t8 = t166 * t11 + (t173 * t77 - t44 - t77) * t158 + (t196 * t125 - t14) * m(3) + (-t77 * t166 * t196 - m(3) * t220 - mrSges(3,2) * t17) * t190;
t7 = t166 * t10 + (t172 * t76 - t43 - t76) * t157 + (t194 * t124 - t13) * m(3) + (-t76 * t166 * t194 - m(3) * t221 - mrSges(3,2) * t16) * t188;
t6 = -t117 * t18 - t123 * t12 + t166 * t15 + 0.2e1 * t48 * t299 + (-t126 * t156 - t219 * t153) * t198 + t192 * (-t126 * t153 + t219 * t156) + (t120 * t250 - 0.2e1 * t138 * t174 + t138) * t78;
t5 = -t116 * t17 - t122 * t11 + t166 * t14 + 0.2e1 * t47 * t300 + (-t125 * t156 - t220 * t152) * t196 + t190 * (-t125 * t152 + t220 * t156) + (t119 * t254 - 0.2e1 * t137 * t173 + t137) * t77;
t4 = -t115 * t16 - t121 * t10 + t166 * t13 + 0.2e1 * t46 * t301 + (-t124 * t156 - t221 * t151) * t194 + t188 * (-t124 * t151 + t221 * t156) + (t118 * t258 - 0.2e1 * t136 * t172 + t136) * t76;
t3 = -t96 * t18 - t117 * t12 + 0.4e1 * (t36 - t299 / 0.2e1) * t317 + (mrSges(3,2) * t337 - t154 * t48) * t48 * t198 - 0.2e1 * (t36 - t299) * t81 + t155 * t199 + t193 * (g(3) * t275 + t200) + (-mrSges(3,2) * t15 + 0.2e1 * (-t120 * t48 + t166 * t66) * t81 * t198 - t45 * t150 + t193 * g(3) * t153) * t192 + ((-t153 * t192 - mrSges(1,1) - t275) * t199 + t193 * t167) * t129;
t2 = -t95 * t17 - t116 * t11 + 0.4e1 * (t35 - t300 / 0.2e1) * t318 + (mrSges(3,2) * t338 - t154 * t47) * t47 * t196 - 0.2e1 * (t35 - t300) * t80 + t155 * t197 + t191 * (g(3) * t276 + t200) + (-mrSges(3,2) * t14 + 0.2e1 * (-t119 * t47 + t166 * t65) * t80 * t196 - t44 * t149 + t191 * g(3) * t152) * t190 + ((-t152 * t190 - mrSges(1,1) - t276) * t197 + t191 * t167) * t128;
t1 = -t94 * t16 - t115 * t10 + 0.4e1 * (t34 - t301 / 0.2e1) * t319 + (mrSges(3,2) * t339 - t154 * t46) * t46 * t194 - 0.2e1 * (t34 - t301) * t79 + t155 * t195 + t189 * (g(3) * t277 + t200) + (-mrSges(3,2) * t13 + 0.2e1 * (-t118 * t46 + t166 * t64) * t79 * t194 - t43 * t148 + t189 * g(3) * t151) * t188 + ((-t151 * t188 - mrSges(1,1) - t277) * t195 + t189 * t167) * t127;
t37 = [t9 * t291 + t7 * t308 + t8 * t304 - g(1) * m(4) + (t60 * t290 + t59 * t292 + t58 * t307) * t186 + (t60 * t287 + t59 * t288 + t58 * t289) * t185 + (t60 * t291 + t59 * t304 + t58 * t308 + m(4)) * t187 + ((t33 * t263 + t54 * t303) * t187 + (-t33 * t269 + t54 * t302) * t186 + (-t193 * t33 + t54 * t227) * t185 + t3 * t263 + t6 * t303) * t135 + ((t32 * t265 + t53 * t306) * t187 + (-t32 * t271 + t53 * t305) * t186 + (-t191 * t32 + t228 * t53) * t185 + t2 * t265 + t5 * t306) * t134 + ((t31 * t267 + t52 * t310) * t187 + (-t31 * t273 + t52 * t309) * t186 + (-t189 * t31 + t52 * t229) * t185 + t1 * t267 + t4 * t310) * t133; t8 * t292 + t9 * t290 + t7 * t307 - g(2) * m(4) + (t57 * t291 + t56 * t304 + t55 * t308) * t187 + (t57 * t287 + t56 * t288 + t55 * t289) * t185 + (t57 * t290 + t56 * t292 + t55 * t307 + m(4)) * t186 + ((t30 * t263 + t51 * t303) * t187 + (-t30 * t269 + t51 * t302) * t186 + (-t193 * t30 + t51 * t227) * t185 - t3 * t269 + t6 * t302) * t135 + ((t29 * t265 + t50 * t306) * t187 + (-t29 * t271 + t50 * t305) * t186 + (-t191 * t29 + t228 * t50) * t185 - t2 * t271 + t5 * t305) * t134 + ((t28 * t267 + t49 * t310) * t187 + (-t28 * t273 + t49 * t309) * t186 + (-t189 * t28 + t49 * t229) * t185 - t1 * t273 + t4 * t309) * t133; t7 * t289 + t8 * t288 + t9 * t287 - g(3) * m(4) + (t72 * t291 + t71 * t304 + t70 * t308) * t187 + (t72 * t290 + t71 * t292 + t70 * t307) * t186 + (t72 * t287 + t71 * t288 + t70 * t289 + m(4)) * t185 + ((t63 * t263 + t69 * t303) * t187 + (-t63 * t269 + t69 * t302) * t186 + (-t193 * t63 + t69 * t227) * t185 - t193 * t3 + t6 * t227) * t135 + ((t62 * t265 + t68 * t306) * t187 + (-t62 * t271 + t68 * t305) * t186 + (-t191 * t62 + t228 * t68) * t185 - t191 * t2 + t5 * t228) * t134 + ((t61 * t267 + t67 * t310) * t187 + (-t61 * t273 + t67 * t309) * t186 + (-t189 * t61 + t67 * t229) * t185 - t189 * t1 + t4 * t229) * t133;];
tauX  = t37;
