% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V1G3A0
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
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:09:56
% EndTime: 2020-08-06 19:10:03
% DurationCPUTime: 6.98s
% Computational Cost: add. (30168->477), mult. (45402->784), div. (6867->6), fcn. (37389->18), ass. (0->305)
t204 = pkin(1) + pkin(2);
t193 = cos(qJ(2,3));
t249 = t193 * t204;
t187 = sin(qJ(2,3));
t290 = qJ(3,3) * t187;
t139 = t249 + t290;
t188 = sin(qJ(1,3));
t194 = cos(qJ(1,3));
t106 = pkin(4) * t194 + t139 * t188;
t201 = xDP(3);
t202 = xDP(2);
t203 = xDP(1);
t130 = 0.1e1 / t139;
t206 = 0.1e1 / qJ(3,3);
t276 = t130 * t206;
t257 = t187 * t194;
t332 = pkin(4) * t188;
t127 = qJ(3,3) * t257 - t332;
t181 = legFrame(3,2);
t159 = sin(t181);
t162 = cos(t181);
t171 = t193 ^ 2;
t247 = t194 * t204;
t256 = t187 * t204;
t273 = t159 * qJ(3,3);
t82 = (t162 * t247 - t273) * t171 + (t127 * t162 + t159 * t256) * t193 + t273;
t267 = t162 * qJ(3,3);
t85 = (-t159 * t247 - t267) * t171 + (-t127 * t159 + t162 * t256) * t193 + t267;
t46 = (-t106 * t193 * t201 + t202 * t85 + t203 * t82) * t276;
t309 = t204 * t46;
t195 = cos(qJ(2,2));
t246 = t195 * t204;
t189 = sin(qJ(2,2));
t292 = qJ(3,2) * t189;
t140 = t246 + t292;
t190 = sin(qJ(1,2));
t196 = cos(qJ(1,2));
t107 = pkin(4) * t196 + t140 * t190;
t131 = 0.1e1 / t140;
t208 = 0.1e1 / qJ(3,2);
t275 = t131 * t208;
t254 = t189 * t196;
t331 = pkin(4) * t190;
t128 = qJ(3,2) * t254 - t331;
t182 = legFrame(2,2);
t160 = sin(t182);
t163 = cos(t182);
t172 = t195 ^ 2;
t244 = t196 * t204;
t253 = t189 * t204;
t271 = t160 * qJ(3,2);
t83 = (t163 * t244 - t271) * t172 + (t128 * t163 + t160 * t253) * t195 + t271;
t265 = t163 * qJ(3,2);
t86 = (-t160 * t244 - t265) * t172 + (-t128 * t160 + t163 * t253) * t195 + t265;
t47 = (-t107 * t195 * t201 + t202 * t86 + t203 * t83) * t275;
t308 = t204 * t47;
t197 = cos(qJ(2,1));
t243 = t197 * t204;
t191 = sin(qJ(2,1));
t294 = qJ(3,1) * t191;
t141 = t243 + t294;
t192 = sin(qJ(1,1));
t198 = cos(qJ(1,1));
t108 = pkin(4) * t198 + t141 * t192;
t132 = 0.1e1 / t141;
t210 = 0.1e1 / qJ(3,1);
t274 = t132 * t210;
t251 = t191 * t198;
t330 = pkin(4) * t192;
t129 = qJ(3,1) * t251 - t330;
t183 = legFrame(1,2);
t161 = sin(t183);
t164 = cos(t183);
t173 = t197 ^ 2;
t241 = t198 * t204;
t250 = t191 * t204;
t269 = t161 * qJ(3,1);
t84 = (t164 * t241 - t269) * t173 + (t129 * t164 + t161 * t250) * t197 + t269;
t263 = t164 * qJ(3,1);
t87 = (-t161 * t241 - t263) * t173 + (-t129 * t161 + t164 * t250) * t197 + t263;
t48 = (-t108 * t197 * t201 + t202 * t87 + t203 * t84) * t274;
t307 = t204 * t48;
t341 = t139 * t194 - t332;
t340 = t140 * t196 - t331;
t339 = t141 * t198 - t330;
t285 = t106 * t206;
t289 = qJ(3,3) * t193;
t218 = -t256 + t289;
t100 = -t341 * t159 - t218 * t162;
t288 = t100 * t206;
t97 = -t218 * t159 + t341 * t162;
t304 = t206 * t97;
t64 = -t201 * t285 + t202 * t288 + t203 * t304;
t338 = 0.2e1 * t64;
t284 = t107 * t208;
t291 = qJ(3,2) * t195;
t219 = -t253 + t291;
t101 = -t340 * t160 - t219 * t163;
t287 = t101 * t208;
t98 = -t219 * t160 + t340 * t163;
t301 = t208 * t98;
t65 = -t201 * t284 + t202 * t287 + t203 * t301;
t337 = 0.2e1 * t65;
t283 = t108 * t210;
t293 = qJ(3,1) * t197;
t220 = -t250 + t293;
t102 = -t339 * t161 - t220 * t164;
t286 = t102 * t210;
t99 = -t220 * t161 + t339 * t164;
t298 = t210 * t99;
t66 = -t201 * t283 + t202 * t286 + t203 * t298;
t336 = 0.2e1 * t66;
t79 = (-t194 * t201 + (t159 * t202 - t162 * t203) * t188) * t130;
t335 = pkin(4) * t79;
t80 = (-t196 * t201 + (t160 * t202 - t163 * t203) * t190) * t131;
t334 = pkin(4) * t80;
t81 = (-t198 * t201 + (t161 * t202 - t164 * t203) * t192) * t132;
t333 = pkin(4) * t81;
t326 = Ifges(2,1) + Ifges(3,1);
t325 = Ifges(2,6) - Ifges(3,6);
t324 = mrSges(3,2) * t187;
t323 = mrSges(3,2) * t189;
t322 = mrSges(3,2) * t191;
t321 = mrSges(3,3) * qJ(3,1);
t320 = mrSges(3,3) * qJ(3,2);
t319 = mrSges(3,3) * qJ(3,3);
t318 = t130 * t79;
t317 = t131 * t80;
t316 = t132 * t81;
t315 = t171 * t79;
t314 = t172 * t80;
t313 = t173 * t81;
t22 = -t64 + t309;
t312 = t187 * t22;
t23 = -t65 + t308;
t311 = t189 * t23;
t24 = -t66 + t307;
t310 = t191 * t24;
t306 = t206 * t82;
t305 = t206 * t85;
t303 = t208 * t83;
t302 = t208 * t86;
t300 = t210 * t84;
t299 = t210 * t87;
t156 = m(3) * qJ(3,3) + mrSges(3,3);
t297 = t64 * t156;
t157 = m(3) * qJ(3,2) + mrSges(3,3);
t296 = t65 * t157;
t158 = m(3) * qJ(3,1) + mrSges(3,3);
t295 = t66 * t158;
t148 = mrSges(3,2) * qJ(3,3) + t325;
t154 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t112 = t148 * t193 - t187 * t154;
t282 = t112 * t206;
t149 = mrSges(3,2) * qJ(3,2) + t325;
t113 = t149 * t195 - t189 * t154;
t281 = t113 * t208;
t150 = qJ(3,1) * mrSges(3,2) + t325;
t114 = t150 * t197 - t191 * t154;
t280 = t114 * t210;
t167 = 0.2e1 * t319;
t205 = qJ(3,3) ^ 2;
t213 = pkin(1) ^ 2;
t240 = 0.2e1 * mrSges(3,1) * pkin(1);
t224 = Ifges(3,2) + Ifges(2,3) + t240;
t118 = (t205 + t213) * m(3) + t167 + t224;
t279 = t118 * t206;
t168 = 0.2e1 * t320;
t207 = qJ(3,2) ^ 2;
t119 = (t207 + t213) * m(3) + t168 + t224;
t278 = t119 * t208;
t169 = 0.2e1 * t321;
t209 = qJ(3,1) ^ 2;
t120 = (t209 + t213) * m(3) + t169 + t224;
t277 = t120 * t210;
t272 = t159 * t188;
t270 = t160 * t190;
t268 = t161 * t192;
t266 = t162 * t188;
t264 = t163 * t190;
t262 = t164 * t192;
t165 = m(3) * pkin(1) + mrSges(3,1);
t261 = t165 * t206;
t260 = t165 * t208;
t259 = t165 * t210;
t258 = t187 * t193;
t255 = t189 * t195;
t252 = t191 * t197;
t248 = t193 * t206;
t245 = t195 * t208;
t242 = t197 * t210;
t239 = 0.2e1 * t258;
t238 = 0.2e1 * t255;
t237 = 0.2e1 * t252;
t236 = Ifges(1,3) + t326;
t235 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t234 = t188 * t324;
t233 = t206 * t324;
t232 = t190 * t323;
t231 = t208 * t323;
t230 = t192 * t322;
t229 = t210 * t322;
t228 = t106 * t248;
t227 = t107 * t245;
t226 = t108 * t242;
t225 = -t213 + (-0.2e1 * pkin(1) - pkin(2)) * pkin(2);
t124 = g(1) * t162 - g(2) * t159;
t223 = g(3) * t188 - t124 * t194;
t125 = g(1) * t163 - g(2) * t160;
t222 = g(3) * t190 - t125 * t196;
t126 = g(1) * t164 - g(2) * t161;
t221 = g(3) * t192 - t126 * t198;
t151 = -mrSges(2,2) + t156;
t155 = mrSges(2,1) + t165;
t217 = t151 * t187 + t155 * t193;
t152 = -mrSges(2,2) + t157;
t216 = t152 * t189 + t155 * t195;
t153 = -mrSges(2,2) + t158;
t215 = t153 * t191 + t155 * t197;
t214 = m(3) * t213 + Ifges(2,2) + Ifges(3,3) + t240 - t326;
t211 = pkin(4) ^ 2;
t199 = g(3) * mrSges(1,1);
t186 = xDDP(1);
t185 = xDDP(2);
t184 = xDDP(3);
t179 = m(3) * t209;
t177 = m(3) * t207;
t175 = m(3) * t205;
t166 = mrSges(1,2) - mrSges(3,2) - mrSges(2,3);
t135 = qJ(3,1) * t165 + t235;
t134 = qJ(3,2) * t165 + t235;
t133 = qJ(3,3) * t165 + t235;
t123 = g(1) * t161 + g(2) * t164;
t122 = g(1) * t160 + g(2) * t163;
t121 = g(1) * t159 + g(2) * t162;
t117 = -t179 + t214 - 0.2e1 * t321;
t116 = -t177 + t214 - 0.2e1 * t320;
t115 = -t175 + t214 - 0.2e1 * t319;
t90 = t117 * t173 + t135 * t237 + t169 + t179 + t236;
t89 = t116 * t172 + t134 * t238 + t168 + t177 + t236;
t88 = t115 * t171 + t133 * t239 + t167 + t175 + t236;
t78 = t81 ^ 2;
t77 = t80 ^ 2;
t76 = t79 ^ 2;
t75 = t191 * t333;
t74 = t189 * t334;
t73 = t187 * t335;
t72 = -m(3) * t283 + (-mrSges(3,2) * t251 + t165 * t226) * t132;
t71 = -m(3) * t284 + (-mrSges(3,2) * t254 + t165 * t227) * t131;
t70 = -m(3) * t285 + (-mrSges(3,2) * t257 + t165 * t228) * t130;
t69 = t108 * t259 + (-t114 * t198 - t120 * t226) * t132;
t68 = t107 * t260 + (-t113 * t196 - t119 * t227) * t131;
t67 = t106 * t261 + (-t112 * t194 - t118 * t228) * t130;
t63 = -t108 * t229 + (-t114 * t226 - t198 * t90) * t132;
t62 = -t107 * t231 + (-t113 * t227 - t196 * t89) * t131;
t61 = -t106 * t233 + (-t112 * t228 - t194 * t88) * t130;
t60 = m(3) * t298 + (-t164 * t230 - t84 * t259) * t132;
t59 = m(3) * t301 + (-t163 * t232 - t83 * t260) * t131;
t58 = m(3) * t304 + (-t162 * t234 - t82 * t261) * t130;
t57 = m(3) * t286 + (t161 * t230 - t87 * t259) * t132;
t56 = m(3) * t287 + (t160 * t232 - t86 * t260) * t131;
t55 = m(3) * t288 + (t159 * t234 - t85 * t261) * t130;
t54 = -t99 * t259 + (-t114 * t262 + t84 * t277) * t132;
t53 = -t98 * t260 + (-t113 * t264 + t83 * t278) * t131;
t52 = -t97 * t261 + (-t112 * t266 + t82 * t279) * t130;
t51 = -t102 * t259 + (t114 * t268 + t87 * t277) * t132;
t50 = -t101 * t260 + (t113 * t270 + t86 * t278) * t131;
t49 = -t100 * t261 + (t112 * t272 + t85 * t279) * t130;
t45 = t48 ^ 2;
t44 = t47 ^ 2;
t43 = t46 ^ 2;
t36 = t135 * t48;
t35 = t134 * t47;
t34 = t133 * t46;
t33 = t99 * t229 + (-t90 * t262 + t84 * t280) * t132;
t32 = t98 * t231 + (-t89 * t264 + t83 * t281) * t131;
t31 = t97 * t233 + (-t88 * t266 + t82 * t282) * t130;
t30 = t102 * t229 + (t90 * t268 + t87 * t280) * t132;
t29 = t101 * t231 + (t89 * t270 + t86 * t281) * t131;
t28 = t100 * t233 + (t88 * t272 + t85 * t282) * t130;
t27 = t75 + t307;
t26 = t74 + t308;
t25 = t73 + t309;
t21 = (-t48 * t293 + t310) * pkin(4) + ((qJ(3,1) + t204) * (-qJ(3,1) + t204) * t173 + t204 * qJ(3,1) * t237 + t209 + t211) * t81;
t20 = (-t47 * t291 + t311) * pkin(4) + ((qJ(3,2) + t204) * (-qJ(3,2) + t204) * t172 + t204 * qJ(3,2) * t238 + t207 + t211) * t80;
t19 = (-t46 * t289 + t312) * pkin(4) + ((qJ(3,3) + t204) * (-qJ(3,3) + t204) * t171 + t204 * qJ(3,3) * t239 + t205 + t211) * t79;
t18 = (-t333 + (t220 + t293) * t48 + (t336 - t307) * t191) * t316;
t17 = (-t334 + (t219 + t291) * t47 + (t337 - t308) * t189) * t317;
t16 = (-t335 + (t218 + t289) * t46 + (t338 - t309) * t187) * t318;
t15 = (((-t209 + t225) * t48 + t204 * t66) * t48 + t27 * t66 + (t220 * t48 * pkin(4) - t21) * t81) * t210;
t14 = (((-t207 + t225) * t47 + t204 * t65) * t47 + t26 * t65 + (t219 * t47 * pkin(4) - t20) * t80) * t208;
t13 = (((-t205 + t225) * t46 + t204 * t64) * t46 + t25 * t64 + (t218 * t46 * pkin(4) - t19) * t79) * t206;
t12 = -t21 * t242 * t316 + ((-(t75 + t24) * t243 + (pkin(4) * t313 - t310) * qJ(3,1)) * t48 + (t197 * t27 + t48 * t294) * t66) * t274;
t11 = -t20 * t245 * t317 + ((-(t74 + t23) * t246 + (pkin(4) * t314 - t311) * qJ(3,2)) * t47 + (t195 * t26 + t47 * t292) * t65) * t275;
t10 = -t19 * t248 * t318 + ((-(t73 + t22) * t249 + (pkin(4) * t315 - t312) * qJ(3,3)) * t46 + (t193 * t25 + t46 * t290) * t64) * t276;
t9 = t165 * t12 + (t173 * t78 - t45 - t78) * t158 + (t197 * t123 - t15) * m(3) + (-t78 * t165 * t197 + m(3) * t221 - mrSges(3,2) * t18) * t191;
t8 = t165 * t11 + (t172 * t77 - t44 - t77) * t157 + (t195 * t122 - t14) * m(3) + (-t77 * t165 * t195 + m(3) * t222 - mrSges(3,2) * t17) * t189;
t7 = t165 * t10 + (t171 * t76 - t43 - t76) * t156 + (t193 * t121 - t13) * m(3) + (-t76 * t165 * t193 + m(3) * t223 - mrSges(3,2) * t16) * t187;
t6 = -t114 * t18 - t120 * t12 + t165 * t15 + 0.2e1 * t48 * t295 + (-t123 * t155 + t221 * t153) * t197 + t191 * (-t123 * t153 - t221 * t155) + (t117 * t252 - 0.2e1 * t135 * t173 + t135) * t78;
t5 = -t113 * t17 - t119 * t11 + t165 * t14 + 0.2e1 * t47 * t296 + (-t122 * t155 + t222 * t152) * t195 + t189 * (-t122 * t152 - t222 * t155) + (t116 * t255 - 0.2e1 * t134 * t172 + t134) * t77;
t4 = -t112 * t16 - t118 * t10 + t165 * t13 + 0.2e1 * t46 * t297 + (-t121 * t155 + t223 * t151) * t193 + t187 * (-t121 * t151 - t223 * t155) + (t115 * t258 - 0.2e1 * t133 * t171 + t133) * t76;
t3 = -t90 * t18 - t114 * t12 + 0.4e1 * (t36 - t295 / 0.2e1) * t313 + t48 * (mrSges(3,2) * t336 - t154 * t48) * t197 - 0.2e1 * (t36 - t295) * t81 + t199 * t198 + (-mrSges(3,2) * t15 + 0.2e1 * (-t117 * t48 + t165 * t66) * t81 * t197 - t45 * t150) * t191 + (-t166 * t192 + t215 * t198) * g(3) + (t166 * t198 + (mrSges(1,1) + t215) * t192) * t126;
t2 = -t89 * t17 - t113 * t11 + 0.4e1 * (t35 - t296 / 0.2e1) * t314 + t47 * (mrSges(3,2) * t337 - t154 * t47) * t195 - 0.2e1 * (t35 - t296) * t80 + t199 * t196 + (-mrSges(3,2) * t14 + 0.2e1 * (-t116 * t47 + t165 * t65) * t80 * t195 - t44 * t149) * t189 + (-t166 * t190 + t216 * t196) * g(3) + (t166 * t196 + (mrSges(1,1) + t216) * t190) * t125;
t1 = -t88 * t16 - t112 * t10 + 0.4e1 * (t34 - t297 / 0.2e1) * t315 + t46 * (mrSges(3,2) * t338 - t154 * t46) * t193 - 0.2e1 * (t34 - t297) * t79 + t199 * t194 + (-mrSges(3,2) * t13 + 0.2e1 * (-t115 * t46 + t165 * t64) * t79 * t193 - t43 * t148) * t187 + (-t166 * t188 + t217 * t194) * g(3) + (t166 * t194 + (mrSges(1,1) + t217) * t188) * t124;
t37 = [t7 * t304 + t8 * t301 + t9 * t298 - g(1) * m(4) + (t60 * t286 + t59 * t287 + t58 * t288) * t185 + (-t60 * t283 - t59 * t284 - t58 * t285) * t184 + (t60 * t298 + t59 * t301 + t58 * t304 + m(4)) * t186 + ((-t33 * t262 + t54 * t300) * t186 + (t33 * t268 + t54 * t299) * t185 + (-t198 * t33 - t54 * t226) * t184 - t3 * t262 + t6 * t300) * t132 + ((-t32 * t264 + t53 * t303) * t186 + (t32 * t270 + t53 * t302) * t185 + (-t196 * t32 - t53 * t227) * t184 - t2 * t264 + t5 * t303) * t131 + ((-t31 * t266 + t52 * t306) * t186 + (t31 * t272 + t52 * t305) * t185 + (-t194 * t31 - t52 * t228) * t184 - t1 * t266 + t4 * t306) * t130; t7 * t288 + t8 * t287 + t9 * t286 - g(2) * m(4) + (t57 * t298 + t56 * t301 + t55 * t304) * t186 + (-t57 * t283 - t56 * t284 - t55 * t285) * t184 + (t57 * t286 + t56 * t287 + t55 * t288 + m(4)) * t185 + ((-t30 * t262 + t51 * t300) * t186 + (t30 * t268 + t51 * t299) * t185 + (-t198 * t30 - t51 * t226) * t184 + t3 * t268 + t6 * t299) * t132 + ((-t29 * t264 + t50 * t303) * t186 + (t29 * t270 + t50 * t302) * t185 + (-t196 * t29 - t50 * t227) * t184 + t2 * t270 + t5 * t302) * t131 + ((-t28 * t266 + t49 * t306) * t186 + (t28 * t272 + t49 * t305) * t185 + (-t194 * t28 - t49 * t228) * t184 + t1 * t272 + t4 * t305) * t130; -t7 * t285 - t8 * t284 - t9 * t283 - g(3) * m(4) + (t72 * t298 + t71 * t301 + t70 * t304) * t186 + (t72 * t286 + t71 * t287 + t70 * t288) * t185 + (-t72 * t283 - t71 * t284 - t70 * t285 + m(4)) * t184 + ((-t63 * t262 + t69 * t300) * t186 + (t63 * t268 + t69 * t299) * t185 + (-t198 * t63 - t69 * t226) * t184 - t198 * t3 - t6 * t226) * t132 + ((-t62 * t264 + t68 * t303) * t186 + (t62 * t270 + t68 * t302) * t185 + (-t196 * t62 - t68 * t227) * t184 - t196 * t2 - t5 * t227) * t131 + ((-t61 * t266 + t67 * t306) * t186 + (t61 * t272 + t67 * t305) * t185 + (-t194 * t61 - t67 * t228) * t184 - t194 * t1 - t4 * t228) * t130;];
tauX  = t37;
