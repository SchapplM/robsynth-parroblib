% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V2G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:14:53
% EndTime: 2020-08-06 19:15:03
% DurationCPUTime: 10.75s
% Computational Cost: add. (64614->630), mult. (78054->995), div. (9303->6), fcn. (67950->18), ass. (0->357)
t235 = pkin(5) - pkin(6);
t236 = pkin(2) + pkin(3);
t421 = t235 * t236;
t420 = 2 * pkin(1);
t419 = -0.2e1 * t236 * (-pkin(5) / 0.2e1 + pkin(6) / 0.2e1);
t241 = (qJ(3,1) ^ 2);
t244 = (pkin(2) ^ 2);
t418 = ((t241 - t244) * m(3));
t239 = (qJ(3,2) ^ 2);
t417 = ((t239 - t244) * m(3));
t237 = (qJ(3,3) ^ 2);
t416 = ((t237 - t244) * m(3));
t224 = cos(qJ(2,3));
t204 = t224 ^ 2;
t211 = legFrame(3,3);
t178 = sin(t211);
t181 = cos(t211);
t219 = sin(qJ(1,3));
t225 = cos(qJ(1,3));
t110 = t178 * t219 - t181 * t225;
t111 = t178 * t225 + t181 * t219;
t218 = sin(qJ(2,3));
t347 = qJ(3,3) * t218;
t157 = pkin(1) + t347;
t311 = t224 * t236;
t137 = 0.1e1 / (t157 + t311);
t233 = xDP(2);
t234 = xDP(1);
t91 = (-t110 * t233 - t111 * t234) * t137;
t415 = t204 * t91;
t226 = cos(qJ(2,2));
t206 = t226 ^ 2;
t212 = legFrame(2,3);
t179 = sin(t212);
t182 = cos(t212);
t221 = sin(qJ(1,2));
t227 = cos(qJ(1,2));
t112 = t179 * t221 - t182 * t227;
t113 = t179 * t227 + t182 * t221;
t220 = sin(qJ(2,2));
t350 = qJ(3,2) * t220;
t159 = pkin(1) + t350;
t309 = t236 * t226;
t138 = 0.1e1 / (t159 + t309);
t92 = (-t112 * t233 - t113 * t234) * t138;
t414 = t206 * t92;
t228 = cos(qJ(2,1));
t208 = t228 ^ 2;
t213 = legFrame(1,3);
t180 = sin(t213);
t183 = cos(t213);
t223 = sin(qJ(1,1));
t229 = cos(qJ(1,1));
t114 = t180 * t223 - t183 * t229;
t115 = t180 * t229 + t183 * t223;
t222 = sin(qJ(2,1));
t353 = qJ(3,1) * t222;
t161 = pkin(1) + t353;
t308 = t236 * t228;
t139 = 0.1e1 / (t161 + t308);
t93 = (-t114 * t233 - t115 * t234) * t139;
t413 = t208 * t93;
t232 = xDP(3);
t238 = 0.1e1 / qJ(3,3);
t168 = t235 * t225;
t117 = t157 * t219 - t168;
t165 = t219 * t235;
t275 = t157 * t225 + t165;
t82 = t110 * t311 + t117 * t178 - t275 * t181;
t83 = t111 * t311 + t117 * t181 + t178 * t275;
t64 = (t218 * t232 + (t233 * t83 - t234 * t82) * t224 * t137) * t238;
t34 = t235 * t64;
t412 = t236 * t34;
t240 = 0.1e1 / qJ(3,2);
t169 = t235 * t227;
t119 = t159 * t221 - t169;
t166 = t221 * t235;
t273 = t159 * t227 + t166;
t84 = t112 * t309 + t119 * t179 - t273 * t182;
t85 = t113 * t309 + t119 * t182 + t179 * t273;
t65 = (t220 * t232 + (t233 * t85 - t234 * t84) * t226 * t138) * t240;
t35 = t235 * t65;
t411 = t236 * t35;
t242 = 0.1e1 / qJ(3,1);
t170 = t235 * t229;
t121 = t161 * t223 - t170;
t167 = t223 * t235;
t271 = t161 * t229 + t167;
t86 = t114 * t308 + t121 * t180 - t271 * t183;
t87 = t115 * t308 + t121 * t183 + t180 * t271;
t66 = (t222 * t232 + (t233 * t87 - t234 * t86) * t228 * t139) * t242;
t36 = t235 * t66;
t410 = t236 * t36;
t366 = t236 * t64;
t365 = t236 * t65;
t364 = t236 * t66;
t333 = t137 * t238;
t332 = t138 * t240;
t331 = t139 * t242;
t185 = m(3) * pkin(2) + mrSges(3,1);
t297 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t409 = qJ(3,3) * t185 + t297;
t408 = qJ(3,2) * t185 + t297;
t407 = qJ(3,1) * t185 + t297;
t405 = -2 * pkin(1);
t403 = -0.2e1 * t204;
t402 = -0.2e1 * t206;
t401 = -0.2e1 * t208;
t400 = -3 * t237;
t399 = -3 * t239;
t398 = -3 * t241;
t397 = m(3) * pkin(5);
t396 = m(2) + m(3);
t395 = (mrSges(3,1) * pkin(2));
t175 = m(3) * qJ(3,3) + mrSges(3,3);
t171 = -mrSges(2,2) + t175;
t394 = pkin(1) * t171;
t176 = m(3) * qJ(3,2) + mrSges(3,3);
t172 = -mrSges(2,2) + t176;
t393 = pkin(1) * t172;
t177 = qJ(3,1) * m(3) + mrSges(3,3);
t173 = -mrSges(2,2) + t177;
t392 = pkin(1) * t173;
t174 = mrSges(2,1) + t185;
t391 = g(3) * t174;
t390 = t174 * pkin(1);
t214 = mrSges(3,2) + mrSges(2,3);
t389 = Ifges(2,1) + Ifges(3,1);
t388 = Ifges(3,6) - Ifges(2,6);
t387 = mrSges(3,3) * qJ(3,1);
t386 = mrSges(3,3) * qJ(3,2);
t385 = mrSges(3,3) * qJ(3,3);
t384 = t218 * t91;
t383 = t220 * t92;
t382 = t222 * t93;
t357 = t91 * t235;
t286 = t218 * t357;
t381 = t224 * (t286 - t366);
t380 = t224 * t82;
t379 = t224 * t83;
t195 = 0.2e1 * t385;
t201 = 2 * t395;
t296 = Ifges(3,2) + Ifges(2,3) + t201;
t128 = m(3) * (t237 + t244) + t195 + t296;
t316 = t218 * t236;
t346 = qJ(3,3) * t224;
t146 = t316 - t346;
t97 = (t128 * t218 - t146 * t185) * t238;
t378 = t224 * t97;
t356 = t92 * t235;
t285 = t220 * t356;
t377 = t226 * (t285 - t365);
t376 = t226 * t84;
t375 = t226 * t85;
t197 = 0.2e1 * t386;
t129 = m(3) * (t239 + t244) + t197 + t296;
t314 = t220 * t236;
t349 = qJ(3,2) * t226;
t147 = t314 - t349;
t98 = (t129 * t220 - t147 * t185) * t240;
t374 = t226 * t98;
t355 = t93 * t235;
t284 = t222 * t355;
t373 = t228 * (t284 - t364);
t372 = t228 * t86;
t371 = t228 * t87;
t199 = 0.2e1 * t387;
t130 = m(3) * (t241 + t244) + t199 + t296;
t312 = t222 * t236;
t352 = qJ(3,1) * t228;
t148 = t312 - t352;
t99 = (t130 * t222 - t148 * t185) * t242;
t370 = t228 * t99;
t363 = t236 * t91;
t362 = t236 * t92;
t361 = t236 * t93;
t162 = pkin(1) * t218 + qJ(3,3);
t106 = t162 * t219 - t218 * t168;
t156 = pkin(1) + 0.2e1 * t347;
t116 = t156 * t219 - t168;
t259 = t162 * t225 + t218 * t165;
t276 = t156 * t225 + t165;
t320 = (qJ(3,3) + t236) * (-qJ(3,3) + t236);
t280 = t204 * t320;
t73 = -t110 * t280 - (t178 * t116 - t276 * t181) * t311 - (t178 * t106 - t259 * t181) * qJ(3,3);
t74 = t111 * t280 + (t116 * t181 + t276 * t178) * t311 + (t106 * t181 + t259 * t178) * qJ(3,3);
t46 = t146 * t238 * t232 + (t233 * t74 + t234 * t73) * t333;
t360 = t46 * t175;
t163 = pkin(1) * t220 + qJ(3,2);
t107 = t163 * t221 - t220 * t169;
t158 = pkin(1) + 0.2e1 * t350;
t118 = t158 * t221 - t169;
t258 = t163 * t227 + t220 * t166;
t274 = t158 * t227 + t166;
t319 = (qJ(3,2) + t236) * (-qJ(3,2) + t236);
t279 = t206 * t319;
t75 = -t112 * t279 - (t179 * t118 - t274 * t182) * t309 - (t179 * t107 - t258 * t182) * qJ(3,2);
t76 = t113 * t279 + (t118 * t182 + t274 * t179) * t309 + (t107 * t182 + t258 * t179) * qJ(3,2);
t47 = t147 * t240 * t232 + (t233 * t76 + t234 * t75) * t332;
t359 = t47 * t176;
t164 = pkin(1) * t222 + qJ(3,1);
t108 = t164 * t223 - t222 * t170;
t160 = pkin(1) + 0.2e1 * t353;
t120 = t160 * t223 - t170;
t257 = t164 * t229 + t222 * t167;
t272 = t160 * t229 + t167;
t318 = (qJ(3,1) + t236) * (-qJ(3,1) + t236);
t278 = t208 * t318;
t77 = -t114 * t278 - (t180 * t120 - t272 * t183) * t308 - (t180 * t108 - t257 * t183) * qJ(3,1);
t78 = t115 * t278 + (t120 * t183 + t272 * t180) * t308 + (t108 * t183 + t257 * t180) * qJ(3,1);
t48 = t148 * t242 * t232 + (t233 * t78 + t234 * t77) * t331;
t358 = t48 * t177;
t253 = pkin(2) * mrSges(3,2) + pkin(5) * t174 - Ifges(3,4) - Ifges(2,5);
t268 = -qJ(3,3) * mrSges(3,2) + t388;
t100 = t224 * (pkin(5) * t171 - t268) - t253 * t218;
t345 = t100 * t224;
t269 = -qJ(3,2) * mrSges(3,2) + t388;
t101 = t226 * (pkin(5) * t172 - t269) - t253 * t220;
t344 = t101 * t226;
t270 = -qJ(3,1) * mrSges(3,2) + t388;
t102 = t228 * (pkin(5) * t173 - t270) - t253 * t222;
t343 = t102 * t228;
t216 = xDDP(2);
t342 = t110 * t216;
t217 = xDDP(1);
t341 = t111 * t217;
t340 = t112 * t216;
t339 = t113 * t217;
t338 = t114 * t216;
t337 = t115 * t217;
t336 = t128 * t224;
t335 = t129 * t226;
t334 = t130 * t228;
t245 = pkin(1) ^ 2;
t303 = pkin(5) ^ 2 + t245;
t150 = (-0.2e1 * pkin(5) + pkin(6)) * pkin(6) + t303;
t330 = t150 * t236;
t329 = t171 * t218;
t328 = t172 * t220;
t327 = t173 * t222;
t184 = mrSges(3,2) + t397;
t326 = t184 * t218;
t325 = t184 * t220;
t324 = t184 * t222;
t323 = t185 * t224;
t322 = t185 * t226;
t321 = t185 * t228;
t317 = t218 * t224;
t315 = t220 * t226;
t313 = t222 * t228;
t307 = -t235 ^ 2 - t245;
t301 = 0.2e1 * t236;
t300 = pkin(1) * qJ(3,1) * t93;
t299 = pkin(1) * qJ(3,2) * t92;
t298 = pkin(1) * qJ(3,3) * t91;
t295 = t46 * t333;
t294 = t64 * t333;
t293 = t91 * t333;
t292 = t47 * t332;
t291 = t65 * t332;
t290 = t92 * t332;
t289 = t48 * t331;
t288 = t66 * t331;
t287 = t93 * t331;
t283 = qJ(3,1) * t208 * t235;
t282 = qJ(3,2) * t206 * t235;
t281 = qJ(3,3) * t204 * t235;
t277 = (-Ifges(2,2) - Ifges(3,3) + t389);
t202 = -2 * t395;
t267 = t202 + t277;
t22 = t48 - t364;
t21 = t22 * t236 - t241 * t66;
t23 = t46 - t366;
t19 = t23 * t236 - t237 * t64;
t24 = t47 - t365;
t20 = t236 * t24 - t239 * t65;
t131 = -g(1) * t178 + g(2) * t181;
t134 = g(1) * t181 + g(2) * t178;
t265 = t131 * t219 + t134 * t225;
t132 = -g(1) * t179 + g(2) * t182;
t135 = g(1) * t182 + g(2) * t179;
t264 = t132 * t221 + t135 * t227;
t133 = -g(1) * t180 + g(2) * t183;
t136 = g(1) * t183 + g(2) * t180;
t263 = t133 * t223 + t136 * t229;
t262 = t277 + t416;
t261 = t277 + t417;
t260 = t277 + t418;
t155 = t396 * pkin(1) + mrSges(1,1);
t256 = t174 * t224 + t155 + t329;
t255 = t174 * t226 + t155 + t328;
t254 = t174 * t228 + t155 + t327;
t252 = t303 * m(2) + 0.2e1 * pkin(5) * t214 + Ifges(1,3) + t389;
t231 = mrSges(2,2) * pkin(5);
t215 = xDDP(3);
t210 = t236 ^ 2;
t207 = t228 * t208;
t205 = t226 * t206;
t203 = t224 * t204;
t154 = 0.2e1 * t390;
t149 = t396 * pkin(5) - mrSges(1,2) + t214;
t109 = (-mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1) * pkin(5) + Ifges(3,4) / 0.2e1 + Ifges(2,5) / 0.2e1 + (-t397 / 0.2e1 - mrSges(3,2) / 0.2e1) * pkin(2);
t105 = (m(3) * t148 - t185 * t222) * t242;
t104 = (m(3) * t147 - t185 * t220) * t240;
t103 = (m(3) * t146 - t185 * t218) * t238;
t90 = t93 ^ 2;
t89 = t92 ^ 2;
t88 = t91 ^ 2;
t81 = (t201 - t260 - 0.2e1 * t387) * t208 + t154 * t228 + t327 * t420 + (t241 + t303) * m(3) + t199 + 0.2e1 * t407 * t313 + t252;
t80 = (t201 - t261 - 0.2e1 * t386) * t206 + t154 * t226 + t328 * t420 + (t239 + t303) * m(3) + t197 + 0.2e1 * t408 * t315 + t252;
t79 = (t201 - t262 - 0.2e1 * t385) * t204 + t154 * t224 + t329 * t420 + (t237 + t303) * m(3) + t195 + 0.2e1 * t409 * t317 + t252;
t63 = t66 ^ 2;
t62 = t65 ^ 2;
t61 = t64 ^ 2;
t60 = t407 * t66;
t59 = t408 * t65;
t58 = t409 * t64;
t54 = (-t115 * t324 + (m(3) * t77 + t86 * t321) * t242) * t139;
t53 = (-t114 * t324 + (m(3) * t78 - t87 * t321) * t242) * t139;
t52 = (-t113 * t325 + (m(3) * t75 + t84 * t322) * t240) * t138;
t51 = (-t112 * t325 + (m(3) * t76 - t85 * t322) * t240) * t138;
t50 = (-t111 * t326 + (m(3) * t73 + t82 * t323) * t238) * t137;
t49 = (-t110 * t326 + (m(3) * t74 - t83 * t323) * t238) * t137;
t45 = (-t102 * t115 + (-t185 * t77 - t86 * t334) * t242) * t139;
t44 = (-t102 * t114 + (-t185 * t78 + t87 * t334) * t242) * t139;
t43 = (-t101 * t113 + (-t185 * t75 - t84 * t335) * t240) * t138;
t42 = (-t101 * t112 + (-t185 * t76 + t85 * t335) * t240) * t138;
t41 = (-t100 * t111 + (-t185 * t73 - t82 * t336) * t238) * t137;
t40 = (-t100 * t110 + (-t185 * t74 + t83 * t336) * t238) * t137;
t39 = t48 * t235;
t38 = t47 * t235;
t37 = t46 * t235;
t30 = (-t115 * t81 + (t77 * t324 - t86 * t343) * t242) * t139;
t29 = (-t114 * t81 + (t78 * t324 + t87 * t343) * t242) * t139;
t28 = (-t113 * t80 + (t75 * t325 - t84 * t344) * t240) * t138;
t27 = (-t112 * t80 + (t76 * t325 + t85 * t344) * t240) * t138;
t26 = (-t111 * t79 + (t73 * t326 - t82 * t345) * t238) * t137;
t25 = (-t110 * t79 + (t74 * t326 + t83 * t345) * t238) * t137;
t18 = (t355 + (-t148 + t352) * t66 + (0.2e1 * t48 - t364) * t222) * t93 * t139;
t17 = (t356 + (-t147 + t349) * t65 + (0.2e1 * t47 - t365) * t220) * t92 * t138;
t16 = (t357 + (-t146 + t346) * t64 + (0.2e1 * t46 - t366) * t218) * t91 * t137;
t15 = (-t66 * t283 + (t421 * t66 - t39) * t313 + (-t207 * t318 + (t353 * t405 - t241 + t307) * t228 - t161 * t208 * t301) * t93) * t287 + (-t93 * t283 + (t22 + t284) * t308 + t22 * t161) * t288 + (t161 * t66 - t373) * t289;
t14 = (-t65 * t282 + (t421 * t65 - t38) * t315 + (-t205 * t319 + (t350 * t405 - t239 + t307) * t226 - t159 * t206 * t301) * t92) * t290 + (-t92 * t282 + (t24 + t285) * t309 + t24 * t159) * t291 + (t159 * t65 - t377) * t292;
t13 = (-t64 * t281 + (t421 * t64 - t37) * t317 + (-t203 * t320 + (t347 * t405 - t237 + t307) * t224 - t157 * t204 * t301) * t91) * t293 + (-t91 * t281 + (t23 + t286) * t311 + t23 * t157) * t294 + (t157 * t64 - t381) * t295;
t12 = (-(t210 + t398) * t207 * t361 + (-0.3e1 * (-t241 / 0.3e1 + t210) * t353 + (t241 - t210) * t420) * t413 + ((-t36 * t241 + (-0.4e1 * t300 - t39 + t410) * t236) * t222 + t361 * t398 - t93 * t330) * t228) * t287 + ((t21 * t236 + t284 * t318) * t228 + t21 * pkin(1)) * t288 + (pkin(1) * t66 - t373) * t236 * t289 + (((t39 - 0.2e1 * t410) * t208 + (-t150 - t241) * t382 - 0.2e1 * t300 - t39 + t66 * t419) * t287 + (t21 * t222 + (t401 + 0.1e1) * t93 * t421) * t288 + (t66 * t312 + (t208 - 0.1e1) * t355) * t289) * qJ(3,1);
t11 = (-(t210 + t399) * t205 * t362 + (-0.3e1 * (-t239 / 0.3e1 + t210) * t350 + (t239 - t210) * t420) * t414 + ((-t35 * t239 + (-0.4e1 * t299 - t38 + t411) * t236) * t220 + t362 * t399 - t92 * t330) * t226) * t290 + ((t20 * t236 + t285 * t319) * t226 + t20 * pkin(1)) * t291 + (pkin(1) * t65 - t377) * t236 * t292 + (((t38 - 0.2e1 * t411) * t206 + (-t150 - t239) * t383 - 0.2e1 * t299 - t38 + t65 * t419) * t290 + (t20 * t220 + (t402 + 0.1e1) * t92 * t421) * t291 + (t65 * t314 + (t206 - 0.1e1) * t356) * t292) * qJ(3,2);
t10 = (-(t210 + t400) * t203 * t363 + (-0.3e1 * (-t237 / 0.3e1 + t210) * t347 + (t237 - t210) * t420) * t415 + ((-t34 * t237 + (-0.4e1 * t298 - t37 + t412) * t236) * t218 + t363 * t400 - t91 * t330) * t224) * t293 + ((t19 * t236 + t286 * t320) * t224 + t19 * pkin(1)) * t294 + (pkin(1) * t64 - t381) * t236 * t295 + (((t37 - 0.2e1 * t412) * t204 + (-t150 - t237) * t384 - 0.2e1 * t298 - t37 + t64 * t419) * t293 + (t19 * t218 + (t403 + 0.1e1) * t91 * t421) * t294 + (t64 * t316 + (t204 - 0.1e1) * t357) * t295) * qJ(3,3);
t9 = t185 * t15 + (t208 * t90 - t63 - t90) * t177 + (t228 * g(3) - t12) * m(3) + (-t90 * t321 - t184 * t18 + (-pkin(1) * t90 - t263) * m(3)) * t222;
t8 = t185 * t14 + (t206 * t89 - t62 - t89) * t176 + (t226 * g(3) - t11) * m(3) + (-t89 * t322 - t184 * t17 + (-pkin(1) * t89 - t264) * m(3)) * t220;
t7 = t185 * t13 + (t204 * t88 - t61 - t88) * t175 + (t224 * g(3) - t10) * m(3) + (-t88 * t323 - t184 * t16 + (-pkin(1) * t88 - t265) * m(3)) * t218;
t6 = -t102 * t18 - t130 * t15 + t185 * t12 + 0.2e1 * t66 * t358 + (-t263 * t173 - t391) * t228 + (-g(3) * t173 + t263 * t174) * t222 + (t407 * t401 - ((t199 + t267 + t418) * t222 + t392) * t228 + t222 * t390 + t407) * t90;
t5 = -t101 * t17 - t129 * t14 + t185 * t11 + 0.2e1 * t65 * t359 + (-t264 * t172 - t391) * t226 + (-g(3) * t172 + t264 * t174) * t220 + (t408 * t402 - ((t197 + t267 + t417) * t220 + t393) * t226 + t220 * t390 + t408) * t89;
t4 = -t100 * t16 - t128 * t13 + t185 * t10 + 0.2e1 * t64 * t360 + (-t265 * t171 - t391) * t224 + (-g(3) * t171 + t265 * t174) * t218 + (t409 * t403 - ((t195 + t267 + t416) * t218 + t394) * t224 + t218 * t390 + t409) * t88;
t3 = -t81 * t18 - t102 * t15 - t12 * t324 + 0.4e1 * (t60 - t358 / 0.2e1) * t413 + 0.2e1 * (((t199 + t202 + t260) * t66 + t48 * t185) * t382 + t66 * (t109 * t66 + t48 * t184 + t93 * t392)) * t228 + ((-t177 * pkin(5) + t231 + t270) * t63 + (m(3) * t48 - t174 * t66) * t93 * t420) * t222 - 0.2e1 * (t60 - t358) * t93 + (-t254 * t133 - t136 * t149) * t229 + t223 * (-t133 * t149 + t254 * t136);
t2 = -t80 * t17 - t101 * t14 - t11 * t325 + 0.4e1 * (t59 - t359 / 0.2e1) * t414 + 0.2e1 * (((t197 + t202 + t261) * t65 + t47 * t185) * t383 + t65 * (t109 * t65 + t47 * t184 + t92 * t393)) * t226 + ((-t176 * pkin(5) + t231 + t269) * t62 + (m(3) * t47 - t174 * t65) * t92 * t420) * t220 - 0.2e1 * (t59 - t359) * t92 + (-t255 * t132 - t135 * t149) * t227 + t221 * (-t132 * t149 + t255 * t135);
t1 = -t79 * t16 - t100 * t13 - t10 * t326 + 0.4e1 * (t58 - t360 / 0.2e1) * t415 + 0.2e1 * (((t195 + t202 + t262) * t64 + t46 * t185) * t384 + t64 * (t109 * t64 + t46 * t184 + t91 * t394)) * t224 + ((-t175 * pkin(5) + t231 + t268) * t61 + (m(3) * t46 - t174 * t64) * t91 * t420) * t218 - 0.2e1 * (t58 - t360) * t91 + (-t256 * t131 - t134 * t149) * t225 + t219 * (-t131 * t149 + t256 * t134);
t31 = [(-g(1) + t217) * m(4) + ((t148 * t54 + t222 * t45) * t242 + (t147 * t52 + t220 * t43) * t240 + (t146 * t50 + t218 * t41) * t238) * t215 + (-t30 * t338 + (-t30 * t217 - t3) * t115 + ((-t45 * t372 + t54 * t77) * t217 + (t45 * t371 + t54 * t78) * t216 - t6 * t372 + t77 * t9) * t242) * t139 + (-t28 * t340 + (-t28 * t217 - t2) * t113 + ((-t43 * t376 + t52 * t75) * t217 + (t43 * t375 + t52 * t76) * t216 - t5 * t376 + t75 * t8) * t240) * t138 + (-t26 * t342 + (-t26 * t217 - t1) * t111 + ((-t41 * t380 + t50 * t73) * t217 + (t41 * t379 + t50 * t74) * t216 - t4 * t380 + t73 * t7) * t238) * t137; (-g(2) + t216) * m(4) + ((t148 * t53 + t222 * t44) * t242 + (t147 * t51 + t220 * t42) * t240 + (t146 * t49 + t218 * t40) * t238) * t215 + (-t29 * t337 + (-t29 * t216 - t3) * t114 + ((-t44 * t372 + t53 * t77) * t217 + (t44 * t371 + t53 * t78) * t216 + t6 * t371 + t78 * t9) * t242) * t139 + (-t27 * t339 + (-t27 * t216 - t2) * t112 + ((-t42 * t376 + t51 * t75) * t217 + (t42 * t375 + t51 * t76) * t216 + t5 * t375 + t76 * t8) * t240) * t138 + (-t25 * t341 + (-t25 * t216 - t1) * t110 + ((-t40 * t380 + t49 * t73) * t217 + (t40 * t379 + t49 * t74) * t216 + t4 * t379 + t74 * t7) * t238) * t137; -g(3) * m(4) + (t148 * t9 + t222 * t6) * t242 + (t147 * t8 + t220 * t5) * t240 + (t146 * t7 + t218 * t4) * t238 + (m(4) + (t105 * t148 + t222 * t99) * t242 + (t104 * t147 + t220 * t98) * t240 + (t103 * t146 + t218 * t97) * t238) * t215 + ((-t337 - t338) * (t148 * t184 + t102) * t222 + (t105 * t77 - t86 * t370) * t217 + (t105 * t78 + t87 * t370) * t216) * t331 + ((-t339 - t340) * (t147 * t184 + t101) * t220 + (t104 * t75 - t84 * t374) * t217 + (t104 * t76 + t85 * t374) * t216) * t332 + ((-t341 - t342) * (t146 * t184 + t100) * t218 + (t103 * t73 - t82 * t378) * t217 + (t103 * t74 + t83 * t378) * t216) * t333;];
tauX  = t31;
