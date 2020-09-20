% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V2G2A0
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
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:20:45
% EndTime: 2020-08-06 19:20:56
% DurationCPUTime: 11.09s
% Computational Cost: add. (86121->655), mult. (101859->1049), div. (11898->6), fcn. (76671->18), ass. (0->365)
t213 = sin(qJ(2,3));
t396 = pkin(1) * t213;
t157 = qJ(3,3) + t396;
t214 = sin(qJ(1,3));
t230 = pkin(5) - pkin(6);
t160 = t214 * t230;
t170 = t213 * qJ(3,3);
t220 = cos(qJ(1,3));
t219 = cos(qJ(2,3));
t201 = t219 ^ 2;
t231 = pkin(2) + pkin(3);
t322 = (qJ(3,3) + t231) * (-qJ(3,3) + t231);
t277 = t201 * t322;
t303 = t231 * t219;
t103 = t220 * t277 + ((0.2e1 * t170 + pkin(1)) * t220 + t160) * t303 + qJ(3,3) * (t157 * t220 + t213 * t160);
t227 = xDP(3);
t228 = xDP(2);
t229 = xDP(1);
t154 = t170 + pkin(1);
t419 = t154 + t303;
t131 = 0.1e1 / t419;
t233 = 0.1e1 / qJ(3,3);
t358 = t131 * t233;
t263 = pkin(1) * t214 - t230 * t220;
t317 = t214 * qJ(3,3);
t112 = t263 * t213 + t317;
t268 = t213 * t317;
t115 = t263 + 0.2e1 * t268;
t276 = t213 * t322;
t400 = pkin(1) * qJ(3,3);
t128 = -t276 + t400;
t207 = legFrame(3,2);
t173 = sin(t207);
t176 = cos(t207);
t275 = t214 * t322;
t406 = -0.2e1 * t231;
t292 = qJ(3,3) * t406;
t337 = t176 * t231;
t349 = t173 * t231;
t352 = t173 * qJ(3,3);
t76 = (-t173 * t275 + t176 * t292) * t201 + (-t115 * t349 - t176 * t128) * t219 - t112 * t352 + t157 * t337;
t340 = t176 * qJ(3,3);
t77 = (t173 * t292 + t176 * t275) * t201 + (t115 * t337 - t173 * t128) * t219 + t112 * t340 + t157 * t349;
t34 = (t103 * t227 + t228 * t76 + t229 * t77) * t358;
t109 = t419 * t220 + t160;
t361 = t109 * t219;
t116 = t263 + t268;
t315 = t214 * t231;
t318 = t213 * t231;
t88 = (-t173 * t315 - t340) * t201 + (-t116 * t173 + t176 * t318) * t219 + t176 * t157;
t89 = (t176 * t315 - t352) * t201 + (t116 * t176 + t173 * t318) * t219 + t173 * t157;
t61 = (t227 * t361 + t228 * t88 + t229 * t89) * t358;
t383 = t231 * t61;
t26 = t34 - t383;
t215 = sin(qJ(2,2));
t395 = pkin(1) * t215;
t158 = qJ(3,2) + t395;
t216 = sin(qJ(1,2));
t161 = t216 * t230;
t171 = t215 * qJ(3,2);
t222 = cos(qJ(1,2));
t221 = cos(qJ(2,2));
t202 = t221 ^ 2;
t321 = (qJ(3,2) + t231) * (-qJ(3,2) + t231);
t274 = t202 * t321;
t302 = t231 * t221;
t104 = t222 * t274 + ((0.2e1 * t171 + pkin(1)) * t222 + t161) * t302 + qJ(3,2) * (t158 * t222 + t215 * t161);
t155 = t171 + pkin(1);
t418 = t155 + t302;
t132 = 0.1e1 / t418;
t235 = 0.1e1 / qJ(3,2);
t357 = t132 * t235;
t262 = pkin(1) * t216 - t230 * t222;
t312 = t216 * qJ(3,2);
t113 = t262 * t215 + t312;
t267 = t215 * t312;
t117 = t262 + 0.2e1 * t267;
t273 = t215 * t321;
t401 = pkin(1) * qJ(3,2);
t129 = -t273 + t401;
t208 = legFrame(2,2);
t174 = sin(t208);
t177 = cos(t208);
t272 = t216 * t321;
t293 = qJ(3,2) * t406;
t333 = t177 * t231;
t345 = t174 * t231;
t348 = t174 * qJ(3,2);
t78 = (-t174 * t272 + t177 * t293) * t202 + (-t117 * t345 - t177 * t129) * t221 - t113 * t348 + t158 * t333;
t336 = t177 * qJ(3,2);
t79 = (t174 * t293 + t177 * t272) * t202 + (t117 * t333 - t174 * t129) * t221 + t113 * t336 + t158 * t345;
t35 = (t104 * t227 + t228 * t78 + t229 * t79) * t357;
t110 = t418 * t222 + t161;
t360 = t110 * t221;
t118 = t262 + t267;
t310 = t216 * t231;
t313 = t215 * t231;
t90 = (-t174 * t310 - t336) * t202 + (-t118 * t174 + t177 * t313) * t221 + t177 * t158;
t91 = (t177 * t310 - t348) * t202 + (t118 * t177 + t174 * t313) * t221 + t174 * t158;
t62 = (t227 * t360 + t228 * t90 + t229 * t91) * t357;
t382 = t231 * t62;
t27 = t35 - t382;
t217 = sin(qJ(2,1));
t394 = pkin(1) * t217;
t159 = qJ(3,1) + t394;
t218 = sin(qJ(1,1));
t162 = t218 * t230;
t172 = t217 * qJ(3,1);
t224 = cos(qJ(1,1));
t223 = cos(qJ(2,1));
t203 = t223 ^ 2;
t320 = (qJ(3,1) + t231) * (-qJ(3,1) + t231);
t271 = t203 * t320;
t301 = t231 * t223;
t105 = t224 * t271 + ((0.2e1 * t172 + pkin(1)) * t224 + t162) * t301 + qJ(3,1) * (t159 * t224 + t217 * t162);
t156 = t172 + pkin(1);
t417 = t156 + t301;
t133 = 0.1e1 / t417;
t237 = 0.1e1 / qJ(3,1);
t356 = t133 * t237;
t261 = pkin(1) * t218 - t230 * t224;
t307 = t218 * qJ(3,1);
t114 = t261 * t217 + t307;
t266 = t217 * t307;
t119 = t261 + 0.2e1 * t266;
t270 = t217 * t320;
t402 = pkin(1) * qJ(3,1);
t130 = -t270 + t402;
t209 = legFrame(1,2);
t175 = sin(t209);
t178 = cos(t209);
t269 = t218 * t320;
t294 = qJ(3,1) * t406;
t329 = t178 * t231;
t341 = t175 * t231;
t344 = t175 * qJ(3,1);
t80 = (-t175 * t269 + t178 * t294) * t203 + (-t119 * t341 - t178 * t130) * t223 - t114 * t344 + t159 * t329;
t332 = t178 * qJ(3,1);
t81 = (t175 * t294 + t178 * t269) * t203 + (t119 * t329 - t175 * t130) * t223 + t114 * t332 + t159 * t341;
t36 = (t105 * t227 + t228 * t80 + t229 * t81) * t356;
t111 = t417 * t224 + t162;
t359 = t111 * t223;
t120 = t261 + t266;
t305 = t218 * t231;
t308 = t217 * t231;
t92 = (-t175 * t305 - t332) * t203 + (-t120 * t175 + t178 * t308) * t223 + t178 * t159;
t93 = (t178 * t305 - t344) * t203 + (t120 * t178 + t175 * t308) * t223 + t175 * t159;
t63 = (t227 * t359 + t228 * t92 + t229 * t93) * t356;
t381 = t231 * t63;
t25 = t36 - t381;
t424 = 0.2e1 * pkin(1);
t423 = 0.2e1 * t231;
t236 = qJ(3,1) ^ 2;
t240 = pkin(2) ^ 2;
t422 = (t236 - t240) * m(3);
t234 = qJ(3,2) ^ 2;
t421 = (t234 - t240) * m(3);
t232 = qJ(3,3) ^ 2;
t420 = (t232 - t240) * m(3);
t180 = m(3) * pkin(2) + mrSges(3,1);
t288 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t416 = qJ(3,3) * t180 + t288;
t415 = qJ(3,2) * t180 + t288;
t414 = qJ(3,1) * t180 + t288;
t412 = -0.2e1 * t201;
t411 = -0.2e1 * t202;
t410 = -0.2e1 * t203;
t409 = 0.2e1 * t219;
t408 = 0.2e1 * t221;
t407 = 0.2e1 * t223;
t405 = m(3) * pkin(5);
t404 = m(2) + m(3);
t403 = mrSges(3,1) * pkin(2);
t167 = m(3) * qJ(3,3) + mrSges(3,3);
t163 = -mrSges(2,2) + t167;
t399 = pkin(1) * t163;
t168 = m(3) * qJ(3,2) + mrSges(3,3);
t164 = -mrSges(2,2) + t168;
t398 = pkin(1) * t164;
t169 = qJ(3,1) * m(3) + mrSges(3,3);
t165 = -mrSges(2,2) + t169;
t397 = pkin(1) * t165;
t206 = mrSges(3,2) + mrSges(2,3);
t393 = Ifges(2,1) + Ifges(3,1);
t392 = Ifges(3,6) - Ifges(2,6);
t23 = t230 * t26;
t100 = (-t214 * t227 + (-t173 * t228 + t176 * t229) * t220) * t131;
t289 = t100 * t400;
t391 = 0.2e1 * t289 + t23;
t24 = t230 * t27;
t101 = (-t216 * t227 + (-t174 * t228 + t177 * t229) * t222) * t132;
t290 = t101 * t401;
t390 = 0.2e1 * t290 + t24;
t22 = t230 * t25;
t102 = (-t218 * t227 + (-t175 * t228 + t178 * t229) * t224) * t133;
t291 = t102 * t402;
t389 = 0.2e1 * t291 + t22;
t388 = mrSges(3,3) * qJ(3,1);
t387 = mrSges(3,3) * qJ(3,2);
t386 = mrSges(3,3) * qJ(3,3);
t319 = t213 * t230;
t283 = t100 * t319;
t385 = t219 * (t283 - t383);
t314 = t215 * t230;
t282 = t101 * t314;
t384 = t221 * (t282 - t382);
t380 = t232 * t61;
t379 = t234 * t62;
t378 = t236 * t63;
t377 = t34 * t167;
t376 = t35 * t168;
t375 = t36 * t169;
t309 = t217 * t230;
t281 = t102 * t309;
t374 = (t281 - t381) * t223;
t370 = t100 * t201;
t369 = t100 * t213;
t368 = t100 * t230;
t367 = t101 * t202;
t366 = t101 * t215;
t365 = t101 * t230;
t364 = t102 * t203;
t363 = t102 * t217;
t362 = t102 * t230;
t355 = t163 * t213;
t354 = t164 * t215;
t353 = t165 * t217;
t211 = xDDP(2);
t351 = t173 * t211;
t350 = t173 * t220;
t347 = t174 * t211;
t346 = t174 * t222;
t343 = t175 * t211;
t342 = t175 * t224;
t212 = xDDP(1);
t339 = t176 * t212;
t338 = t176 * t220;
t335 = t177 * t212;
t334 = t177 * t222;
t331 = t178 * t212;
t330 = t178 * t224;
t179 = mrSges(3,2) + t405;
t328 = t179 * t213;
t327 = t179 * t215;
t326 = t179 * t217;
t325 = t180 * t219;
t324 = t180 * t221;
t323 = t180 * t223;
t210 = xDDP(3);
t316 = t214 * t210;
t311 = t216 * t210;
t306 = t218 * t210;
t304 = t230 * t231;
t297 = pkin(1) ^ 2 + pkin(5) ^ 2;
t199 = 0.2e1 * t403;
t287 = Ifges(3,2) + Ifges(2,3) + t199;
t286 = t219 * qJ(3,3) * t61;
t285 = t221 * qJ(3,2) * t62;
t284 = t63 * t223 * qJ(3,1);
t280 = t220 * t328;
t279 = t222 * t327;
t278 = t224 * t326;
t265 = -Ifges(2,2) - Ifges(3,3) + t393;
t260 = t297 + (-0.2e1 * pkin(5) + pkin(6)) * pkin(6);
t259 = -qJ(3,1) * mrSges(3,2) + t392;
t258 = -qJ(3,2) * mrSges(3,2) + t392;
t257 = -qJ(3,3) * mrSges(3,2) + t392;
t200 = -0.2e1 * t403;
t256 = t200 + t265;
t137 = g(1) * t176 - g(2) * t173;
t255 = g(3) * t220 + t137 * t214;
t138 = g(1) * t177 - g(2) * t174;
t254 = g(3) * t222 + t138 * t216;
t139 = g(1) * t178 - g(2) * t175;
t253 = g(3) * t224 + t139 * t218;
t19 = t231 * t26 - t380;
t20 = t231 * t27 - t379;
t21 = t231 * t25 - t378;
t166 = mrSges(2,1) + t180;
t252 = t166 * t219 + t355;
t251 = t166 * t221 + t354;
t250 = t166 * t223 + t353;
t249 = t265 + t420;
t248 = t265 + t421;
t247 = t265 + t422;
t246 = pkin(2) * mrSges(3,2) + pkin(5) * t166 - Ifges(3,4) - Ifges(2,5);
t245 = t297 * m(2) + 0.2e1 * t206 * pkin(5) + Ifges(1,3) + t393;
t226 = mrSges(2,2) * pkin(5);
t205 = t231 ^ 2;
t197 = 0.2e1 * t388;
t195 = 0.2e1 * t387;
t193 = 0.2e1 * t386;
t153 = t404 * pkin(1) + mrSges(1,1);
t152 = t166 * t424;
t151 = t153 * g(3);
t150 = t404 * pkin(5) - mrSges(1,2) + t206;
t149 = t236 + t260;
t148 = t234 + t260;
t147 = t232 + t260;
t146 = t150 * g(3);
t136 = g(1) * t175 + g(2) * t178;
t135 = g(1) * t174 + g(2) * t177;
t134 = g(1) * t173 + g(2) * t176;
t127 = m(3) * (t236 + t240) + t197 + t287;
t126 = m(3) * (t234 + t240) + t195 + t287;
t125 = m(3) * (t232 + t240) + t193 + t287;
t121 = (-mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1) * pkin(5) + Ifges(3,4) / 0.2e1 + Ifges(2,5) / 0.2e1 + (-t405 / 0.2e1 - mrSges(3,2) / 0.2e1) * pkin(2);
t108 = (pkin(5) * t165 - t259) * t223 - t217 * t246;
t107 = (pkin(5) * t164 - t258) * t221 - t215 * t246;
t106 = (pkin(5) * t163 - t257) * t219 - t213 * t246;
t99 = t102 ^ 2;
t98 = t101 ^ 2;
t97 = t100 ^ 2;
t84 = (t199 - t247 - 0.2e1 * t388) * t203 + t152 * t223 + t353 * t424 + (t236 + t297) * m(3) + t197 + t414 * t217 * t407 + t245;
t83 = (t199 - t248 - 0.2e1 * t387) * t202 + t152 * t221 + t354 * t424 + (t234 + t297) * m(3) + t195 + t415 * t215 * t408 + t245;
t82 = (t199 - t249 - 0.2e1 * t386) * t201 + t152 * t219 + t355 * t424 + (t232 + t297) * m(3) + t193 + t416 * t213 * t409 + t245;
t69 = (-t218 * t326 + (m(3) * t105 - t111 * t323) * t237) * t133;
t68 = (-t216 * t327 + (m(3) * t104 - t110 * t324) * t235) * t132;
t67 = (-t214 * t328 + (m(3) * t103 - t109 * t325) * t233) * t131;
t66 = (-t108 * t218 + (-t105 * t180 + t127 * t359) * t237) * t133;
t65 = (-t107 * t216 + (-t104 * t180 + t126 * t360) * t235) * t132;
t64 = (-t106 * t214 + (-t103 * t180 + t125 * t361) * t233) * t131;
t60 = t63 ^ 2;
t59 = t62 ^ 2;
t58 = t61 ^ 2;
t57 = t414 * t63;
t56 = t415 * t62;
t55 = t416 * t61;
t54 = (t178 * t278 + (m(3) * t81 - t180 * t93) * t237) * t133;
t53 = (t177 * t279 + (m(3) * t79 - t180 * t91) * t235) * t132;
t52 = (t176 * t280 + (m(3) * t77 - t180 * t89) * t233) * t131;
t51 = (-t175 * t278 + (m(3) * t80 - t180 * t92) * t237) * t133;
t50 = (-t174 * t279 + (m(3) * t78 - t180 * t90) * t235) * t132;
t49 = (-t173 * t280 + (m(3) * t76 - t180 * t88) * t233) * t131;
t45 = (t108 * t330 + (t127 * t93 - t180 * t81) * t237) * t133;
t44 = (t107 * t334 + (t126 * t91 - t180 * t79) * t235) * t132;
t43 = (t106 * t338 + (t125 * t89 - t180 * t77) * t233) * t131;
t42 = (-t108 * t342 + (t127 * t92 - t180 * t80) * t237) * t133;
t41 = (-t107 * t346 + (t126 * t90 - t180 * t78) * t235) * t132;
t40 = (-t106 * t350 + (t125 * t88 - t180 * t76) * t233) * t131;
t33 = (t84 * t330 + (t108 * t93 + t81 * t326) * t237) * t133;
t32 = (t83 * t334 + (t107 * t91 + t79 * t327) * t235) * t132;
t31 = (t82 * t338 + (t106 * t89 + t77 * t328) * t233) * t131;
t30 = (-t84 * t342 + (t108 * t92 + t80 * t326) * t237) * t133;
t29 = (-t83 * t346 + (t107 * t90 + t78 * t327) * t235) * t132;
t28 = (-t82 * t350 + (t106 * t88 + t76 * t328) * t233) * t131;
t18 = (0.2e1 * t25 * t217 + 0.2e1 * t284 + t362) * t102 * t133;
t17 = (0.2e1 * t27 * t215 + 0.2e1 * t285 + t365) * t101 * t132;
t16 = (0.2e1 * t26 * t213 + 0.2e1 * t286 + t368) * t100 * t131;
t15 = (-t223 * (t230 * t284 + t389 * t217 + (t156 * t223 * t423 + t149 + t271) * t102) * t102 + (-qJ(3,1) * t203 * t362 + (t25 + t281) * t301 + t25 * t156) * t63 + (t63 * t156 - t374) * t36) * t356;
t14 = (-t221 * (t230 * t285 + t390 * t215 + (t155 * t221 * t423 + t148 + t274) * t101) * t101 + (-qJ(3,2) * t202 * t365 + (t27 + t282) * t302 + t27 * t155) * t62 + (t155 * t62 - t384) * t35) * t357;
t13 = (-t219 * (t230 * t286 + t391 * t213 + (t154 * t219 * t423 + t147 + t277) * t100) * t100 + (-qJ(3,3) * t201 * t368 + (t26 + t283) * t303 + t26 * t154) * t61 + (t154 * t61 - t385) * t34) * t358;
t12 = ((-(t205 - 0.3e1 * t236) * t301 * t364 + (t230 * (-t423 * t63 + t36) * qJ(3,1) + (-0.3e1 * (-t236 / 0.3e1 + t205) * t172 + (t236 - t205) * t424) * t102) * t203 + (-t309 * t378 + ((-0.4e1 * t291 - t22) * t217 - t102 * (0.3e1 * t236 + t260)) * t231) * t223 - (t149 * t363 + t389) * qJ(3,1)) * t102 + ((t231 * t21 + t270 * t362) * t223 + t21 * pkin(1) + (t21 * t217 + (t410 + 0.1e1) * t102 * t304) * qJ(3,1)) * t63 + ((pkin(1) * t63 - t374) * t231 + (t63 * t308 + (t203 - 0.1e1) * t362) * qJ(3,1)) * t36) * t356;
t11 = ((-(t205 - 0.3e1 * t234) * t302 * t367 + (t230 * (-t423 * t62 + t35) * qJ(3,2) + (-0.3e1 * (-t234 / 0.3e1 + t205) * t171 + (t234 - t205) * t424) * t101) * t202 + (-t314 * t379 + ((-0.4e1 * t290 - t24) * t215 - t101 * (0.3e1 * t234 + t260)) * t231) * t221 - (t148 * t366 + t390) * qJ(3,2)) * t101 + ((t231 * t20 + t273 * t365) * t221 + t20 * pkin(1) + (t20 * t215 + (t411 + 0.1e1) * t101 * t304) * qJ(3,2)) * t62 + ((pkin(1) * t62 - t384) * t231 + (t62 * t313 + (t202 - 0.1e1) * t365) * qJ(3,2)) * t35) * t357;
t10 = ((-(t205 - 0.3e1 * t232) * t303 * t370 + (t230 * (-t423 * t61 + t34) * qJ(3,3) + (-0.3e1 * (-t232 / 0.3e1 + t205) * t170 + (t232 - t205) * t424) * t100) * t201 + (-t319 * t380 + ((-0.4e1 * t289 - t23) * t213 - t100 * (0.3e1 * t232 + t260)) * t231) * t219 - (t147 * t369 + t391) * qJ(3,3)) * t100 + ((t231 * t19 + t276 * t368) * t219 + t19 * pkin(1) + (t19 * t213 + (t412 + 0.1e1) * t100 * t304) * qJ(3,3)) * t61 + ((pkin(1) * t61 - t385) * t231 + (t61 * t318 + (t201 - 0.1e1) * t368) * qJ(3,3)) * t34) * t358;
t9 = t180 * t15 + (t99 * t203 - t60 - t99) * t169 + (t223 * t136 - t12) * m(3) + (-t99 * t323 - t179 * t18 + (-pkin(1) * t99 - t253) * m(3)) * t217;
t8 = t180 * t14 + (t98 * t202 - t59 - t98) * t168 + (t221 * t135 - t11) * m(3) + (-t98 * t324 - t179 * t17 + (-pkin(1) * t98 - t254) * m(3)) * t215;
t7 = t180 * t13 + (t97 * t201 - t58 - t97) * t167 + (t134 * t219 - t10) * m(3) + (-t97 * t325 - t179 * t16 + (-pkin(1) * t97 - t255) * m(3)) * t213;
t6 = -t108 * t18 - t127 * t15 + t180 * t12 + 0.2e1 * t63 * t375 + (-t136 * t166 - t253 * t165) * t223 + t217 * (-t136 * t165 + t253 * t166) + (t414 * t410 - ((t197 + t256 + t422) * t217 + t397) * t223 + t166 * t394 + t414) * t99;
t5 = -t107 * t17 - t126 * t14 + t180 * t11 + 0.2e1 * t62 * t376 + (-t135 * t166 - t254 * t164) * t221 + t215 * (-t135 * t164 + t254 * t166) + (t415 * t411 - ((t195 + t256 + t421) * t215 + t398) * t221 + t166 * t395 + t415) * t98;
t4 = -t106 * t16 - t125 * t13 + t180 * t10 + 0.2e1 * t61 * t377 + (-t134 * t166 - t255 * t163) * t219 + t213 * (-t134 * t163 + t255 * t166) + (t416 * t412 - ((t193 + t256 + t420) * t213 + t399) * t219 + t166 * t396 + t416) * t97;
t3 = -t84 * t18 - t108 * t15 - t12 * t326 + 0.4e1 * (t57 - t375 / 0.2e1) * t364 + (((t197 + t200 + t247) * t63 + t36 * t180) * t363 + (t102 * t397 + t121 * t63 + t179 * t36) * t63) * t407 + ((-t169 * pkin(5) + t226 + t259) * t60 + (m(3) * t36 - t166 * t63) * t102 * t424) * t217 - 0.2e1 * t102 * (t57 - t375) + (-t146 + (-t153 - t250) * t139) * t224 + (t250 * g(3) - t139 * t150 + t151) * t218;
t2 = -t83 * t17 - t107 * t14 - t11 * t327 + 0.4e1 * (t56 - t376 / 0.2e1) * t367 + (((t195 + t200 + t248) * t62 + t35 * t180) * t366 + (t101 * t398 + t121 * t62 + t179 * t35) * t62) * t408 + ((-t168 * pkin(5) + t226 + t258) * t59 + (m(3) * t35 - t166 * t62) * t101 * t424) * t215 - 0.2e1 * t101 * (t56 - t376) + (-t146 + (-t153 - t251) * t138) * t222 + (t251 * g(3) - t138 * t150 + t151) * t216;
t1 = -t82 * t16 - t106 * t13 - t10 * t328 + 0.4e1 * (t55 - t377 / 0.2e1) * t370 + (((t193 + t200 + t249) * t61 + t34 * t180) * t369 + (t100 * t399 + t121 * t61 + t179 * t34) * t61) * t409 + ((-t167 * pkin(5) + t226 + t257) * t58 + (m(3) * t34 - t166 * t61) * t100 * t424) * t213 - 0.2e1 * t100 * (t55 - t377) + (-t146 + (-t153 - t252) * t137) * t220 + (t252 * g(3) - t137 * t150 + t151) * t214;
t37 = [(-g(1) + t212) * m(4) + (-t33 * t306 + (-t33 * t343 + (t212 * t33 + t3) * t178) * t224 + ((t45 * t93 + t54 * t81) * t212 + (t45 * t92 + t54 * t80) * t211 + (t105 * t54 + t45 * t359) * t210 + t93 * t6 + t81 * t9) * t237) * t133 + (-t32 * t311 + (-t32 * t347 + (t212 * t32 + t2) * t177) * t222 + ((t44 * t91 + t53 * t79) * t212 + (t44 * t90 + t53 * t78) * t211 + (t104 * t53 + t44 * t360) * t210 + t91 * t5 + t79 * t8) * t235) * t132 + (-t31 * t316 + (-t31 * t351 + (t212 * t31 + t1) * t176) * t220 + ((t43 * t89 + t52 * t77) * t212 + (t43 * t88 + t52 * t76) * t211 + (t103 * t52 + t43 * t361) * t210 + t89 * t4 + t77 * t7) * t233) * t131; (-g(2) + t211) * m(4) + (-t30 * t306 + (t30 * t331 + (-t211 * t30 - t3) * t175) * t224 + ((t42 * t93 + t51 * t81) * t212 + (t42 * t92 + t51 * t80) * t211 + (t105 * t51 + t42 * t359) * t210 + t92 * t6 + t80 * t9) * t237) * t133 + (-t29 * t311 + (t29 * t335 + (-t211 * t29 - t2) * t174) * t222 + ((t41 * t91 + t50 * t79) * t212 + (t41 * t90 + t50 * t78) * t211 + (t104 * t50 + t41 * t360) * t210 + t90 * t5 + t78 * t8) * t235) * t132 + (-t28 * t316 + (t28 * t339 + (-t211 * t28 - t1) * t173) * t220 + ((t40 * t89 + t49 * t77) * t212 + (t40 * t88 + t49 * t76) * t211 + (t103 * t49 + t40 * t361) * t210 + t88 * t4 + t76 * t7) * t233) * t131; (-g(3) + t210) * m(4) + (-t218 * t3 + (-t306 + (t331 - t343) * t224) * (-t218 * t84 + (t105 * t326 + t108 * t359) * t237) * t133 + ((t66 * t93 + t69 * t81) * t212 + (t66 * t92 + t69 * t80) * t211 + (t105 * t69 + t66 * t359) * t210 + t6 * t359 + t105 * t9) * t237) * t133 + (-t216 * t2 + (-t311 + (t335 - t347) * t222) * (-t216 * t83 + (t104 * t327 + t107 * t360) * t235) * t132 + ((t65 * t91 + t68 * t79) * t212 + (t65 * t90 + t68 * t78) * t211 + (t104 * t68 + t65 * t360) * t210 + t5 * t360 + t104 * t8) * t235) * t132 + (-t214 * t1 + (-t316 + (t339 - t351) * t220) * (-t214 * t82 + (t103 * t328 + t106 * t361) * t233) * t131 + ((t64 * t89 + t67 * t77) * t212 + (t64 * t88 + t67 * t76) * t211 + (t103 * t67 + t64 * t361) * t210 + t4 * t361 + t103 * t7) * t233) * t131;];
tauX  = t37;
