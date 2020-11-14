% Calculate vector of inverse dynamics forces for parallel robot
% P3RRPRR12V2G3A0
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
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:27:12
% EndTime: 2020-08-06 19:27:23
% DurationCPUTime: 11.63s
% Computational Cost: add. (86121->655), mult. (101859->1052), div. (11898->6), fcn. (76671->18), ass. (0->365)
t213 = sin(qJ(2,3));
t399 = pkin(1) * t213;
t157 = qJ(3,3) + t399;
t170 = t213 * qJ(3,3);
t214 = sin(qJ(1,3));
t219 = cos(qJ(2,3));
t201 = t219 ^ 2;
t231 = pkin(2) + pkin(3);
t325 = (qJ(3,3) + t231) * (-qJ(3,3) + t231);
t271 = t201 * t325;
t303 = t231 * t219;
t220 = cos(qJ(1,3));
t230 = pkin(5) - pkin(6);
t307 = t230 * t220;
t103 = -t214 * t271 - ((0.2e1 * t170 + pkin(1)) * t214 - t307) * t303 - qJ(3,3) * (t157 * t214 - t213 * t307);
t227 = xDP(3);
t228 = xDP(2);
t229 = xDP(1);
t154 = t170 + pkin(1);
t422 = t154 + t303;
t131 = 0.1e1 / t422;
t233 = 0.1e1 / qJ(3,3);
t364 = t131 * t233;
t300 = pkin(1) * t220 + t214 * t230;
t112 = qJ(3,3) * t220 + t300 * t213;
t322 = t213 * t220;
t278 = qJ(3,3) * t322;
t115 = 0.2e1 * t278 + t300;
t270 = t213 * t325;
t403 = pkin(1) * qJ(3,3);
t128 = -t270 + t403;
t207 = legFrame(3,2);
t173 = sin(t207);
t176 = cos(t207);
t269 = t220 * t325;
t409 = -0.2e1 * t231;
t289 = qJ(3,3) * t409;
t343 = t176 * t231;
t355 = t173 * t231;
t358 = t173 * qJ(3,3);
t76 = (-t173 * t269 + t176 * t289) * t201 + (-t115 * t355 - t128 * t176) * t219 - t112 * t358 + t157 * t343;
t346 = t176 * qJ(3,3);
t77 = (t173 * t289 + t176 * t269) * t201 + (t115 * t343 - t128 * t173) * t219 + t112 * t346 + t157 * t355;
t34 = (t103 * t227 + t228 * t76 + t229 * t77) * t364;
t109 = t422 * t214 - t307;
t367 = t109 * t219;
t116 = t278 + t300;
t312 = t220 * t231;
t320 = t213 * t231;
t88 = (-t173 * t312 - t346) * t201 + (-t116 * t173 + t176 * t320) * t219 + t176 * t157;
t89 = (t176 * t312 - t358) * t201 + (t116 * t176 + t173 * t320) * t219 + t173 * t157;
t61 = (-t227 * t367 + t228 * t88 + t229 * t89) * t364;
t388 = t231 * t61;
t26 = t34 - t388;
t215 = sin(qJ(2,2));
t398 = pkin(1) * t215;
t158 = qJ(3,2) + t398;
t171 = t215 * qJ(3,2);
t216 = sin(qJ(1,2));
t221 = cos(qJ(2,2));
t202 = t221 ^ 2;
t324 = (qJ(3,2) + t231) * (-qJ(3,2) + t231);
t268 = t202 * t324;
t302 = t231 * t221;
t222 = cos(qJ(1,2));
t306 = t230 * t222;
t104 = -t216 * t268 - ((0.2e1 * t171 + pkin(1)) * t216 - t306) * t302 - qJ(3,2) * (t158 * t216 - t215 * t306);
t155 = t171 + pkin(1);
t421 = t155 + t302;
t132 = 0.1e1 / t421;
t235 = 0.1e1 / qJ(3,2);
t363 = t132 * t235;
t299 = pkin(1) * t222 + t216 * t230;
t113 = qJ(3,2) * t222 + t299 * t215;
t319 = t215 * t222;
t279 = qJ(3,2) * t319;
t117 = 0.2e1 * t279 + t299;
t267 = t215 * t324;
t404 = pkin(1) * qJ(3,2);
t129 = -t267 + t404;
t208 = legFrame(2,2);
t174 = sin(t208);
t177 = cos(t208);
t266 = t222 * t324;
t290 = qJ(3,2) * t409;
t339 = t177 * t231;
t351 = t174 * t231;
t354 = t174 * qJ(3,2);
t78 = (-t174 * t266 + t177 * t290) * t202 + (-t117 * t351 - t129 * t177) * t221 - t113 * t354 + t158 * t339;
t342 = t177 * qJ(3,2);
t79 = (t174 * t290 + t177 * t266) * t202 + (t117 * t339 - t129 * t174) * t221 + t113 * t342 + t158 * t351;
t35 = (t104 * t227 + t228 * t78 + t229 * t79) * t363;
t110 = t421 * t216 - t306;
t366 = t110 * t221;
t118 = t279 + t299;
t310 = t222 * t231;
t317 = t215 * t231;
t90 = (-t174 * t310 - t342) * t202 + (-t118 * t174 + t177 * t317) * t221 + t177 * t158;
t91 = (t177 * t310 - t354) * t202 + (t118 * t177 + t174 * t317) * t221 + t174 * t158;
t62 = (-t227 * t366 + t228 * t90 + t229 * t91) * t363;
t387 = t231 * t62;
t27 = t35 - t387;
t217 = sin(qJ(2,1));
t397 = pkin(1) * t217;
t159 = qJ(3,1) + t397;
t172 = t217 * qJ(3,1);
t218 = sin(qJ(1,1));
t223 = cos(qJ(2,1));
t203 = t223 ^ 2;
t323 = (qJ(3,1) + t231) * (-qJ(3,1) + t231);
t265 = t203 * t323;
t301 = t231 * t223;
t224 = cos(qJ(1,1));
t305 = t230 * t224;
t105 = -t218 * t265 - ((0.2e1 * t172 + pkin(1)) * t218 - t305) * t301 - qJ(3,1) * (t159 * t218 - t217 * t305);
t156 = t172 + pkin(1);
t420 = t156 + t301;
t133 = 0.1e1 / t420;
t237 = 0.1e1 / qJ(3,1);
t362 = t133 * t237;
t298 = pkin(1) * t224 + t218 * t230;
t114 = qJ(3,1) * t224 + t298 * t217;
t316 = t217 * t224;
t280 = qJ(3,1) * t316;
t119 = 0.2e1 * t280 + t298;
t264 = t217 * t323;
t405 = pkin(1) * qJ(3,1);
t130 = -t264 + t405;
t209 = legFrame(1,2);
t175 = sin(t209);
t178 = cos(t209);
t263 = t224 * t323;
t291 = qJ(3,1) * t409;
t335 = t178 * t231;
t347 = t175 * t231;
t350 = t175 * qJ(3,1);
t80 = (-t175 * t263 + t178 * t291) * t203 + (-t119 * t347 - t130 * t178) * t223 - t114 * t350 + t159 * t335;
t338 = t178 * qJ(3,1);
t81 = (t175 * t291 + t178 * t263) * t203 + (t119 * t335 - t130 * t175) * t223 + t114 * t338 + t159 * t347;
t36 = (t105 * t227 + t228 * t80 + t229 * t81) * t362;
t111 = t420 * t218 - t305;
t365 = t111 * t223;
t120 = t280 + t298;
t308 = t224 * t231;
t314 = t217 * t231;
t92 = (-t175 * t308 - t338) * t203 + (-t120 * t175 + t178 * t314) * t223 + t178 * t159;
t93 = (t178 * t308 - t350) * t203 + (t120 * t178 + t175 * t314) * t223 + t175 * t159;
t63 = (-t227 * t365 + t228 * t92 + t229 * t93) * t362;
t386 = t231 * t63;
t25 = t36 - t386;
t427 = 0.2e1 * pkin(1);
t426 = 0.2e1 * t231;
t236 = qJ(3,1) ^ 2;
t240 = pkin(2) ^ 2;
t425 = (t236 - t240) * m(3);
t234 = qJ(3,2) ^ 2;
t424 = (t234 - t240) * m(3);
t232 = qJ(3,3) ^ 2;
t423 = (t232 - t240) * m(3);
t180 = m(3) * pkin(2) + mrSges(3,1);
t285 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t419 = t180 * qJ(3,3) + t285;
t418 = t180 * qJ(3,2) + t285;
t417 = t180 * qJ(3,1) + t285;
t415 = -0.2e1 * t201;
t414 = -0.2e1 * t202;
t413 = -0.2e1 * t203;
t412 = 0.2e1 * t219;
t411 = 0.2e1 * t221;
t410 = 0.2e1 * t223;
t408 = m(3) * pkin(5);
t407 = m(2) + m(3);
t406 = mrSges(3,1) * pkin(2);
t167 = m(3) * qJ(3,3) + mrSges(3,3);
t163 = -mrSges(2,2) + t167;
t402 = pkin(1) * t163;
t168 = m(3) * qJ(3,2) + mrSges(3,3);
t164 = -mrSges(2,2) + t168;
t401 = pkin(1) * t164;
t169 = qJ(3,1) * m(3) + mrSges(3,3);
t165 = -mrSges(2,2) + t169;
t400 = pkin(1) * t165;
t206 = mrSges(3,2) + mrSges(2,3);
t396 = Ifges(2,1) + Ifges(3,1);
t395 = Ifges(3,6) - Ifges(2,6);
t23 = t230 * t26;
t100 = (-t220 * t227 + (t173 * t228 - t176 * t229) * t214) * t131;
t286 = t100 * t403;
t394 = 0.2e1 * t286 + t23;
t24 = t230 * t27;
t101 = (-t222 * t227 + (t174 * t228 - t177 * t229) * t216) * t132;
t287 = t101 * t404;
t393 = 0.2e1 * t287 + t24;
t22 = t230 * t25;
t102 = (-t224 * t227 + (t175 * t228 - t178 * t229) * t218) * t133;
t288 = t102 * t405;
t392 = 0.2e1 * t288 + t22;
t391 = mrSges(3,3) * qJ(3,1);
t390 = mrSges(3,3) * qJ(3,2);
t389 = mrSges(3,3) * qJ(3,3);
t385 = t232 * t61;
t384 = t234 * t62;
t383 = t236 * t63;
t382 = t34 * t167;
t381 = t35 * t168;
t380 = t36 * t169;
t321 = t213 * t230;
t277 = t100 * t321;
t379 = (t277 - t388) * t219;
t318 = t215 * t230;
t276 = t101 * t318;
t378 = (t276 - t387) * t221;
t315 = t217 * t230;
t275 = t102 * t315;
t377 = (t275 - t386) * t223;
t376 = t100 * t201;
t375 = t100 * t213;
t374 = t100 * t230;
t373 = t101 * t202;
t372 = t101 * t215;
t371 = t101 * t230;
t370 = t102 * t203;
t369 = t102 * t217;
t368 = t102 * t230;
t361 = t163 * t213;
t360 = t164 * t215;
t359 = t165 * t217;
t211 = xDDP(2);
t357 = t173 * t211;
t356 = t173 * t214;
t353 = t174 * t211;
t352 = t174 * t216;
t349 = t175 * t211;
t348 = t175 * t218;
t212 = xDDP(1);
t345 = t176 * t212;
t344 = t176 * t214;
t341 = t177 * t212;
t340 = t177 * t216;
t337 = t178 * t212;
t336 = t178 * t218;
t179 = mrSges(3,2) + t408;
t334 = t179 * t213;
t333 = t179 * t215;
t332 = t179 * t217;
t328 = t180 * t219;
t327 = t180 * t221;
t326 = t180 * t223;
t210 = xDDP(3);
t313 = t220 * t210;
t311 = t222 * t210;
t309 = t224 * t210;
t304 = t230 * t231;
t294 = pkin(1) ^ 2 + pkin(5) ^ 2;
t199 = 0.2e1 * t406;
t284 = Ifges(3,2) + Ifges(2,3) + t199;
t283 = t219 * qJ(3,3) * t61;
t282 = t221 * qJ(3,2) * t62;
t281 = t63 * t223 * qJ(3,1);
t274 = t214 * t334;
t273 = t216 * t333;
t272 = t218 * t332;
t262 = -Ifges(2,2) - Ifges(3,3) + t396;
t260 = t294 + (-0.2e1 * pkin(5) + pkin(6)) * pkin(6);
t259 = -qJ(3,1) * mrSges(3,2) + t395;
t258 = -qJ(3,2) * mrSges(3,2) + t395;
t257 = -qJ(3,3) * mrSges(3,2) + t395;
t200 = -0.2e1 * t406;
t256 = t200 + t262;
t137 = t176 * g(1) - t173 * g(2);
t255 = g(3) * t214 - t137 * t220;
t138 = t177 * g(1) - t174 * g(2);
t254 = g(3) * t216 - t138 * t222;
t139 = t178 * g(1) - t175 * g(2);
t253 = g(3) * t218 - t139 * t224;
t19 = t231 * t26 - t385;
t20 = t231 * t27 - t384;
t21 = t231 * t25 - t383;
t166 = mrSges(2,1) + t180;
t252 = t166 * t219 + t361;
t251 = t166 * t221 + t360;
t250 = t166 * t223 + t359;
t249 = t262 + t423;
t248 = t262 + t424;
t247 = t262 + t425;
t246 = pkin(2) * mrSges(3,2) + t166 * pkin(5) - Ifges(3,4) - Ifges(2,5);
t245 = t294 * m(2) + 0.2e1 * t206 * pkin(5) + Ifges(1,3) + t396;
t226 = mrSges(2,2) * pkin(5);
t205 = t231 ^ 2;
t197 = 0.2e1 * t391;
t195 = 0.2e1 * t390;
t193 = 0.2e1 * t389;
t153 = t407 * pkin(1) + mrSges(1,1);
t152 = t166 * t427;
t151 = g(3) * t153;
t150 = t407 * pkin(5) - mrSges(1,2) + t206;
t149 = t236 + t260;
t148 = t234 + t260;
t147 = t232 + t260;
t146 = g(3) * t150;
t136 = t175 * g(1) + t178 * g(2);
t135 = t174 * g(1) + t177 * g(2);
t134 = t173 * g(1) + t176 * g(2);
t127 = m(3) * (t236 + t240) + t197 + t284;
t126 = m(3) * (t234 + t240) + t195 + t284;
t125 = m(3) * (t232 + t240) + t193 + t284;
t121 = (-mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1) * pkin(5) + Ifges(3,4) / 0.2e1 + Ifges(2,5) / 0.2e1 + (-t408 / 0.2e1 - mrSges(3,2) / 0.2e1) * pkin(2);
t108 = t223 * (t165 * pkin(5) - t259) - t217 * t246;
t107 = t221 * (t164 * pkin(5) - t258) - t215 * t246;
t106 = t219 * (t163 * pkin(5) - t257) - t213 * t246;
t99 = t102 ^ 2;
t98 = t101 ^ 2;
t97 = t100 ^ 2;
t84 = (t199 - t247 - 0.2e1 * t391) * t203 + t152 * t223 + t359 * t427 + (t236 + t294) * m(3) + t197 + t417 * t217 * t410 + t245;
t83 = (t199 - t248 - 0.2e1 * t390) * t202 + t152 * t221 + t360 * t427 + (t234 + t294) * m(3) + t195 + t418 * t215 * t411 + t245;
t82 = (t199 - t249 - 0.2e1 * t389) * t201 + t152 * t219 + t361 * t427 + (t232 + t294) * m(3) + t193 + t419 * t213 * t412 + t245;
t69 = (-t179 * t316 + (m(3) * t105 + t111 * t326) * t237) * t133;
t68 = (-t179 * t319 + (m(3) * t104 + t110 * t327) * t235) * t132;
t67 = (-t179 * t322 + (m(3) * t103 + t109 * t328) * t233) * t131;
t66 = (-t108 * t224 + (-t105 * t180 - t127 * t365) * t237) * t133;
t65 = (-t107 * t222 + (-t104 * t180 - t126 * t366) * t235) * t132;
t64 = (-t106 * t220 + (-t103 * t180 - t125 * t367) * t233) * t131;
t60 = t63 ^ 2;
t59 = t62 ^ 2;
t58 = t61 ^ 2;
t57 = t417 * t63;
t56 = t418 * t62;
t55 = t419 * t61;
t54 = (-t178 * t272 + (m(3) * t81 - t180 * t93) * t237) * t133;
t53 = (-t177 * t273 + (m(3) * t79 - t180 * t91) * t235) * t132;
t52 = (-t176 * t274 + (m(3) * t77 - t180 * t89) * t233) * t131;
t51 = (t175 * t272 + (m(3) * t80 - t180 * t92) * t237) * t133;
t50 = (t174 * t273 + (m(3) * t78 - t180 * t90) * t235) * t132;
t49 = (t173 * t274 + (m(3) * t76 - t180 * t88) * t233) * t131;
t45 = (-t108 * t336 + (t127 * t93 - t180 * t81) * t237) * t133;
t44 = (-t107 * t340 + (t126 * t91 - t180 * t79) * t235) * t132;
t43 = (-t106 * t344 + (t125 * t89 - t180 * t77) * t233) * t131;
t42 = (t108 * t348 + (t127 * t92 - t180 * t80) * t237) * t133;
t41 = (t107 * t352 + (t126 * t90 - t180 * t78) * t235) * t132;
t40 = (t106 * t356 + (t125 * t88 - t180 * t76) * t233) * t131;
t33 = (-t84 * t336 + (t108 * t93 + t81 * t332) * t237) * t133;
t32 = (-t83 * t340 + (t107 * t91 + t79 * t333) * t235) * t132;
t31 = (-t82 * t344 + (t106 * t89 + t77 * t334) * t233) * t131;
t30 = (t84 * t348 + (t108 * t92 + t80 * t332) * t237) * t133;
t29 = (t83 * t352 + (t107 * t90 + t78 * t333) * t235) * t132;
t28 = (t82 * t356 + (t106 * t88 + t76 * t334) * t233) * t131;
t18 = (0.2e1 * t25 * t217 + 0.2e1 * t281 + t368) * t102 * t133;
t17 = (0.2e1 * t27 * t215 + 0.2e1 * t282 + t371) * t101 * t132;
t16 = (0.2e1 * t26 * t213 + 0.2e1 * t283 + t374) * t100 * t131;
t15 = (-(t230 * t281 + t392 * t217 + (t156 * t223 * t426 + t149 + t265) * t102) * t223 * t102 + (-qJ(3,1) * t203 * t368 + (t25 + t275) * t301 + t25 * t156) * t63 + (t63 * t156 - t377) * t36) * t362;
t14 = (-(t230 * t282 + t393 * t215 + (t155 * t221 * t426 + t148 + t268) * t101) * t221 * t101 + (-qJ(3,2) * t202 * t371 + (t27 + t276) * t302 + t27 * t155) * t62 + (t62 * t155 - t378) * t35) * t363;
t13 = (-(t230 * t283 + t394 * t213 + (t154 * t219 * t426 + t147 + t271) * t100) * t219 * t100 + (-qJ(3,3) * t201 * t374 + (t26 + t277) * t303 + t26 * t154) * t61 + (t61 * t154 - t379) * t34) * t364;
t12 = ((-(t205 - 0.3e1 * t236) * t301 * t370 + (t230 * (-t426 * t63 + t36) * qJ(3,1) + (-0.3e1 * (-t236 / 0.3e1 + t205) * t172 + (t236 - t205) * t427) * t102) * t203 + (-t315 * t383 + ((-0.4e1 * t288 - t22) * t217 - t102 * (0.3e1 * t236 + t260)) * t231) * t223 - qJ(3,1) * (t149 * t369 + t392)) * t102 + ((t231 * t21 + t264 * t368) * t223 + t21 * pkin(1) + (t21 * t217 + (t413 + 0.1e1) * t102 * t304) * qJ(3,1)) * t63 + ((pkin(1) * t63 - t377) * t231 + (t63 * t314 + (t203 - 0.1e1) * t368) * qJ(3,1)) * t36) * t362;
t11 = ((-(t205 - 0.3e1 * t234) * t302 * t373 + (t230 * (-t426 * t62 + t35) * qJ(3,2) + (-0.3e1 * (-t234 / 0.3e1 + t205) * t171 + (t234 - t205) * t427) * t101) * t202 + (-t318 * t384 + ((-0.4e1 * t287 - t24) * t215 - t101 * (0.3e1 * t234 + t260)) * t231) * t221 - qJ(3,2) * (t148 * t372 + t393)) * t101 + ((t231 * t20 + t267 * t371) * t221 + t20 * pkin(1) + (t20 * t215 + (t414 + 0.1e1) * t101 * t304) * qJ(3,2)) * t62 + ((pkin(1) * t62 - t378) * t231 + (t62 * t317 + (t202 - 0.1e1) * t371) * qJ(3,2)) * t35) * t363;
t10 = ((-(t205 - 0.3e1 * t232) * t303 * t376 + (t230 * (-t426 * t61 + t34) * qJ(3,3) + (-0.3e1 * (-t232 / 0.3e1 + t205) * t170 + (t232 - t205) * t427) * t100) * t201 + (-t321 * t385 + ((-0.4e1 * t286 - t23) * t213 - t100 * (0.3e1 * t232 + t260)) * t231) * t219 - qJ(3,3) * (t147 * t375 + t394)) * t100 + ((t231 * t19 + t270 * t374) * t219 + t19 * pkin(1) + (t19 * t213 + (t415 + 0.1e1) * t100 * t304) * qJ(3,3)) * t61 + ((pkin(1) * t61 - t379) * t231 + (t61 * t320 + (t201 - 0.1e1) * t374) * qJ(3,3)) * t34) * t364;
t9 = t180 * t15 + (t99 * t203 - t60 - t99) * t169 + (t223 * t136 - t12) * m(3) + (-t99 * t326 - t179 * t18 + (-pkin(1) * t99 + t253) * m(3)) * t217;
t8 = t180 * t14 + (t98 * t202 - t59 - t98) * t168 + (t221 * t135 - t11) * m(3) + (-t98 * t327 - t179 * t17 + (-pkin(1) * t98 + t254) * m(3)) * t215;
t7 = t180 * t13 + (t97 * t201 - t58 - t97) * t167 + (t219 * t134 - t10) * m(3) + (-t97 * t328 - t179 * t16 + (-pkin(1) * t97 + t255) * m(3)) * t213;
t6 = -t108 * t18 - t127 * t15 + t180 * t12 + 0.2e1 * t63 * t380 + (-t136 * t166 + t253 * t165) * t223 + t217 * (-t136 * t165 - t253 * t166) + (t417 * t413 - ((t197 + t256 + t425) * t217 + t400) * t223 + t166 * t397 + t417) * t99;
t5 = -t107 * t17 - t126 * t14 + t180 * t11 + 0.2e1 * t62 * t381 + (-t135 * t166 + t254 * t164) * t221 + t215 * (-t135 * t164 - t254 * t166) + (t418 * t414 - ((t195 + t256 + t424) * t215 + t401) * t221 + t166 * t398 + t418) * t98;
t4 = -t106 * t16 - t125 * t13 + t180 * t10 + 0.2e1 * t61 * t382 + (-t134 * t166 + t255 * t163) * t219 + t213 * (-t134 * t163 - t255 * t166) + (t419 * t415 - ((t193 + t256 + t423) * t213 + t402) * t219 + t166 * t399 + t419) * t97;
t3 = -t84 * t18 - t108 * t15 - t12 * t332 + 0.4e1 * (t57 - t380 / 0.2e1) * t370 + (((t197 + t200 + t247) * t63 + t36 * t180) * t369 + t63 * (t102 * t400 + t121 * t63 + t36 * t179)) * t410 + ((-t169 * pkin(5) + t226 + t259) * t60 + (m(3) * t36 - t166 * t63) * t102 * t427) * t217 - 0.2e1 * t102 * (t57 - t380) + (g(3) * t250 - t150 * t139 + t151) * t224 + (t146 + (t153 + t250) * t139) * t218;
t2 = -t83 * t17 - t107 * t14 - t11 * t333 + 0.4e1 * (t56 - t381 / 0.2e1) * t373 + (((t195 + t200 + t248) * t62 + t35 * t180) * t372 + t62 * (t101 * t401 + t121 * t62 + t35 * t179)) * t411 + ((-t168 * pkin(5) + t226 + t258) * t59 + (m(3) * t35 - t166 * t62) * t101 * t427) * t215 - 0.2e1 * t101 * (t56 - t381) + (g(3) * t251 - t150 * t138 + t151) * t222 + (t146 + (t153 + t251) * t138) * t216;
t1 = -t82 * t16 - t106 * t13 - t10 * t334 + 0.4e1 * (t55 - t382 / 0.2e1) * t376 + (((t193 + t200 + t249) * t61 + t34 * t180) * t375 + t61 * (t100 * t402 + t121 * t61 + t34 * t179)) * t412 + ((-t167 * pkin(5) + t226 + t257) * t58 + (m(3) * t34 - t166 * t61) * t100 * t427) * t213 - 0.2e1 * t100 * (t55 - t382) + (g(3) * t252 - t150 * t137 + t151) * t220 + (t146 + (t153 + t252) * t137) * t214;
t37 = [(-g(1) + t212) * m(4) + (-t33 * t309 + (t33 * t349 + (-t212 * t33 - t3) * t178) * t218 + ((t45 * t93 + t54 * t81) * t212 + (t45 * t92 + t54 * t80) * t211 + (t105 * t54 - t45 * t365) * t210 + t93 * t6 + t81 * t9) * t237) * t133 + (-t32 * t311 + (t32 * t353 + (-t212 * t32 - t2) * t177) * t216 + ((t44 * t91 + t53 * t79) * t212 + (t44 * t90 + t53 * t78) * t211 + (t104 * t53 - t44 * t366) * t210 + t91 * t5 + t79 * t8) * t235) * t132 + (-t31 * t313 + (t31 * t357 + (-t212 * t31 - t1) * t176) * t214 + ((t43 * t89 + t52 * t77) * t212 + (t43 * t88 + t52 * t76) * t211 + (t103 * t52 - t43 * t367) * t210 + t89 * t4 + t77 * t7) * t233) * t131; (-g(2) + t211) * m(4) + (-t30 * t309 + (-t30 * t337 + (t211 * t30 + t3) * t175) * t218 + ((t42 * t93 + t51 * t81) * t212 + (t42 * t92 + t51 * t80) * t211 + (t105 * t51 - t42 * t365) * t210 + t92 * t6 + t80 * t9) * t237) * t133 + (-t29 * t311 + (-t29 * t341 + (t211 * t29 + t2) * t174) * t216 + ((t41 * t91 + t50 * t79) * t212 + (t41 * t90 + t50 * t78) * t211 + (t104 * t50 - t41 * t366) * t210 + t90 * t5 + t78 * t8) * t235) * t132 + (-t28 * t313 + (-t28 * t345 + (t211 * t28 + t1) * t173) * t214 + ((t40 * t89 + t49 * t77) * t212 + (t40 * t88 + t49 * t76) * t211 + (t103 * t49 - t40 * t367) * t210 + t88 * t4 + t76 * t7) * t233) * t131; (-g(3) + t210) * m(4) + (-t224 * t3 + (-t309 + (-t337 + t349) * t218) * (-t224 * t84 + (t105 * t332 - t108 * t365) * t237) * t133 + ((t66 * t93 + t69 * t81) * t212 + (t66 * t92 + t69 * t80) * t211 + (t105 * t69 - t66 * t365) * t210 - t6 * t365 + t105 * t9) * t237) * t133 + (-t222 * t2 + (-t311 + (-t341 + t353) * t216) * (-t222 * t83 + (t104 * t333 - t107 * t366) * t235) * t132 + ((t65 * t91 + t68 * t79) * t212 + (t65 * t90 + t68 * t78) * t211 + (t104 * t68 - t65 * t366) * t210 - t5 * t366 + t104 * t8) * t235) * t132 + (-t220 * t1 + (-t313 + (-t345 + t357) * t214) * (-t220 * t82 + (t103 * t334 - t106 * t367) * t233) * t131 + ((t64 * t89 + t67 * t77) * t212 + (t64 * t88 + t67 * t76) * t211 + (t103 * t67 - t64 * t367) * t210 - t4 * t367 + t103 * t7) * t233) * t131;];
tauX  = t37;
