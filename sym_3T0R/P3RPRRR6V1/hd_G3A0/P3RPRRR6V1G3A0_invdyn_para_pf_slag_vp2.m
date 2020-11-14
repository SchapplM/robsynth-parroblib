% Calculate vector of inverse dynamics forces for parallel robot
% P3RPRRR6V1G3A0
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-08-06 18:43
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX = P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G3A0_invdyn_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:41:33
% EndTime: 2020-08-06 18:41:44
% DurationCPUTime: 10.21s
% Computational Cost: add. (47637->636), mult. (53652->971), div. (8043->16), fcn. (36882->94), ass. (0->436)
t267 = xDDP(1);
t536 = t267 / 0.2e1;
t257 = cos(pkin(7));
t200 = t257 * pkin(1);
t162 = t200 + pkin(2);
t535 = mrSges(3,1) * t162;
t256 = sin(pkin(7));
t489 = pkin(1) * t256;
t350 = pkin(5) + t489;
t534 = t350 * mrSges(3,1);
t533 = t350 * mrSges(3,2);
t287 = -pkin(6) - pkin(5);
t419 = t256 * t287;
t138 = -pkin(1) + t419;
t275 = cos(qJ(1,3));
t119 = t138 * t275;
t269 = sin(qJ(1,3));
t414 = t269 * t287;
t422 = t256 * t269;
t532 = t119 + pkin(2) * t422 - (pkin(2) * t275 - t414) * t257;
t277 = cos(qJ(1,2));
t120 = t138 * t277;
t271 = sin(qJ(1,2));
t413 = t271 * t287;
t421 = t256 * t271;
t531 = t120 + pkin(2) * t421 - (pkin(2) * t277 - t413) * t257;
t279 = cos(qJ(1,1));
t121 = t138 * t279;
t273 = sin(qJ(1,1));
t412 = t273 * t287;
t420 = t256 * t273;
t530 = t121 + pkin(2) * t420 - (pkin(2) * t279 - t412) * t257;
t229 = qJ(1,1) + pkin(7);
t192 = cos(t229);
t265 = xDDP(3);
t264 = legFrame(1,2);
t154 = t264 + t229;
t155 = -t264 + t229;
t99 = -sin(t154) - sin(t155);
t529 = -t192 * t265 + t99 * t536;
t228 = qJ(1,2) + pkin(7);
t191 = cos(t228);
t263 = legFrame(2,2);
t152 = t263 + t228;
t153 = -t263 + t228;
t98 = -sin(t152) - sin(t153);
t528 = -t191 * t265 + t98 * t536;
t225 = qJ(1,3) + pkin(7);
t188 = cos(t225);
t262 = legFrame(3,2);
t150 = t262 + t225;
t151 = -t262 + t225;
t97 = -sin(t150) - sin(t151);
t527 = -t188 * t265 + t97 * t536;
t526 = 0.4e1 * t287;
t266 = xDDP(2);
t525 = t266 / 0.2e1;
t100 = cos(t151) - cos(t150);
t521 = t100 * t525 + t527;
t101 = cos(t153) - cos(t152);
t520 = t101 * t525 + t528;
t102 = cos(t155) - cos(t154);
t519 = t102 * t525 + t529;
t518 = 0.2e1 * pkin(1);
t306 = 0.2e1 * pkin(2);
t517 = 0.4e1 * pkin(2);
t459 = t97 / 0.2e1;
t458 = t98 / 0.2e1;
t457 = t99 / 0.2e1;
t448 = t100 / 0.2e1;
t447 = t101 / 0.2e1;
t446 = t102 / 0.2e1;
t285 = xDP(2);
t514 = t285 / 0.2e1;
t286 = xDP(1);
t513 = t286 / 0.2e1;
t207 = sin(t263);
t270 = sin(qJ(3,2));
t250 = 0.1e1 / t270;
t276 = cos(qJ(3,2));
t214 = t276 * pkin(3);
t169 = t214 + pkin(2);
t309 = 0.1e1 / pkin(3);
t438 = t169 * t277;
t453 = t309 * ((-t413 + t438) * t257 - t120 - t169 * t421);
t363 = t250 * t453;
t338 = t207 * t363;
t123 = 0.1e1 / (t200 + t169);
t441 = t123 * t285;
t70 = t338 * t441;
t284 = xDP(3);
t221 = pkin(7) + qJ(3,2);
t195 = qJ(1,2) + t221;
t158 = sin(t195);
t223 = -pkin(7) + qJ(3,2);
t196 = qJ(1,2) - t223;
t159 = sin(t196);
t176 = sin(t221);
t178 = sin(t223);
t182 = sin(t228);
t302 = 0.2e1 * qJ(3,2);
t235 = sin(t302);
t485 = pkin(3) * t235;
t486 = pkin(2) * t270;
t496 = 0.2e1 * t287;
t483 = (t191 * t496 + t271 * t518 + t182 * t306 + (t158 + t159) * pkin(3)) / (0.2e1 * t486 + t485 + (t176 + t178) * pkin(1));
t373 = t309 * t483;
t71 = t284 * t373;
t479 = t70 + t71;
t296 = 0.2e1 * pkin(7);
t218 = cos(t296);
t311 = pkin(1) ^ 2;
t174 = t311 * t218;
t308 = pkin(3) ^ 2;
t509 = -t174 - t308 / 0.2e1;
t217 = sin(t296);
t408 = t311 * t217;
t382 = 0.2e1 * t408;
t394 = -0.4e1 * t200;
t508 = t287 * t394 + t382;
t504 = -0.4e1 * pkin(3);
t206 = sin(t262);
t209 = cos(t262);
t274 = cos(qJ(3,3));
t213 = t274 * pkin(3);
t168 = t213 + pkin(2);
t122 = 0.1e1 / (t200 + t168);
t268 = sin(qJ(3,3));
t249 = 0.1e1 / t268;
t443 = t122 * t249;
t390 = pkin(7) + qJ(3,3);
t193 = qJ(1,3) + t390;
t156 = sin(t193);
t392 = -pkin(7) + qJ(3,3);
t194 = qJ(1,3) - t392;
t157 = sin(t194);
t179 = sin(t225);
t299 = 0.2e1 * qJ(3,3);
t232 = sin(t299);
t484 = (t188 * t496 + t269 * t518 + t179 * t306 + (t156 + t157) * pkin(3)) / (t268 * t306 + pkin(3) * t232 + (sin(t390) + sin(t392)) * pkin(1));
t439 = t168 * t275;
t73 = (-t414 + t439) * t257 - t119 - t168 * t422;
t32 = (t284 * t484 + (t206 * t285 - t209 * t286) * t73 * t443) * t309;
t503 = t32 ^ 2;
t210 = cos(t263);
t337 = t210 * t363;
t331 = t286 * t337;
t322 = t123 * t331;
t33 = -t322 + t479;
t502 = t33 ^ 2;
t208 = sin(t264);
t211 = cos(t264);
t278 = cos(qJ(3,1));
t215 = t278 * pkin(3);
t171 = t215 + pkin(2);
t124 = 0.1e1 / (t200 + t171);
t272 = sin(qJ(3,1));
t251 = 0.1e1 / t272;
t440 = t124 * t251;
t391 = pkin(7) + qJ(3,1);
t197 = qJ(1,1) + t391;
t160 = sin(t197);
t393 = -pkin(7) + qJ(3,1);
t198 = qJ(1,1) - t393;
t161 = sin(t198);
t183 = sin(t229);
t305 = 0.2e1 * qJ(3,1);
t238 = sin(t305);
t482 = (t192 * t496 + t273 * t518 + t183 * t306 + (t160 + t161) * pkin(3)) / (t272 * t306 + pkin(3) * t238 + (sin(t391) + sin(t393)) * pkin(1));
t437 = t171 * t279;
t77 = (-t412 + t437) * t257 - t121 - t171 * t420;
t34 = (t284 * t482 + (t208 * t285 - t211 * t286) * t77 * t440) * t309;
t501 = t34 ^ 2;
t500 = -2 * Ifges(3,4);
t388 = pkin(2) * t200;
t310 = pkin(2) ^ 2;
t291 = 0.2e1 * t310;
t404 = t291 + t311;
t103 = 0.4e1 * t388 + t174 + t404;
t499 = -0.2e1 * t308 - 0.2e1 * t103;
t173 = m(3) * pkin(5) - mrSges(2,2) + mrSges(3,3);
t498 = -0.2e1 * t173;
t212 = m(3) * pkin(2) + mrSges(2,1);
t497 = 0.2e1 * t212;
t495 = -0.6e1 * t308;
t295 = -0.2e1 * t311;
t494 = pkin(1) * pkin(3);
t289 = m(2) + m(3);
t493 = mrSges(3,2) * pkin(2);
t491 = g(3) * mrSges(1,2);
t490 = t309 / 0.2e1;
t147 = pkin(6) + t350;
t488 = pkin(2) * t147;
t301 = 0.3e1 * qJ(3,2);
t487 = pkin(2) * sin(t301);
t255 = t287 ^ 2;
t433 = t191 * t284;
t353 = t123 * t433;
t84 = t123 * t286 * t458;
t85 = t441 * t447;
t478 = t84 + t85;
t65 = -t353 + t478;
t63 = t311 * t65;
t481 = t63 + 0.3e1 * (t310 + t255 / 0.3e1) * t65;
t480 = 0.2e1 * t479;
t477 = mrSges(3,1) * t268;
t476 = mrSges(3,1) * t270;
t475 = mrSges(3,1) * t272;
t116 = g(1) * t209 - g(2) * t206;
t474 = mrSges(3,2) * t116;
t117 = g(1) * t210 - g(2) * t207;
t473 = mrSges(3,2) * t117;
t118 = g(1) * t211 - g(2) * t208;
t472 = mrSges(3,2) * t118;
t471 = mrSges(3,2) * t272;
t241 = cos(t299);
t341 = -0.2e1 * t287 * t200 + t408;
t93 = t341 + 0.2e1 * t488;
t470 = t241 * t93;
t247 = cos(t305);
t469 = t247 * t93;
t465 = t268 * t32;
t464 = t270 * t33;
t463 = t272 * t34;
t72 = (t168 * t269 + t275 * t287) * t257 - t138 * t269 + t256 * t439;
t462 = t274 * t72;
t74 = (t169 * t271 + t277 * t287) * t257 - t138 * t271 + t256 * t438;
t461 = t276 * t74;
t76 = (t171 * t273 + t279 * t287) * t257 - t138 * t273 + t256 * t437;
t460 = t278 * t76;
t455 = t287 * t65;
t454 = t309 * t73;
t452 = t309 * t77;
t252 = t274 ^ 2;
t399 = 0.2e1 * pkin(3);
t425 = t249 * t274;
t64 = (t100 * t514 - t188 * t284 + t97 * t513) * t122;
t378 = pkin(1) * t419;
t140 = -0.2e1 * t378;
t405 = t255 + t311;
t95 = t140 + t310 + 0.2e1 * t388 + t405;
t451 = -(pkin(3) * (-t162 - t213) * t503 * t249 + (-(t252 * t308 + t95) * t64 + (-t162 * t274 * t64 + t147 * t465) * t399) * t64 * t425) * t122 - g(1) * t206 - g(2) * t209;
t253 = t276 ^ 2;
t424 = t250 * t276;
t450 = -(pkin(3) * (-t162 - t214) * t502 * t250 + (-(t253 * t308 + t95) * t65 + (-t162 * t276 * t65 + t147 * t464) * t399) * t65 * t424) * t123 - g(1) * t207 - g(2) * t210;
t254 = t278 ^ 2;
t423 = t251 * t278;
t66 = (t102 * t514 - t192 * t284 + t99 * t513) * t124;
t449 = (pkin(3) * (-t162 - t215) * t501 * t251 + (-(t254 * t308 + t95) * t66 + (-t162 * t278 * t66 + t147 * t463) * t399) * t66 * t423) * t124 + g(1) * t208 + g(2) * t211;
t445 = t103 * t241;
t444 = t103 * t247;
t442 = t123 * t250;
t436 = t173 * t256;
t431 = t206 * t268;
t430 = t207 * t270;
t429 = t208 * t272;
t428 = t209 * t268;
t427 = t210 * t270;
t426 = t211 * t272;
t261 = Ifges(3,1) - Ifges(3,2);
t418 = t261 * t268;
t417 = t261 * t270;
t416 = t261 * t272;
t415 = t265 * t309;
t297 = 0.4e1 * qJ(3,3);
t411 = t308 * cos(t297);
t300 = 0.4e1 * qJ(3,2);
t410 = t308 * cos(t300);
t303 = 0.4e1 * qJ(3,1);
t409 = t308 * cos(t303);
t406 = t255 / 0.2e1 + t311;
t403 = -0.2e1 * t493;
t290 = 0.3e1 * t310;
t294 = 0.2e1 * t311;
t402 = (0.6e1 * t388 + t140 + t294 + t290 + t255 - t509) * t504;
t244 = cos(t302);
t397 = t244 * t517;
t396 = t162 * t504;
t395 = pkin(3) * t295;
t389 = t287 * t494;
t387 = pkin(3) * t465;
t386 = pkin(3) * t464;
t385 = pkin(3) * t463;
t384 = -0.2e1 * t147 * t308;
t383 = t162 * t495;
t381 = 0.3e1 / 0.4e1 * t84 + 0.3e1 / 0.4e1 * t85 - 0.3e1 / 0.4e1 * t353;
t167 = t213 + t306;
t380 = t167 * t200;
t170 = t215 + t306;
t379 = t170 * t200;
t377 = pkin(3) * (-t257 * t275 + t422) * t252;
t376 = pkin(3) * (-t257 * t277 + t421) * t253;
t375 = pkin(3) * (-t257 * t279 + t420) * t254;
t374 = t309 * t484;
t372 = t309 * t482;
t371 = t206 * t454;
t370 = t207 * t453;
t369 = t208 * t452;
t368 = t209 * t454;
t367 = t210 * t453;
t366 = t211 * t452;
t365 = t72 * t425;
t364 = t74 * t424;
t362 = t76 * t423;
t361 = t265 * t462;
t360 = t265 * t461;
t359 = t265 * t460;
t356 = -t70 / 0.2e1 - t71 / 0.2e1;
t355 = -t70 / 0.3e1 - t71 / 0.3e1;
t354 = -t70 / 0.6e1 - t71 / 0.6e1;
t220 = pkin(7) + t302;
t184 = cos(t220);
t222 = -pkin(7) + t302;
t186 = cos(t222);
t352 = t184 + t186 - 0.2e1 * t257;
t351 = -0.2e1 * t310 + t509;
t243 = cos(t301);
t348 = pkin(3) * (t243 - t276);
t347 = t147 * t399;
t345 = t252 * t500 + Ifges(3,4);
t344 = t253 * t500 + Ifges(3,4);
t343 = t254 * t500 + Ifges(3,4);
t342 = -t287 + t489;
t129 = Ifges(3,6) - t533;
t321 = Ifges(3,5) - t534;
t86 = t129 * t274 + t268 * t321;
t340 = t249 * t86 * t454;
t88 = t129 * t278 + t272 * t321;
t339 = t251 * t88 * t452;
t336 = -0.2e1 * t174 - 0.4e1 * t310 + t295 - t308;
t335 = 0.3e1 / 0.8e1 * t308 + t310 / 0.2e1 + t406;
t334 = -t311 + t351;
t332 = 0.2e1 * t348;
t327 = 0.2e1 * t352;
t185 = cos(t221);
t187 = cos(t223);
t226 = t301 + pkin(7);
t227 = t301 - pkin(7);
t326 = (-t185 - t187 + cos(t226) + cos(t227)) * pkin(3);
t325 = Ifges(3,6) / 0.2e1 - t533 / 0.2e1;
t324 = -Ifges(3,5) / 0.2e1 + t534 / 0.2e1;
t323 = mrSges(3,2) * t200 + t493;
t298 = 0.3e1 * qJ(3,3);
t240 = cos(t298);
t317 = pkin(3) * (-pkin(2) * t274 + t162 * t240);
t304 = 0.3e1 * qJ(3,1);
t246 = cos(t304);
t316 = pkin(3) * (-pkin(2) * t278 + t162 * t246);
t315 = -0.8e1 * (0.3e1 / 0.4e1 * t308 + t290 + t405) * t200 - 0.8e1 * (t335 - t378) * t306 + 0.8e1 * (-pkin(2) * t218 + t217 * t287) * t311;
t314 = (pkin(5) ^ 2 + t310) * m(3) + 0.2e1 * mrSges(3,3) * pkin(5) + t289 * t311 + Ifges(3,1) + Ifges(1,3) + Ifges(2,3);
t313 = t70 / 0.4e1 + t71 / 0.4e1 - t322 / 0.4e1;
t307 = pkin(3) * t308;
t282 = mrSges(3,1) * g(3);
t281 = mrSges(3,2) * g(3);
t237 = sin(t304);
t236 = sin(t303);
t233 = sin(t300);
t231 = sin(t298);
t230 = sin(t297);
t224 = -0.2e1 * pkin(7) + qJ(3,2);
t219 = t296 + qJ(3,2);
t204 = t270 * mrSges(3,2);
t203 = t268 * mrSges(3,2);
t181 = sin(t227);
t180 = sin(t226);
t177 = sin(t222);
t175 = sin(t220);
t165 = 0.2e1 * t223;
t164 = 0.2e1 * t221;
t149 = g(3) * t497;
t148 = t289 * pkin(1) + mrSges(1,1);
t146 = 0.2e1 * g(3) * t173;
t145 = cos(t165);
t144 = cos(t164);
t143 = sin(t165);
t142 = sin(t164);
t137 = mrSges(3,2) * t162;
t136 = g(3) * t148;
t135 = mrSges(3,2) * t278 + t475;
t134 = mrSges(3,2) * t276 + t476;
t133 = mrSges(3,2) * t274 + t477;
t132 = mrSges(3,1) * t276 - t204;
t131 = mrSges(3,1) * t274 - t203;
t130 = -t278 * mrSges(3,1) + t471;
t106 = mrSges(3,1) * t118;
t105 = mrSges(3,1) * t117;
t104 = mrSges(3,1) * t116;
t87 = t129 * t276 + t270 * t321;
t69 = -t261 * t254 + 0.2e1 * (Ifges(3,4) * t272 + t535) * t278 + t272 * t403 + (-(-t212 + t471) * t257 + t436) * t518 + t314;
t68 = -t261 * t253 + 0.2e1 * (Ifges(3,4) * t270 + t535) * t276 + t270 * t403 + (-(t204 - t212) * t257 + t436) * t518 + t314;
t67 = -t261 * t252 + 0.2e1 * (Ifges(3,4) * t268 + t535) * t274 + t268 * t403 + (-(t203 - t212) * t257 + t436) * t518 + t314;
t62 = -t211 * t375 + (pkin(3) * t429 - t530 * t211) * t278 + t162 * t429;
t61 = -t210 * t376 + (pkin(3) * t430 - t531 * t210) * t276 + t162 * t430;
t60 = -t209 * t377 + (pkin(3) * t431 - t532 * t209) * t274 + t162 * t431;
t59 = t208 * t375 + (pkin(3) * t426 + t530 * t208) * t278 + t162 * t426;
t58 = t207 * t376 + (pkin(3) * t427 + t531 * t207) * t276 + t162 * t427;
t57 = t206 * t377 + (pkin(3) * t428 + t532 * t206) * t274 + t162 * t428;
t55 = -t124 * t289 * t362 - t130 * t372;
t54 = -t123 * t289 * t364 + t132 * t373;
t53 = -t122 * t289 * t365 + t131 * t374;
t49 = (t130 * t366 + t289 * t62) * t440;
t48 = (-t132 * t367 + t289 * t61) * t442;
t47 = (-t131 * t368 + t289 * t60) * t443;
t46 = (-t130 * t369 + t289 * t59) * t440;
t45 = (t132 * t370 + t289 * t58) * t442;
t44 = (t131 * t371 + t289 * t57) * t443;
t43 = Ifges(3,3) * t372 + (t130 * t362 - t192 * t88) * t124;
t42 = Ifges(3,3) * t373 + (-t132 * t364 - t191 * t87) * t123;
t41 = Ifges(3,3) * t374 + (-t131 * t365 - t188 * t86) * t122;
t31 = (t88 * t446 + (Ifges(3,3) * t369 - t130 * t59) * t251) * t124;
t30 = (t87 * t447 + (Ifges(3,3) * t370 + t132 * t58) * t250) * t123;
t29 = (t86 * t448 + (Ifges(3,3) * t371 + t131 * t57) * t249) * t122;
t28 = (t88 * t457 + (-Ifges(3,3) * t366 - t130 * t62) * t251) * t124;
t27 = (t87 * t458 + (-Ifges(3,3) * t367 + t132 * t61) * t250) * t123;
t26 = (t86 * t459 + (-Ifges(3,3) * t368 + t131 * t60) * t249) * t122;
t25 = (-0.2e1 * t331 - t433) * t123 + t478 + t480;
t24 = (0.2e1 * t331 - t433) * t123 + t478 - t480;
t21 = t342 * t65 - t386;
t20 = t342 * t64 - t387;
t19 = t342 * t66 - t385;
t18 = (t21 - t386) * t65 * t123;
t17 = (t20 - t387) * t64 * t122;
t16 = (t19 - t385) * t66 * t124;
t12 = ((t246 * t384 + (t147 * t170 + t341 - t469) * t399) * t34 + (-t307 * t236 + t237 * t383 + t238 * t402 + t315 * t272) * t66) / (t444 + t409 / 0.2e1 - 0.2e1 * t379 + 0.2e1 * t316 + t334) * t66 * t490 + (t19 * t517 + t385 * t394 + (-0.2e1 * t469 + (-t246 + t278) * t347 + t508) * t66 + (-t308 * t236 + t237 * t396 + t238 * t499) * t34) / (0.4e1 * t316 + t336 - 0.4e1 * t379 + t409 + 0.2e1 * t444) * t34;
t11 = ((t240 * t384 + (t147 * t167 + t341 - t470) * t399) * t32 + (-t307 * t230 + t231 * t383 + t232 * t402 + t315 * t268) * t64) / (t445 + t411 / 0.2e1 - 0.2e1 * t380 + 0.2e1 * t317 + t334) * t64 * t490 + (t20 * t517 + t387 * t394 + (-0.2e1 * t470 + (-t240 + t274) * t347 + t508) * t64 + (-t308 * t230 + t231 * t396 + t232 * t499) * t32) / (0.4e1 * t317 + t336 - 0.4e1 * t380 + t411 + 0.2e1 * t445) * t32;
t10 = t12 * t130 - t135 * t501 - t449 * t289;
t9 = -t11 * t131 - t133 * t503 + t451 * t289;
t8 = -t88 * t16 - Ifges(3,3) * t12 + t66 ^ 2 * ((t137 - t416) * t278 + t162 * t475 + t343) + (-g(3) * t183 + t118 * t192) * t135 + t449 * t130;
t7 = -t86 * t17 - Ifges(3,3) * t11 + t64 ^ 2 * ((t137 - t418) * t274 + t162 * t477 + t345) + (-g(3) * t179 + t116 * t188) * t133 + t451 * t131;
t6 = -t69 * t16 - t88 * t12 - 0.2e1 * ((t325 * t272 + t324 * t278) * t34 + ((t323 - t416) * t278 + t272 * t535 + t343) * t66) * t34 + (t282 - t472) * cos(t198) / 0.2e1 + (t106 + t281) * t161 / 0.2e1 + (t282 + t472) * cos(t197) / 0.2e1 + (t106 - t281) * t160 / 0.2e1 + (t118 * t498 + t149) * t192 / 0.2e1 + (t118 * t497 + t146) * t183 / 0.2e1 + (mrSges(1,2) * t118 + t136) * t279 + t273 * (t118 * t148 - t491);
t5 = -t67 * t17 - t86 * t11 - 0.2e1 * ((t325 * t268 + t324 * t274) * t32 + ((t323 - t418) * t274 + t268 * t535 + t345) * t64) * t32 + (t282 - t474) * cos(t194) / 0.2e1 + (t104 + t281) * t157 / 0.2e1 + (t282 + t474) * cos(t193) / 0.2e1 + (t104 - t281) * t156 / 0.2e1 + (t116 * t498 + t149) * t188 / 0.2e1 + (t116 * t497 + t146) * t179 / 0.2e1 + (mrSges(1,2) * t116 + t136) * t275 + t269 * (t116 * t148 - t491);
t4 = (-0.4e1 * (((t313 + t381) * t308 + t481) * t178 + ((-t313 + t381) * t308 + t481) * t176) * pkin(1) - 0.3e1 * (((-t433 + t331 / 0.3e1) * t123 + t355 + t478) * t181 + ((-t433 - t331 / 0.3e1) * t123 - t355 + t478) * t180) * pkin(1) * t308 + (cos(t224) - cos(t219)) * t63 * t526 + (t142 * t395 + 0.4e1 * t186 * t389) * ((-t433 - t331 / 0.2e1) * t123 - t356 + t478) + (t143 * t395 - 0.4e1 * t184 * t389) * ((-t433 + t331 / 0.2e1) * t123 + t356 + t478) + (-0.8e1 * (t308 / 0.4e1 + 0.3e1 / 0.2e1 * t310 + t406) * t485 + t487 * t495 - t307 * t233 - 0.16e2 * t335 * t486) * t65 + ((t382 + 0.4e1 * t488 + (-t200 - t214 / 0.2e1) * t526) * pkin(3) + (pkin(3) * t397 + 0.2e1 * t308 * t243) * t287) * t33 + (-0.12e2 * (((-t433 + t331 / 0.6e1) * t123 + t354 + t478) * t177 + ((-t433 - t331 / 0.6e1) * t123 - t354 + t478) * t175) * t494 - 0.4e1 * (sin(t224) + sin(t219)) * t63 + 0.8e1 * (-t185 + t187) * pkin(1) * t455) * pkin(2)) / (t291 * t244 + t410 / 0.2e1 + (t145 / 0.2e1 + t144 / 0.2e1 + t244 - 0.1e1) * t311 + pkin(2) * t332 + (pkin(2) * t327 + t326) * pkin(1) + t351) * t65 * t490 + (t21 * t517 + ((t478 - t479) * t143 - (t478 + t479) * t142 + ((t331 - t433) * t143 - (-t331 - t433) * t142) * t123) * t311 + (-0.2e1 * (t308 + t404) * t235 + t487 * t504 - t308 * t233) * t33 + (t382 + (t332 + t397) * t287) * t65 + ((-t175 * t25 + t177 * t24) * t306 + t327 * t455 + ((-t178 - t180) * t25 + (t176 + t181) * t24) * pkin(3)) * pkin(1)) / ((t294 + 0.4e1 * t310) * t244 + t410 + (t145 + t144) * t311 + t348 * t517 + (t352 * t517 + 0.2e1 * t326) * pkin(1) + t336) * t33;
t3 = -t132 * t4 - t134 * t502 + t450 * t289;
t2 = -t87 * t18 - Ifges(3,3) * t4 + t65 ^ 2 * ((t137 - t417) * t276 + t162 * t476 + t344) + (-g(3) * t182 + t117 * t191) * t134 + t450 * t132;
t1 = -t68 * t18 - t87 * t4 - 0.2e1 * ((t325 * t270 + t324 * t276) * t33 + ((t323 - t417) * t276 + t270 * t535 + t344) * t65) * t33 + (t282 - t473) * cos(t196) / 0.2e1 + (t105 + t281) * t159 / 0.2e1 + (t282 + t473) * cos(t195) / 0.2e1 + (t105 - t281) * t158 / 0.2e1 + (t117 * t498 + t149) * t191 / 0.2e1 + (t117 * t497 + t146) * t182 / 0.2e1 + (mrSges(1,2) * t117 + t136) * t277 + t271 * (t117 * t148 - t491);
t13 = [(-g(1) + t267) * m(4) + (t26 * t484 + t27 * t483 + t28 * t482) * t415 + (t6 * t457 + ((-t28 * t366 + t49 * t62) * t267 + (t28 * t369 + t49 * t59) * t266 - t49 * t359 + t62 * t10 - t8 * t366) * t251 + t519 * (-t211 * t339 + t69 * t457) * t124) * t124 + (t1 * t458 + ((-t27 * t367 + t48 * t61) * t267 + (t27 * t370 + t48 * t58) * t266 - t48 * t360 + t61 * t3 - t2 * t367) * t250 + t520 * (-t87 * t337 + t68 * t458) * t123) * t123 + (t5 * t459 + ((-t26 * t368 + t47 * t60) * t267 + (t26 * t371 + t47 * t57) * t266 - t47 * t361 + t60 * t9 - t7 * t368) * t249 + t521 * (-t209 * t340 + t67 * t459) * t122) * t122; (-g(2) + t266) * m(4) + (t29 * t484 + t30 * t483 + t31 * t482) * t415 + (t6 * t446 + ((-t31 * t366 + t46 * t62) * t267 + (t31 * t369 + t46 * t59) * t266 - t46 * t359 + t59 * t10 + t8 * t369) * t251 + (t266 * t446 + t529) * (t208 * t339 + t69 * t446) * t124) * t124 + (t1 * t447 + ((-t30 * t367 + t45 * t61) * t267 + (t30 * t370 + t45 * t58) * t266 - t45 * t360 + t58 * t3 + t2 * t370) * t250 + (t266 * t447 + t528) * (t87 * t338 + t68 * t447) * t123) * t123 + (t5 * t448 + ((-t29 * t368 + t44 * t60) * t267 + (t29 * t371 + t44 * t57) * t266 - t44 * t361 + t57 * t9 + t7 * t371) * t249 + (t266 * t448 + t527) * (t206 * t340 + t67 * t448) * t122) * t122; (-g(3) + t265) * m(4) + (t2 * t483 + t7 * t484 + t8 * t482 + (t41 * t484 + t42 * t483 + t43 * t482) * t265) * t309 + (-t192 * t6 + t519 * (-t124 * t192 * t69 + t88 * t372) + ((-t43 * t366 + t55 * t62) * t267 + (t43 * t369 + t55 * t59) * t266 + (-t265 * t55 - t10) * t460) * t251) * t124 + (-t191 * t1 + t520 * (-t123 * t191 * t68 + t87 * t373) + ((-t42 * t367 + t54 * t61) * t267 + (t42 * t370 + t54 * t58) * t266 + (-t265 * t54 - t3) * t461) * t250) * t123 + (-t188 * t5 + t521 * (-t122 * t188 * t67 + t86 * t374) + ((-t41 * t368 + t53 * t60) * t267 + (t41 * t371 + t53 * t57) * t266 + (-t265 * t53 - t9) * t462) * t249) * t122;];
tauX  = t13;
