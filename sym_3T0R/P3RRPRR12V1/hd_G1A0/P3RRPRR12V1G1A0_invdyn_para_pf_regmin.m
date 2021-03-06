% Calculate minimal parameter regressor of inverse dynamics forces for
% P3RRPRR12V1G1A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tauX_reg [3x15]
%   minimal parameter regressor of inverse dynamics force vector
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:02
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauX_reg = P3RRPRR12V1G1A0_invdyn_para_pf_regmin(xP, xDP, xDDP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(4,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(xDDP) && all(size(xDDP) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_regmin: xDDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_regmin: pkin has to be [4x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G1A0_invdyn_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_reg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:01:54
% EndTime: 2020-08-06 19:02:07
% DurationCPUTime: 12.72s
% Computational Cost: add. (57231->597), mult. (103806->1072), div. (6153->12), fcn. (66936->18), ass. (0->479)
t287 = sin(qJ(2,3));
t288 = sin(qJ(1,3));
t261 = t288 * g(1);
t294 = cos(qJ(1,3));
t541 = t294 * g(2);
t243 = -t261 + t541;
t244 = g(1) * t294 + g(2) * t288;
t281 = legFrame(3,3);
t252 = sin(t281);
t255 = cos(t281);
t574 = t243 * t252 + t244 * t255;
t587 = t287 * t574;
t289 = sin(qJ(2,2));
t290 = sin(qJ(1,2));
t262 = t290 * g(1);
t296 = cos(qJ(1,2));
t540 = t296 * g(2);
t245 = -t262 + t540;
t246 = g(1) * t296 + g(2) * t290;
t282 = legFrame(2,3);
t253 = sin(t282);
t256 = cos(t282);
t575 = t245 * t253 + t246 * t256;
t586 = t289 * t575;
t291 = sin(qJ(2,1));
t292 = sin(qJ(1,1));
t266 = g(1) * t292;
t298 = cos(qJ(1,1));
t548 = g(2) * t298;
t247 = -t266 + t548;
t248 = g(1) * t298 + g(2) * t292;
t283 = legFrame(1,3);
t254 = sin(t283);
t257 = cos(t283);
t576 = t247 * t254 + t248 * t257;
t585 = t291 * t576;
t180 = -t252 * t288 + t255 * t294;
t263 = t294 * pkin(4);
t521 = qJ(3,3) * t287;
t213 = t288 * t521 + t263;
t552 = pkin(4) * t288;
t214 = t294 * t521 - t552;
t293 = cos(qJ(2,3));
t302 = pkin(1) + pkin(2);
t447 = t293 * t302;
t118 = t180 * t447 - t213 * t252 + t214 * t255;
t181 = t252 * t294 + t255 * t288;
t119 = t181 * t447 + t213 * t255 + t214 * t252;
t300 = xDP(2);
t301 = xDP(1);
t207 = -t288 * t301 + t294 * t300;
t208 = t288 * t300 + t294 * t301;
t160 = t207 * t255 - t208 * t252;
t303 = qJ(3,3) ^ 2;
t312 = pkin(4) ^ 2;
t258 = t303 + t312;
t277 = t293 ^ 2;
t284 = xDDP(3);
t285 = xDDP(2);
t286 = xDDP(1);
t304 = 0.1e1 / qJ(3,3);
t376 = (qJ(3,3) + t302) * (-qJ(3,3) + t302) * t277;
t237 = t447 + t521;
t225 = 0.1e1 / t237;
t455 = t287 * t302;
t299 = xDP(3);
t519 = qJ(3,3) * t299;
t142 = -t519 + (t207 * t252 + t208 * t255) * t302;
t549 = pkin(4) * t301;
t215 = t300 * t521 - t549;
t267 = pkin(4) * t300;
t216 = t301 * t521 + t267;
t440 = t299 * t302;
t372 = (t215 * t288 + t216 * t294) * t255 + (t215 * t294 - t216 * t288) * t252 + t287 * t440;
t106 = t142 * t277 + t372 * t293 + t519;
t514 = t106 * t304;
t563 = -pkin(4) / 0.2e1;
t567 = -0.2e1 * qJ(3,3);
t430 = (t160 * t455 + t514 * t563) * t225 * t567;
t479 = t225 * t293;
t226 = 0.1e1 / t237 ^ 2;
t493 = t160 * t226;
t494 = t160 * t225;
t112 = t142 * t293 + t372;
t454 = t287 * t304;
t515 = t106 * t225;
t67 = (-t302 * t515 + t112) * pkin(4) * t454;
t584 = ((t118 * t286 + t119 * t285) * t479 + t284 * t287 - (-t160 * t376 * t479 + t277 * t430 + (-t258 * t494 + t67) * t293) * t493) * t304;
t182 = -t253 * t290 + t256 * t296;
t264 = t296 * pkin(4);
t524 = qJ(3,2) * t289;
t217 = t290 * t524 + t264;
t551 = pkin(4) * t290;
t218 = t296 * t524 - t551;
t295 = cos(qJ(2,2));
t444 = t295 * t302;
t120 = t182 * t444 - t217 * t253 + t218 * t256;
t183 = t253 * t296 + t256 * t290;
t121 = t183 * t444 + t217 * t256 + t218 * t253;
t209 = -t290 * t301 + t296 * t300;
t210 = t290 * t300 + t296 * t301;
t161 = t209 * t256 - t210 * t253;
t306 = qJ(3,2) ^ 2;
t259 = t306 + t312;
t278 = t295 ^ 2;
t307 = 0.1e1 / qJ(3,2);
t375 = (qJ(3,2) + t302) * (-qJ(3,2) + t302) * t278;
t238 = t444 + t524;
t228 = 0.1e1 / t238;
t452 = t289 * t302;
t522 = qJ(3,2) * t299;
t143 = -t522 + (t209 * t253 + t210 * t256) * t302;
t219 = t300 * t524 - t549;
t220 = t301 * t524 + t267;
t371 = (t219 * t290 + t220 * t296) * t256 + (t219 * t296 - t220 * t290) * t253 + t289 * t440;
t107 = t143 * t278 + t371 * t295 + t522;
t512 = t107 * t307;
t568 = -0.2e1 * qJ(3,2);
t431 = (t161 * t452 + t512 * t563) * t228 * t568;
t445 = t295 * t228;
t229 = 0.1e1 / t238 ^ 2;
t491 = t161 * t229;
t492 = t161 * t228;
t113 = t143 * t295 + t371;
t451 = t289 * t307;
t513 = t107 * t228;
t68 = (-t302 * t513 + t113) * pkin(4) * t451;
t583 = ((t120 * t286 + t121 * t285) * t445 + t284 * t289 - (-t161 * t375 * t445 + t278 * t431 + (-t259 * t492 + t68) * t295) * t491) * t307;
t184 = -t254 * t292 + t257 * t298;
t265 = t298 * pkin(4);
t527 = qJ(3,1) * t291;
t223 = t292 * t527 + t265;
t550 = pkin(4) * t292;
t224 = t298 * t527 - t550;
t297 = cos(qJ(2,1));
t442 = t297 * t302;
t122 = t184 * t442 - t223 * t254 + t224 * t257;
t185 = t254 * t298 + t257 * t292;
t123 = t185 * t442 + t223 * t257 + t224 * t254;
t211 = -t292 * t301 + t298 * t300;
t212 = t292 * t300 + t298 * t301;
t162 = t211 * t257 - t212 * t254;
t309 = qJ(3,1) ^ 2;
t260 = t309 + t312;
t279 = t297 ^ 2;
t310 = 0.1e1 / qJ(3,1);
t374 = (qJ(3,1) + t302) * (-qJ(3,1) + t302) * t279;
t239 = t442 + t527;
t231 = 0.1e1 / t239;
t449 = t291 * t302;
t525 = qJ(3,1) * t299;
t144 = -t525 + (t211 * t254 + t212 * t257) * t302;
t221 = t300 * t527 - t549;
t222 = t301 * t527 + t267;
t370 = (t221 * t292 + t222 * t298) * t257 + (t221 * t298 - t222 * t292) * t254 + t291 * t440;
t108 = t144 * t279 + t370 * t297 + t525;
t510 = t108 * t310;
t569 = -0.2e1 * qJ(3,1);
t432 = (t162 * t449 + t510 * t563) * t231 * t569;
t472 = t231 * t297;
t232 = 0.1e1 / t239 ^ 2;
t489 = t162 * t232;
t490 = t162 * t231;
t114 = t144 * t297 + t370;
t448 = t291 * t310;
t511 = t108 * t231;
t69 = (-t302 * t511 + t114) * pkin(4) * t448;
t582 = ((t122 * t286 + t123 * t285) * t472 + t284 * t291 - (-t162 * t374 * t472 + t279 * t432 + (-t260 * t490 + t69) * t297) * t489) * t310;
t547 = g(3) * t287;
t146 = t293 * t574 + t547;
t546 = g(3) * t289;
t148 = t295 * t575 + t546;
t545 = g(3) * t291;
t150 = t297 * t576 + t545;
t190 = t244 * t252;
t581 = t243 * t255 - t190;
t194 = t246 * t253;
t580 = t245 * t256 - t194;
t198 = t248 * t254;
t436 = -t247 * t257 + t198;
t153 = t162 ^ 2;
t496 = t153 * t232;
t141 = pkin(1) * t496;
t311 = 0.1e1 / qJ(3,1) ^ 2;
t470 = t231 * t311;
t401 = t114 * t470;
t564 = 0.2e1 * t108;
t87 = t401 * t564;
t579 = -0.2e1 * t279 * t141 + t141 + t87;
t152 = t161 ^ 2;
t498 = t152 * t229;
t140 = pkin(1) * t498;
t308 = 0.1e1 / qJ(3,2) ^ 2;
t474 = t228 * t308;
t402 = t113 * t474;
t565 = 0.2e1 * t107;
t86 = t402 * t565;
t578 = -0.2e1 * t278 * t140 + t140 + t86;
t151 = t160 ^ 2;
t500 = t151 * t226;
t139 = pkin(1) * t500;
t305 = 0.1e1 / qJ(3,3) ^ 2;
t477 = t225 * t305;
t403 = t112 * t477;
t566 = 0.2e1 * t106;
t85 = t403 * t566;
t577 = -0.2e1 * t277 * t139 + t139 + t85;
t570 = 0.2e1 * pkin(1);
t562 = pkin(1) * g(3);
t429 = pkin(4) * t494;
t130 = t287 * t429;
t478 = t225 * t304;
t406 = t106 * t478;
t109 = t112 * t304;
t98 = pkin(1) * t406;
t80 = t98 - t109;
t426 = pkin(2) * t406 + t80;
t536 = t106 * (-(t130 + t426) * t447 + (t277 * t429 - t287 * t426) * qJ(3,3));
t94 = t302 * t406 + t130;
t73 = t287 * t515 + t293 * t94;
t34 = (-t112 * t225 * t73 - t226 * t536) * t305 + t584;
t561 = t34 * pkin(1);
t428 = pkin(4) * t492;
t131 = t289 * t428;
t475 = t228 * t307;
t405 = t107 * t475;
t100 = pkin(1) * t405;
t110 = t113 * t307;
t82 = t100 - t110;
t410 = pkin(2) * t405 + t82;
t535 = t107 * (-(t131 + t410) * t444 + (t278 * t428 - t289 * t410) * qJ(3,2));
t95 = t302 * t405 + t131;
t74 = t289 * t513 + t295 * t95;
t35 = (-t113 * t228 * t74 - t229 * t535) * t308 + t583;
t560 = t35 * pkin(1);
t427 = pkin(4) * t490;
t132 = t291 * t427;
t471 = t231 * t310;
t404 = t108 * t471;
t102 = pkin(1) * t404;
t111 = t114 * t310;
t84 = t102 - t111;
t373 = pkin(2) * t404 + t84;
t534 = t108 * (-(t132 + t373) * t442 + (t279 * t427 - t291 * t373) * qJ(3,1));
t96 = t302 * t404 + t132;
t75 = t291 * t511 + t297 * t96;
t36 = (-t114 * t231 * t75 - t232 * t534) * t311 + t582;
t559 = t36 * pkin(1);
t558 = 0.2e1 * t277 - 0.1e1;
t557 = 0.2e1 * t278 - 0.1e1;
t556 = 0.2e1 * t279 - 0.1e1;
t555 = pkin(4) * t160;
t554 = pkin(4) * t161;
t553 = pkin(4) * t162;
t544 = g(3) * t293;
t543 = g(3) * t295;
t542 = g(3) * t297;
t233 = t231 * t232;
t469 = t232 * t310;
t380 = t162 * t469;
t351 = t291 * t380;
t437 = t302 * t310;
t526 = qJ(3,1) * t297;
t344 = -t449 + t526;
t466 = t344 * t310;
t480 = t185 * t231;
t481 = t184 * t231;
t54 = -t286 * t480 + t285 * t481 - (t111 * t291 + (-t553 + (-t291 * t437 + t297) * t108) * t231) * t489 - t162 * t233 * t108 * t466 - t114 * t351;
t539 = qJ(3,1) * t54;
t230 = t228 * t229;
t473 = t229 * t307;
t382 = t161 * t473;
t352 = t289 * t382;
t438 = t302 * t307;
t523 = qJ(3,2) * t295;
t343 = -t452 + t523;
t467 = t343 * t307;
t482 = t183 * t228;
t483 = t182 * t228;
t53 = -t286 * t482 + t285 * t483 - (t110 * t289 + (-t554 + (-t289 * t438 + t295) * t107) * t228) * t491 - t161 * t230 * t107 * t467 - t113 * t352;
t538 = qJ(3,2) * t53;
t227 = t225 * t226;
t476 = t226 * t304;
t384 = t160 * t476;
t353 = t287 * t384;
t439 = t302 * t304;
t520 = qJ(3,3) * t293;
t342 = -t455 + t520;
t468 = t342 * t304;
t484 = t181 * t225;
t485 = t180 * t225;
t52 = -t286 * t484 + t285 * t485 - (t109 * t287 + (-t555 + (-t287 * t439 + t293) * t106) * t225) * t493 - t160 * t227 * t106 * t468 - t112 * t353;
t537 = qJ(3,3) * t52;
t533 = t34 * qJ(3,3);
t532 = t35 * qJ(3,2);
t531 = t36 * qJ(3,1);
t103 = t106 ^ 2;
t518 = t103 * t305;
t104 = t107 ^ 2;
t517 = t104 * t308;
t105 = t108 ^ 2;
t516 = t105 * t311;
t174 = t237 * t288 + t263;
t177 = t237 * t294 - t552;
t133 = -t174 * t252 + t177 * t255;
t509 = t133 * t304;
t175 = t238 * t290 + t264;
t178 = t238 * t296 - t551;
t134 = -t175 * t253 + t178 * t256;
t508 = t134 * t307;
t176 = t239 * t292 + t265;
t179 = t239 * t298 - t550;
t135 = -t176 * t254 + t179 * t257;
t507 = t135 * t310;
t136 = t174 * t255 + t177 * t252;
t506 = t136 * t304;
t137 = t175 * t256 + t178 * t253;
t505 = t137 * t307;
t138 = t176 * t257 + t179 * t254;
t504 = t138 * t310;
t145 = -t544 + t587;
t503 = t145 * t304;
t147 = -t543 + t586;
t502 = t147 * t307;
t149 = -t542 + t585;
t501 = t149 * t310;
t499 = t151 * t277;
t497 = t152 * t278;
t495 = t153 * t279;
t201 = g(1) * t252 - g(2) * t255;
t204 = g(1) * t255 + g(2) * t252;
t163 = t201 * t294 + t204 * t288;
t488 = t163 * t287;
t202 = g(1) * t253 - g(2) * t256;
t205 = g(1) * t256 + g(2) * t253;
t164 = t202 * t296 + t205 * t290;
t487 = t164 * t289;
t203 = g(1) * t254 - g(2) * t257;
t206 = g(1) * t257 + g(2) * t254;
t165 = t203 * t298 + t206 * t292;
t486 = t165 * t291;
t274 = t287 ^ 2;
t459 = t274 * t304;
t275 = t289 ^ 2;
t458 = t275 * t307;
t276 = t291 ^ 2;
t457 = t276 * t310;
t456 = t287 * t293;
t453 = t289 * t295;
t450 = t291 * t297;
t446 = t293 * t304;
t443 = t295 * t307;
t441 = t297 * t310;
t314 = pkin(1) ^ 2;
t435 = -t303 + t314;
t434 = -t306 + t314;
t433 = -t309 + t314;
t425 = (0.4e1 * t98 - 0.2e1 * t109) * t494;
t424 = t80 * t494;
t423 = (0.4e1 * t100 - 0.2e1 * t110) * t492;
t422 = t82 * t492;
t421 = (0.4e1 * t102 - 0.2e1 * t111) * t490;
t420 = t84 * t490;
t419 = t277 * t304 * t52;
t418 = t278 * t307 * t53;
t417 = t279 * t310 * t54;
t416 = t52 * t454;
t415 = t53 * t451;
t414 = t54 * t448;
t413 = t52 * t446;
t412 = t53 * t443;
t411 = t54 * t441;
t409 = t226 * t518;
t408 = t229 * t517;
t407 = t232 * t516;
t400 = t118 * t446;
t399 = t119 * t446;
t398 = t120 * t443;
t397 = t121 * t443;
t396 = t122 * t441;
t395 = t123 * t441;
t394 = t146 * t446;
t393 = t148 * t443;
t392 = t150 * t441;
t391 = t287 * t500;
t390 = t293 * t500;
t389 = t289 * t498;
t388 = t295 * t498;
t387 = t291 * t496;
t386 = t297 * t496;
t385 = t106 * t493;
t383 = t107 * t491;
t381 = t108 * t489;
t379 = t225 * t446;
t378 = t228 * t443;
t377 = t231 * t441;
t369 = -t314 + (-0.2e1 * pkin(1) - pkin(2)) * pkin(2);
t240 = -pkin(1) * t287 + t520;
t365 = t240 * t413;
t241 = -pkin(1) * t289 + t523;
t364 = t241 * t412;
t242 = -pkin(1) * t291 + t526;
t363 = t242 * t411;
t362 = t287 * t413;
t361 = t289 * t412;
t360 = t291 * t411;
t359 = t118 * t379;
t358 = t119 * t379;
t357 = t120 * t378;
t356 = t121 * t378;
t355 = t122 * t377;
t354 = t123 * t377;
t115 = t558 * t500;
t116 = t557 * t498;
t117 = t556 * t496;
t347 = t227 * t454 * t499;
t346 = t230 * t451 * t497;
t345 = t233 * t448 * t495;
t335 = -t286 * t509 - t285 * t506 + t284 * t468 + (t293 * t430 + t67 + (-t258 - t376) * t494) * t160 * t478 + (t112 * t439 + ((-t303 + t369) * t514 + t342 * t555) * t225) * t106 * t477 + t94 * t305 * t112;
t334 = -t286 * t508 - t285 * t505 + t284 * t467 + (t295 * t431 + t68 + (-t259 - t375) * t492) * t161 * t475 + (t113 * t438 + ((-t306 + t369) * t512 + t343 * t554) * t228) * t107 * t474 + t95 * t308 * t113;
t333 = -t286 * t507 - t285 * t504 + t284 * t466 + (t297 * t432 + t69 + (-t260 - t374) * t490) * t162 * t471 + (t114 * t437 + ((-t309 + t369) * t510 + t344 * t553) * t231) * t108 * t470 + t96 * t311 * t114;
t329 = t335 + t561;
t328 = t334 + t560;
t327 = t333 + t559;
t326 = -t335 + t544;
t325 = -t334 + t543;
t324 = -t333 + t542;
t168 = -t203 * t292 + t206 * t298;
t167 = -t202 * t290 + t205 * t296;
t166 = -t201 * t288 + t204 * t294;
t90 = (-t153 + t495 - t516) * t232;
t89 = (-t152 + t497 - t517) * t229;
t88 = (-t151 + t499 - t518) * t226;
t51 = 0.2e1 * t539;
t50 = 0.2e1 * t538;
t49 = 0.2e1 * t537;
t45 = t54 * t570 + 0.4e1 * t381;
t44 = t53 * t570 + 0.4e1 * t383;
t43 = t52 * t570 + 0.4e1 * t385;
t42 = t297 * t351 * t564 + t276 * t54;
t41 = t295 * t352 * t565 + t275 * t53;
t40 = t293 * t353 * t566 + t274 * t52;
t39 = t556 * t108 * t380 + t54 * t450;
t38 = t558 * t106 * t384 + t52 * t456;
t37 = t557 * t107 * t382 + t53 * t453;
t33 = t75 * t401 + (-t153 * t450 + t311 * t534) * t232 - t582;
t32 = t74 * t402 + (-t152 * t453 + t308 * t535) * t229 - t583;
t31 = t73 * t403 + (-t151 * t456 + t305 * t536) * t226 - t584;
t30 = t45 * t279 + ((t51 - t421) * t291 + t436) * t297 - 0.2e1 * t381;
t29 = t44 * t278 + ((t50 - t423) * t289 - t580) * t295 - 0.2e1 * t383;
t28 = t43 * t277 + ((t49 - t425) * t287 - t581) * t293 - 0.2e1 * t385;
t27 = -t291 * t407 + t297 * t36;
t26 = t291 * t36 + t297 * t407;
t25 = -t289 * t408 + t295 * t35;
t24 = t289 * t35 + t295 * t408;
t23 = -t287 * t409 + t293 * t34;
t22 = t287 * t34 + t293 * t409;
t21 = (t387 * t569 - t576) * t297 + 0.2e1 * t531 - t545 + t579;
t20 = (t389 * t568 - t575) * t295 + 0.2e1 * t532 - t546 + t578;
t19 = (t391 * t567 - t574) * t293 + 0.2e1 * t533 - t547 + t577;
t18 = (t297 * t45 + t436) * t291 + (t421 - 0.2e1 * t539) * t279 - 0.2e1 * t420 + t51;
t17 = (t295 * t44 - t580) * t289 + (t423 - 0.2e1 * t538) * t278 - 0.2e1 * t422 + t50;
t16 = (t293 * t43 - t581) * t287 + (t425 - 0.2e1 * t537) * t277 - 0.2e1 * t424 + t49;
t15 = (t386 * t570 + t576) * t291 + 0.2e1 * t559 - qJ(3,1) * t117 - t324;
t14 = (t388 * t570 + t575) * t289 + 0.2e1 * t560 - qJ(3,2) * t116 - t325;
t13 = (t390 * t570 + t574) * t287 + 0.2e1 * t561 - qJ(3,3) * t115 - t326;
t12 = -t585 - t559 + (-t105 * t310 + (-pkin(1) * t450 + qJ(3,1) * t279 - qJ(3,1)) * t153) * t232 + t324;
t11 = -t586 - t560 + (-t104 * t307 + (-pkin(1) * t453 + qJ(3,2) * t278 - qJ(3,2)) * t152) * t229 + t325;
t10 = -t587 - t561 + (-t103 * t304 + (-pkin(1) * t456 + qJ(3,3) * t277 - qJ(3,3)) * t151) * t226 + t326;
t9 = (0.4e1 * (t102 - t111 / 0.2e1) * qJ(3,1) * t490 + t433 * t54) * t279 + (0.2e1 * (qJ(3,1) * t381 + (-t420 + t539) * pkin(1)) * t291 + t436 * pkin(1)) * t297 + (((t548 / 0.2e1 - t266 / 0.2e1) * t257 - t198 / 0.2e1) * t291 - t539 / 0.2e1 + t420) * t569;
t8 = (0.4e1 * (t100 - t110 / 0.2e1) * qJ(3,2) * t492 + t434 * t53) * t278 + (0.2e1 * (qJ(3,2) * t383 + (-t422 + t538) * pkin(1)) * t289 - t580 * pkin(1)) * t295 + (((t540 / 0.2e1 - t262 / 0.2e1) * t256 - t194 / 0.2e1) * t289 - t538 / 0.2e1 + t422) * t568;
t7 = (0.4e1 * (t98 - t109 / 0.2e1) * qJ(3,3) * t494 + t435 * t52) * t277 + (0.2e1 * (qJ(3,3) * t385 + (-t424 + t537) * pkin(1)) * t287 - t581 * pkin(1)) * t293 + (((t541 / 0.2e1 - t261 / 0.2e1) * t255 - t190 / 0.2e1) * t287 - t537 / 0.2e1 + t424) * t567;
t6 = (t433 * t387 - t562) * t297 + t36 * t309 + (t327 + t585) * pkin(1) + (-t150 + t579) * qJ(3,1);
t5 = (t434 * t389 - t562) * t295 + t35 * t306 + (t328 + t586) * pkin(1) + (-t148 + t578) * qJ(3,2);
t4 = (t435 * t391 - t562) * t293 + t34 * t303 + (t329 + t587) * pkin(1) + (-t146 + t577) * qJ(3,3);
t3 = (-pkin(1) * t407 + t531 + t87) * t297 + (-t105 * t469 - t327) * t291 - t576;
t2 = (-pkin(1) * t408 + t532 + t86) * t295 + (-t104 * t473 - t328) * t289 - t575;
t1 = (-pkin(1) * t409 + t533 + t85) * t293 + (-t103 * t476 - t329) * t287 - t574;
t46 = [-t54 * t480 - t53 * t482 - t52 * t484, -t163 * t484 - t164 * t482 - t165 * t480, -t166 * t484 - t167 * t482 - t168 * t480, -t118 * t347 - t120 * t346 - t122 * t345 - t40 * t484 - t41 * t482 - t42 * t480, -t115 * t359 - t116 * t357 - t117 * t355 - 0.2e1 * t37 * t482 - 0.2e1 * t38 * t484 - 0.2e1 * t39 * t480, (t122 * t360 - t185 * t26) * t231 + (t120 * t361 - t183 * t24) * t228 + (t118 * t362 - t181 * t22) * t225, (t122 * t417 - t185 * t27) * t231 + (t120 * t418 - t183 * t25) * t228 + (t118 * t419 - t181 * t23) * t225, t34 * t359 + t35 * t357 + t355 * t36, (t122 * t501 - t165 * t185) * t472 + (t120 * t502 + t183 * t580) * t445 + (t118 * t503 + t181 * t581) * t479, (t122 * t392 + t185 * t486) * t231 + (t120 * t393 + t183 * t487) * t228 + (t118 * t394 + t181 * t488) * t225, t31 * t509 + t32 * t508 + t33 * t507 + (t15 * t396 - t185 * t30) * t231 + (t14 * t398 - t183 * t29) * t228 + (t13 * t400 - t181 * t28) * t225, t133 * t416 + t134 * t415 + t135 * t414 + (t122 * t363 - t185 * t3) * t231 + (t120 * t364 - t183 * t2) * t228 + (-t1 * t181 + t118 * t365) * t225, t88 * t509 + t89 * t508 + t90 * t507 + (-t18 * t185 + t21 * t396) * t231 + (-t17 * t183 + t20 * t398) * t228 + (-t16 * t181 + t19 * t400) * t225, t10 * t509 + t11 * t508 + t12 * t507 + (-t185 * t9 + t396 * t6) * t231 + (-t183 * t8 + t398 * t5) * t228 + (-t181 * t7 + t4 * t400) * t225, t286 - g(1); t54 * t481 + t53 * t483 + t52 * t485, t163 * t485 + t164 * t483 + t165 * t481, t166 * t485 + t167 * t483 + t168 * t481, -t119 * t347 - t121 * t346 - t123 * t345 + t40 * t485 + t41 * t483 + t42 * t481, -t115 * t358 - t116 * t356 - t117 * t354 + 0.2e1 * t37 * t483 + 0.2e1 * t38 * t485 + 0.2e1 * t39 * t481, (t123 * t360 + t184 * t26) * t231 + (t121 * t361 + t182 * t24) * t228 + (t119 * t362 + t180 * t22) * t225, (t123 * t417 + t184 * t27) * t231 + (t121 * t418 + t182 * t25) * t228 + (t119 * t419 + t180 * t23) * t225, t34 * t358 + t35 * t356 + t354 * t36, (t123 * t501 + t165 * t184) * t472 + (t121 * t502 - t182 * t580) * t445 + (t119 * t503 - t180 * t581) * t479, (t123 * t392 - t184 * t486) * t231 + (t121 * t393 - t182 * t487) * t228 + (t119 * t394 - t180 * t488) * t225, t31 * t506 + t32 * t505 + t33 * t504 + (t15 * t395 + t184 * t30) * t231 + (t14 * t397 + t182 * t29) * t228 + (t13 * t399 + t180 * t28) * t225, t136 * t416 + t137 * t415 + t138 * t414 + (t123 * t363 + t184 * t3) * t231 + (t121 * t364 + t182 * t2) * t228 + (t1 * t180 + t119 * t365) * t225, t88 * t506 + t89 * t505 + t90 * t504 + (t18 * t184 + t21 * t395) * t231 + (t17 * t182 + t20 * t397) * t228 + (t16 * t180 + t19 * t399) * t225, t10 * t506 + t11 * t505 + t12 * t504 + (t184 * t9 + t395 * t6) * t231 + (t182 * t8 + t397 * t5) * t228 + (t180 * t7 + t399 * t4) * t225, t285 - g(2); 0, 0, 0, -t386 * t457 - t388 * t458 - t390 * t459, -t115 * t454 - t116 * t451 - t117 * t448, t54 * t457 + t458 * t53 + t459 * t52, t360 + t361 + t362, t34 * t454 + t35 * t451 + t36 * t448, t145 * t454 + t147 * t451 + t149 * t448, t146 * t454 + t148 * t451 + t150 * t448, (t15 * t291 - t33 * t344) * t310 + (t14 * t289 - t32 * t343) * t307 + (t13 * t287 - t31 * t342) * t304, (-t344 + t242) * t414 + (-t343 + t241) * t415 + (-t342 + t240) * t416, (t21 * t291 - t344 * t90) * t310 + (t20 * t289 - t343 * t89) * t307 + (t19 * t287 - t342 * t88) * t304, (-t12 * t344 + t291 * t6) * t310 + (-t11 * t343 + t289 * t5) * t307 + (-t10 * t342 + t287 * t4) * t304, t284 - g(3);];
tauX_reg  = t46;
